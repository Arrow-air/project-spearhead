"""Protocol/physics checks independent of firmware; real smoke is a separate CLI."""

from dataclasses import replace
import json
import socket
import struct
import threading

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from spearhead.aircraft import load_aircraft_config
from spearhead.controls import Control
from spearhead.dynamics import dynamics
from spearhead.propulsion import propulsion_model
from spearhead.sitl.actuators import ActuatorMap, Channel
from spearhead.sitl.backend import Session, load_setup, serve, verify_firmware
from spearhead.sitl.physics import Physics, sensor_state
from spearhead.sitl.protocol import ServoFrame, encode_sensor


def packet(count=0, rate=200, channels=16):
    return struct.pack(f"<HHI{channels}H", 18458 if channels == 16 else 29569,
                       rate, count, *([1500] * channels))


@pytest.fixture(scope="module")
def setup():
    return load_setup()


def session(setup):
    physics, mapping, _ = setup
    return Session(Physics(physics.params, physics.state), mapping)


@pytest.mark.parametrize("channels", [16, 32])
def test_packet_formats(channels):
    frame = ServoFrame.decode(packet(23, 200, channels))
    assert (frame.rate, frame.count, len(frame.pwm)) == (200, 23, channels)


@pytest.mark.parametrize("bad", [b"", packet()[:-1], b"\0\0" + packet()[2:], packet(rate=0)])
def test_malformed_packets_rejected_before_physics(bad, setup):
    s = session(setup)
    with pytest.raises(ValueError):
        s.process(bad, ("127.0.0.1", 1234))
    assert s.physics.time == 0


def test_duplicate_wrap_rate_change_and_restart_policy(setup):
    s = session(setup)
    peer = ("127.0.0.1", 1234)
    first = s.process(packet(0xffffffff), peer)
    assert s.process(packet(0xffffffff), peer) == first
    assert s.physics.time == .005
    s.process(packet(0, rate=100), peer)
    assert s.physics.time == .01
    for bad, address in [(packet(0, rate=200), peer), (packet(5), peer),
                         (packet(0xffffffff), peer), (packet(1), ("127.0.0.1", 4321))]:
        with pytest.raises(ValueError):
            s.process(bad, address)
        assert s.physics.time == .01
    s.process(packet(1), peer)
    with pytest.raises(ValueError, match="restart"):
        s.process(packet(0), peer)
    assert s.physics.time == .015


def test_json_framing_and_nonfinite_rejection():
    assert encode_sensor({"x": 1}) == b'\n{"x":1}\n'
    for invalid in [float("nan"), float("inf")]:
        with pytest.raises(ValueError):
            encode_sensor({"imu": {"gyro": [invalid]}})


def test_surface_endpoints_reverse_and_motor_off():
    c = Channel(1)
    assert np.rad2deg(c.surface([1000])) == -5
    assert c.surface([0]) == c.surface([1500]) == 0
    assert np.rad2deg(replace(c, reverse=True).surface([2000])) == -5
    motor = Channel(1, minimum=1100, maximum=1940)
    assert [motor.motor([p]) for p in [0, 1100, 1520, 1940, 2000]] == [0, 0, .5, 1, 1]
    with pytest.raises(ValueError):
        ActuatorMap(pusher=Channel(1))
    with pytest.raises(ValueError):
        Channel(33)


@pytest.mark.parametrize("channel,rotor", [(5, 0), (8, 1), (6, 2), (7, 3)])
def test_quad_x_independent_rotor_wrench(channel, rotor):
    params = load_aircraft_config("pt15", aerodynamics="tail_update_2026_10_07")
    mapping = ActuatorMap(lift=tuple(Channel(i, minimum=1100, maximum=1940) for i in [5, 8, 6, 7]))
    pwm = [0] * 16
    pwm[channel-1] = 1520
    control = mapping.control(pwm)
    assert control.lift_motor_throttles[rotor] == .5
    assert sum(control.lift_motor_throttles) == .5
    force, moment = propulsion_model(control, params)
    arm = np.asarray(params.lift_system.positions_body_m[rotor]) - params.cg_body_xyz
    np.testing.assert_allclose(moment[:2], np.cross(arm, force)[:2])
    assert np.sign(moment[2]) == params.lift_system.reaction_torque_signs[rotor]


def test_freefall_specific_force_and_heading():
    params = replace(load_aircraft_config("reference_design", aerodynamics="tail_update_2026_10_07"),
                     use_aero_database=False, rho=0)
    state = np.zeros(12)
    state[3], state[8], state[9] = 10, np.pi/2, .1
    sensor = sensor_state(0, state, params, Control(0, 0, 0, 0))
    np.testing.assert_allclose(sensor["imu"]["accel_body"], 0, atol=1e-12)
    np.testing.assert_allclose(sensor["velocity"], [0, 10, 0], atol=1e-12)
    assert sensor["imu"]["gyro"] == [.1, 0, 0]


def test_hold_kinematics_and_trim(setup):
    physics, mapping, manifest = setup
    assert manifest["config"]["aerodynamics"] == "tail_update_2026_10_07"
    control = Control(**manifest["trim_control"])
    dx = dynamics(0, physics.state, physics.params, control)
    np.testing.assert_allclose(dx[3:12], 0, atol=1e-8)
    p = Physics(physics.params, physics.state)
    before = p.state.copy()
    sensor = p.advance(.1, control)
    np.testing.assert_allclose(p.state[:3] - before[:3], np.array(sensor["velocity"])*.1, atol=1e-12)
    np.testing.assert_allclose(p.state[3:], before[3:])
    assert np.linalg.norm(sensor["imu"]["accel_body"]) == pytest.approx(p.params.g)


def test_stepped_shared_ode_matches_offline_for_ten_seconds(setup):
    physics, _, manifest = setup
    control = Control(**manifest["trim_control"])
    state = physics.state.copy()
    state[7] += np.deg2rad(.05)
    offline = solve_ivp(lambda t, x: dynamics(t, x, physics.params, control), (0, 10), state,
                        rtol=1e-10, atol=1e-12, max_step=.05)
    coarse = Physics(physics.params, state, hold=False)
    fine = Physics(physics.params, state, hold=False, max_step=.00125)
    for _ in range(2000):
        coarse.advance(.005, control)
        fine.advance(.005, control)
    error = coarse.state - offline.y[:, -1]
    assert np.linalg.norm(error[:3]) < .01
    assert np.linalg.norm(error[3:6]) < .001
    assert max(abs(np.rad2deg(error[6:9]))) < .01
    assert np.linalg.norm(fine.state-offline.y[:, -1]) <= np.linalg.norm(error) * 1.1 + 1e-9


def test_real_loopback_udp_round_trip_and_timeout(setup, tmp_path):
    s = session(setup)
    stop, errors = threading.Event(), []
    with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as backend_socket:
        backend_socket.bind(("127.0.0.1", 0))

        def work():
            try:
                serve(backend_socket, s, tmp_path / "truth.jsonl", stop, timeout_s=.25)
            except TimeoutError as exc:
                errors.append(str(exc))

        worker = threading.Thread(target=work)
        worker.start()
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as client:
                client.settimeout(2)
                client.sendto(packet(), backend_socket.getsockname())
                response, _ = client.recvfrom(4096)
                assert json.loads(response)["timestamp"] == .005
                client.sendto(packet(), backend_socket.getsockname())
                assert client.recvfrom(4096)[0] == response
            worker.join(timeout=2)
            assert not worker.is_alive() and len(errors) == 1
        finally:
            stop.set()
            worker.join(timeout=2)


@pytest.mark.parametrize("identity", ["ArduPlane V4.8.0-dev (a4d4ec06)",
                                      "ArduPlane V4.7.0 (1511f271)", "ArduPlane V4.7.1 (00000000)"])
def test_wrong_firmware_rejected(tmp_path, identity):
    binary = tmp_path / "arduplane"
    binary.write_bytes(identity.encode() + b"\0")
    with pytest.raises(ValueError):
        verify_firmware(binary)


def test_pinned_firmware_identity(tmp_path):
    binary = tmp_path / "arduplane"
    binary.write_bytes(b"ArduPlane V4.7.1 (dbe79216)\0ArduPlane V4.7.1\0")
    assert verify_firmware(binary)["tag"] == "Plane-4.7.1"


def test_optional_mavlink_nan_preserved_as_unavailable():
    from spearhead.sitl.__main__ import telemetry_record
    record = telemetry_record({"optional": float("nan"), "valid": 3.5})
    assert record["message"] == {"optional": None, "valid": 3.5}
    assert record["unavailable_fields"] == ["message.optional"]


def test_release_at_new_frame_boundary_uses_actual_pwm(setup):
    s = session(setup)
    peer = ("127.0.0.1", 1234)
    first = s.process(packet(), peer)
    s.release()
    assert s.process(packet(), peer) == first
    assert s.physics.hold
    before = s.physics.state.copy()
    s.process(packet(1), peer)
    assert not s.physics.hold
    assert s.last_record["control"]["throttle"] == .5
    assert s.last_record["control"]["de"] == 0
    assert not np.array_equal(s.physics.state[3:], before[3:])


@pytest.mark.parametrize("key,value", [("aerodynamics", "adb_v1_1"), ("firmware", "Plane-4.7.0"),
                                       ("cg_body_xyz_m", [float("nan"), 0, 0]),
                                       ("physics_rate_hz", 0), ("max_step_s", .1)])
def test_bad_setup_rejected(tmp_path, key, value):
    import yaml
    from spearhead.sitl.backend import CONFIG_ROOT
    config = yaml.safe_load((CONFIG_ROOT / "fixed_wing_new.yaml").read_text())
    config[key] = value
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError):
        load_setup(path)
