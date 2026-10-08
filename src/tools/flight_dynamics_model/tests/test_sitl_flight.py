"""Independent direction, bootstrap boundary, and acceptance metric checks."""

from dataclasses import replace
import json

import numpy as np
import pytest

from spearhead.controls import Control
from spearhead.dynamics import dynamics
from spearhead.sitl.backend import CONFIG_ROOT, Session, load_setup
from spearhead.sitl.flight import bootstrap_mode, flight_parameters, parameters_match, startup_states, surface_pwm
from spearhead.sitl.metrics import analyse_truth, command_at, step_response
from spearhead.sitl.physics import Physics
from test_sitl import packet


@pytest.fixture(scope="module")
def setup():
    return load_setup(CONFIG_ROOT / "fixed_wing_control.yaml")


def test_calibration_to_air_start_and_release_are_distinct_frame_events(setup):
    physics, mapping, _ = setup
    s = Session(Physics(physics.params, np.zeros(12)), mapping)
    peer = ("127.0.0.1", 1234)
    first = s.process(packet(), peer)
    disturbed = physics.state.copy()
    disturbed[6:8] += np.deg2rad([.5, -.5])
    s.air_start(disturbed)
    assert s.process(packet(), peer) == first
    np.testing.assert_equal(s.physics.state, 0)
    s.process(packet(1), peer)
    assert s.last_record["bootstrap_air_start"] and s.physics.hold
    from spearhead.rotations import euler_body_to_ned
    expected_position = disturbed[:3] + euler_body_to_ned(*disturbed[6:9]) @ disturbed[3:6] * .005
    np.testing.assert_allclose(s.physics.state[:3], expected_position, atol=1e-12)
    before = s.physics.state.copy()
    np.testing.assert_allclose(s.physics.state[6:8], disturbed[6:8])
    s.release()
    cached = s.process(packet(1), peer)
    np.testing.assert_equal(s.physics.state, before)
    s.process(packet(2), peer)
    assert not s.physics.hold and s.last_record["release_applied"]
    with pytest.raises(ValueError):
        s.air_start(physics.state)
    assert cached.startswith(b"\n")


@pytest.mark.parametrize("axis,channel", [("da", "aileron"), ("de", "elevator"), ("dr", "rudder")])
def test_servo_trim_quantization_and_firmware_positive_command_moment(setup, axis, channel):
    physics, mapping, manifest = setup
    trim = Control(**manifest["trim_control"])
    c = getattr(mapping, channel)
    pwm = [1500] * 16
    pwm[c.number-1] = surface_pwm(c, getattr(trim, axis))
    assert abs(c.surface(pwm)-getattr(trim, axis)) <= np.deg2rad(.0051)
    plus = pwm.copy()
    plus[c.number-1] += 20
    delta = c.surface(plus)-c.surface(pwm)
    modified = replace(trim, **{axis: getattr(trim, axis)+delta})
    acceleration = dynamics(0, physics.state, physics.params, modified)[9:12]
    assert acceleration[{"da": 0, "de": 1, "dr": 2}[axis]] > 0
    params = flight_parameters(mapping, physics.state, trim)
    assert params["ARMING_SKIPCHK"] == 0 and params["ARMING_REQUIRE"] == 1
    assert params["AHRS_EKF_TYPE"] == 3
    assert params["SERVO2_TRIM"] == 1801


def test_quantized_roll_command_targets_match_rc_resolution():
    assert command_at(83) == (0, 8, .48)
    assert command_at(101) == (0, -8, -.48)
    assert command_at(10) == (1, 10, .2)
    assert command_at(25) == (None, 0, 0.)


def test_initial_attitude_is_established_before_sensor_alignment(setup):
    physics, _, _ = setup
    stationary, airborne = startup_states(physics.state, .5, -.5)
    np.testing.assert_equal(stationary[3:6], 0)
    np.testing.assert_equal(stationary[6:9], airborne[6:9])
    np.testing.assert_allclose(airborne[6:8]-physics.state[6:8], np.deg2rad([.5, -.5]))
    assert bootstrap_mode("fbwa", None, 24) == 0
    assert bootstrap_mode("fbwa", 25, 27.9) == 0
    assert bootstrap_mode("fbwa", 25, 28) == 5
    assert bootstrap_mode("manual", 25, 40) == 0
    with pytest.raises(ValueError):
        startup_states(physics.state, .6, 0)


def test_readback_accepts_mavlink_float_resolution_but_not_missing_or_skipped_checks():
    expected = {"ARMING_SKIPCHK": 0., "PTCH_RATE_P": .08}
    assert parameters_match({"ARMING_SKIPCHK": 0., "PTCH_RATE_P": .079999998}, expected)
    assert not parameters_match({"ARMING_SKIPCHK": 1., "PTCH_RATE_P": .08}, expected)
    assert not parameters_match({"PTCH_RATE_P": .08}, expected)


def test_step_metrics_detect_slow_unsettled_and_overshooting_responses():
    t = np.linspace(0, 12, 1201)
    fast = step_response(t, .5*(1-np.exp(-t)), .5, .5)
    assert fast["passed"] and fast["settling_s"] == pytest.approx(1.61)
    assert not step_response(t, .5*(1-np.exp(-t/10)), .5, .5)["passed"]
    assert not step_response(t, .5+np.exp(-t)*.3, .5, .5)["passed"]
    assert not step_response(t, np.zeros_like(t), .5, .5)["passed"]
    assert step_response(t, -.5*(1-np.exp(-t)), -.5, -.5)["passed"]


def test_nominal_metrics_exclude_hold_and_reject_sustained_saturation(tmp_path):
    state = [0, 0, -100, 28, 0, 0, 0, 0, 0, 0, 0, 0]
    rows = [{"hold": True, "state": [100]*12}]
    for i in range(401):
        rows.append({"hold": False, "time_s": i*.005, "state": state,
                     "pwm": [1500, 1801, 1254, 1495], "alpha_deg": 0, "beta_deg": 0,
                     "control_coverage_exceeded": False})
    path = tmp_path / "truth.jsonl"
    path.write_text("\n".join(json.dumps(row) for row in rows))
    report = analyse_truth(path, np.array(state), "fbwa", 2)
    assert report["passed"] and report["physics_frames"] == 401
    for row in rows[1:]:
        row["pwm"][0] = 2000
    path.write_text("\n".join(json.dumps(row) for row in rows))
    assert not analyse_truth(path, np.array(state), "fbwa", 2)["passed"]
    # Bootstrap may move the aircraft before release; free-flight altitude
    # metrics must use the release position, not count held climb as flight.
    for row in rows[1:]:
        row["state"] = [*state[:2], -105, *state[3:]]
        row["pwm"][0] = 1500
    path.write_text("\n".join(json.dumps(row) for row in rows))
    report = analyse_truth(path, np.array(state), "fbwa", 2)
    assert report["max_altitude_deviation_m"] == 0
    assert report["max_altitude_deviation_from_air_start_m"] == 5
