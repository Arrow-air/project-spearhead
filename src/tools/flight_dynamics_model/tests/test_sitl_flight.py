"""Independent direction, bootstrap boundary, and acceptance metric checks."""

from dataclasses import replace
import json

import numpy as np
import pytest

from spearhead.controls import Control
from spearhead.dynamics import dynamics
from spearhead.sitl.backend import CONFIG_ROOT, Session, load_setup
from spearhead.sitl.flight import bootstrap_mode, flight_parameters, manual_elevator_rc, parameters_match, startup_states, surface_pwm
from spearhead.sitl.metrics import analyse_truth, analyse_elevator_truth, command_at, elevator_at, step_response
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
    assert abs(c.surface(pwm)-getattr(trim, axis)) <= np.deg2rad(c.travel_deg/1000+.0001)
    plus = pwm.copy()
    plus[c.number-1] += 20
    delta = c.surface(plus)-c.surface(pwm)
    modified = replace(trim, **{axis: getattr(trim, axis)+delta})
    acceleration = dynamics(0, physics.state, physics.params, modified)[9:12]
    assert acceleration[{"da": 0, "de": 1, "dr": 2}[axis]] > 0
    params = flight_parameters(mapping, physics.state, trim)
    assert params["ARMING_SKIPCHK"] == 0 and params["ARMING_REQUIRE"] == 1
    assert params["AHRS_EKF_TYPE"] == 3
    assert params["SERVO2_TRIM"] == 1575


def test_quantized_roll_command_targets_match_rc_resolution():
    assert command_at(83) == (0, 8, .48)
    assert command_at(101) == (0, -8, -.48)
    assert command_at(10) == (1, 10, .2)
    assert command_at(25) == (None, 0, 0.)


def test_manual_elevator_input_accounts_for_asymmetric_servo_trim():
    physics, mapping, manifest = load_setup(CONFIG_ROOT / "elevator_demo.yaml")
    trim = Control(**manifest["trim_control"])
    midpoint = surface_pwm(mapping.elevator, trim.de)
    assert bootstrap_mode("elevator", 25, 50) == 0
    assert elevator_at(9.99) == 0 and elevator_at(10) == 1
    assert elevator_at(11.5) == -1 and elevator_at(13) == 0
    for delta in (0, 1, -1):
        rc = manual_elevator_rc(mapping, trim, midpoint, delta)
        # Official SRV_Channel's MANUAL angle mapping uses the relevant
        # servo endpoint span around SERVO2_TRIM, not an equal +/-500 PWM.
        normalized = (rc-1500)/500
        span = 2000-midpoint if normalized >= 0 else midpoint-1000
        pwm = midpoint+int(normalized*span)
        actual = mapping.elevator.surface([1500, pwm, 1300, 1500])
        # One output PWM count is 0.04 deg with the provisional +/-20 travel.
        assert np.rad2deg(actual-trim.de) == pytest.approx(delta, abs=mapping.elevator.travel_deg/500)
    with pytest.raises(ValueError, match="travel"):
        manual_elevator_rc(mapping, trim, midpoint, -21)


def test_elevator_demo_requires_measured_pulses_and_recovery(tmp_path):
    initial = np.array([0, 0, -100, 32, 0, -.5, 0, np.deg2rad(-1), 0, 0, 0, 0])
    trim = Control(de=np.deg2rad(-3), da=0, dr=0, throttle=.3)
    rows = []
    for t in np.arange(0, 45.01, .01):
        state = initial.copy()
        # Test the evaluator with a visible excursion that returns by 20 s.
        state[7] -= np.deg2rad(4*np.exp(-((t-12)/2)**2))
        row = {"hold": False, "time_s": float(t), "state": state.tolist(),
               "pwm": [1500, 1800, 1300, 1500],
               "control": {"de": trim.de+np.deg2rad(elevator_at(t))},
               "alpha_deg": -1., "beta_deg": 0., "control_coverage_exceeded": False}
        rows.append(row)
    path = tmp_path / "truth.jsonl"
    def save():
        path.write_text("\n".join(json.dumps(row) for row in rows))
    save()
    report = analyse_elevator_truth(path, initial, trim, 45)
    assert report["passed"] and report["pitch_excursion_deg"] == pytest.approx(4)
    assert "attitude_rms_error_deg" not in report
    # A visible animation without the actual surface command must fail.
    for row in rows:
        row["control"]["de"] = trim.de
    save()
    assert not analyse_elevator_truth(path, initial, trim, 45)["passed"]
    for row in rows:
        row["control"]["de"] = trim.de+np.deg2rad(elevator_at(row["time_s"]))
        if row["time_s"] > 40:
            row["state"][7] = initial[7]+np.deg2rad(1)
    save()
    assert not analyse_elevator_truth(path, initial, trim, 45)["passed"]


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


def test_provisional_20_degree_pwm_mapping_and_control_gain_scaling(setup):
    physics, mapping, manifest = setup
    assert manifest['control_model']['travel_calibrated'] is False
    assert manifest['control_model']['deflection_policy'] == 'linear'
    for name in ('elevator', 'aileron', 'rudder'):
        channel = getattr(mapping, name)
        assert channel.travel_deg == 20
    elevator = mapping.elevator
    for pwm, angle in [(1000, 20), (1500, 0), (1575, -3), (1800, -12), (2000, -20)]:
        assert np.rad2deg(elevator.surface([1500, pwm, 1300, 1500])) == pytest.approx(angle)
    trim = Control(**manifest['trim_control'])
    gains = flight_parameters(mapping, physics.state, trim)
    assert gains['PTCH_RATE_FF'] == .45
    assert gains['PTCH_RATE_P'] == .08
    assert gains['PTCH_RATE_I'] == .15
    assert gains['PTCH2SRV_TCONST'] == .2
    assert gains['RLL_RATE_FF'] == pytest.approx(.08625)
