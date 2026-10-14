import json

import numpy as np
import pytest
import yaml

from spearhead.controls import Control
from spearhead.dynamics import dynamics
from spearhead.sitl.backend import CONFIG_ROOT, load_setup
from spearhead.sitl.quad import hover_thrust_setting, quad_parameters, analyse_hover, PARAMETER_RENAMES
from spearhead.gui.sitl_page import sitl_figures


@pytest.mark.parametrize("config_name", ["quad_hover.yaml", "quad_ground.yaml", "quad_transition.yaml"])
def test_prototype_proposed_travel_reaches_physical_limits_without_changing_motor_mapping(config_name):
    physics, mapping, manifest = load_setup(CONFIG_ROOT / config_name)
    neutral = [1500] * 8
    neutral[2] = 1050
    for index, sign, surface in [(0, 1, "da"), (1, -1, "de"), (3, 1, "dr")]:
        for pwm, expected in [(800, -30), (1000, -30), (1500, 0), (2000, 30), (2200, 30)]:
            command = neutral.copy()
            command[index] = pwm
            control = mapping.control(command)
            assert np.rad2deg(getattr(control, surface)) == pytest.approx(sign * expected)
            assert control.throttle == 0
            assert control.lift_motor_throttles == mapping.control(neutral).lift_motor_throttles
    assert physics.params.control_surface_limits_deg == (30, 30, 30)
    assert manifest["control_model"]["travel_calibrated"] is False


def test_hover_setup_balances_actual_four_motor_wrench_with_pusher_off():
    physics, mapping, manifest = load_setup(CONFIG_ROOT / "quad_hover.yaml")
    control = Control(**manifest["trim_control"])
    assert physics.params.aircraft_id == "pt1_as_built" and physics.params.pusher_installed
    assert control.throttle == 0
    assert [channel.number for channel in mapping.lift] == [5, 8, 6, 7]
    assert physics.params.mass == 21.0
    assert control.lift_throttle == pytest.approx(.6663451997065988)
    np.testing.assert_allclose(dynamics(0, physics.state, physics.params, control), 0, atol=1e-12)
    assert manifest["config"]["flight_regime"] == "hover"


def test_unbalanced_hover_cg_is_rejected_instead_of_hidden_compensation(tmp_path):
    config = yaml.safe_load((CONFIG_ROOT / "quad_hover.yaml").read_text())
    config["cg_body_xyz_m"][0] += .1
    config["hover_trim"] = "collective"
    path = tmp_path / "shifted.yaml"; path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="does not balance"):
        load_setup(path)


@pytest.mark.parametrize("cg_x", [-.10, .025, .05])
def test_explicit_per_rotor_hover_trim_balances_shifted_cg(tmp_path, cg_x):
    config = yaml.safe_load((CONFIG_ROOT / "quad_hover.yaml").read_text())
    config["cg_body_xyz_m"][0] = cg_x
    config["hover_trim"] = "independent_rotors"
    path = tmp_path / "shifted.yaml"; path.write_text(yaml.safe_dump(config))
    physics, _, manifest = load_setup(path)
    control = Control(**manifest["trim_control"])
    assert control.lift_motor_throttles is not None
    assert np.ptp(control.lift_motor_throttles) > 0
    assert control.throttle == 0
    np.testing.assert_allclose(dynamics(0, physics.state, physics.params, control), 0, atol=1e-8)
    np.testing.assert_allclose(physics.params.cg_body_xyz, [cg_x, 0, .15])


def test_per_rotor_hover_trim_rejects_unreachable_cg(tmp_path):
    config = yaml.safe_load((CONFIG_ROOT / "quad_hover.yaml").read_text())
    config["cg_body_xyz_m"][0] = 2
    config["hover_trim"] = "independent_rotors"
    path = tmp_path / "unreachable.yaml"; path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="no balanced"):
        load_setup(path)


def test_hover_thrust_setting_inverts_pinned_firmware_expo_and_spin_range():
    actuator = .6663451997065988
    thrust = hover_thrust_setting(actuator)
    restored = .15 + .8*((.75-1)+np.sqrt((1-.75)**2+4*.75*thrust))/(2*.75)
    assert restored == pytest.approx(actuator)
    assert thrust == pytest.approx(.47379424044743007, abs=1e-7)
    with pytest.raises(ValueError, match="usable"):
        hover_thrust_setting(.99)


def test_quad_controller_retains_logged_gains_and_checks_enabled():
    baseline = CONFIG_ROOT / "quad_hover_controller.json"
    logged = json.loads(baseline.read_text())["parameters"]
    params = quad_parameters(.6663451997065988, baseline)
    assert all(params[PARAMETER_RENAMES.get(k, k)] == v for k,v in logged.items())
    assert "Q_P_JERK_NE" not in params and params["Q_P_NE_JERK"] == 2
    assert params["ARMING_REQUIRE"] == 1 and params["ARMING_SKIPCHK"] == 0
    assert params["SCHED_LOOP_RATE"] >= 100 and params["Q_ASSIST_SPEED"] == -1
    assert params["Q_ENABLE"] == 1
    assert [params[f"SERVO{i}_FUNCTION"] for i in range(5,9)] == [33,34,35,36]


def hover_records(path, altitude=20):
    rows = [{"time_s": i/10, "hold": False, "state": [0,0,-altitude,0,0,0,0,0,0,0,0,0],
             "pwm": [1500]*4+[1740]*4, "control": {"de": 0}}
            for i in range(201)]
    path.write_text(''.join(json.dumps(row)+'\n' for row in rows))
    return rows


def test_hover_acceptance_checks_altitude_and_actual_motor_channels(tmp_path):
    path = tmp_path / "truth.jsonl"
    rows = hover_records(path)
    report = analyse_hover(path, 20, 20)
    assert report["passed"] and report["final_lift_pwm_us"] == [1740]*4
    hover_records(path, 16)
    report = analyse_hover(path, 20, 20)
    assert not report["passed"] and not report["checks"]["altitude_hold"]
    figures = sitl_figures({"history": rows[:2], "manifest": {"config": {"flight_regime": "hover"}, "trim_control": {"de": 0}}})
    assert len(figures["PWM"].data) == 4
    assert figures["PWM"].data[0].name == "FR · M1"


def test_quad_manual_throttle_follows_selected_mode_not_bootstrap_target():
    from spearhead.sitl.quad import quad_throttle_pwm
    hover = .47379424044743007
    pwm = quad_throttle_pwm(17, hover)
    assert pwm == quad_throttle_pwm(23, hover) == 1505
    # Decode using the firmware RANGE mapping, including its low-end deadzone.
    assert (pwm-1060)/940 == pytest.approx(hover, abs=.5/940)
    for mode in (18, 19, 20, 21, 22, 25, 5):
        assert quad_throttle_pwm(mode, hover) == 1500


def test_ground_home_rejects_air_start_altitude_and_wrong_position():
    from spearhead.sitl.quad import ground_home_matches
    home = {"latitude": 398582000, "longitude": 325460000, "altitude": 1125000}
    assert ground_home_matches(home)
    assert not ground_home_matches({**home, "altitude": 1145000})
    assert not ground_home_matches({**home, "latitude": 408582000})
    assert not ground_home_matches({})


def test_lift_plot_uses_motor_endpoints_instead_of_servo_limits_or_startup_zero(tmp_path):
    rows = hover_records(tmp_path/'truth.jsonl')[:3]
    rows[0]['hold'] = True; rows[0]['pwm'][4:8] = [0]*4
    rows[1]['pwm'][4:8] = [1100, 1610, 1740, 1940]
    rows[2]['pwm'][4] = 2005
    manifest = {'config': {'flight_regime': 'hover'}, 'trim_control': {'de': 0},
                'parameter_overlay': {'Q_M_PWM_MIN': 1100, 'Q_M_PWM_MAX': 1940,
                                      'SERVO5_MIN': 1100, 'SERVO5_MAX': 1900}}
    figure = sitl_figures({'history': rows, 'manifest': manifest})['PWM']
    assert figure.layout.yaxis.range == (1100, 1940)
    assert figure.layout.yaxis.autorange is False
    assert figure.data[0].y == (None, 1100, 2005)  # preserve actual out-of-range values
    assert [shape.y0 for shape in figure.layout.shapes] == [1100, 1940]
    assert any('1 motor samples outside limits' in a.text for a in figure.layout.annotations)
    assert rows[0]['pwm'][4] == 0  # plotting never changes the physical recording


def test_lift_plot_reads_archived_config_limits_without_guessing_missing_metadata():
    from spearhead.gui.sitl_page import lift_pwm_limits
    assert lift_pwm_limits({'config': {'channels': {'lift': [{'minimum': 1100, 'maximum': 1940}]}}}) == (1100, 1940)
    assert lift_pwm_limits({}) is None
    assert lift_pwm_limits({'parameter_overlay': {'Q_M_PWM_MIN': 1940, 'Q_M_PWM_MAX': 1100}}) is None
