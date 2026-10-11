import json

import numpy as np
import pytest
import yaml

from spearhead.controls import Control
from spearhead.dynamics import dynamics
from spearhead.sitl.backend import CONFIG_ROOT, load_setup
from spearhead.sitl.quad import hover_thrust_setting, quad_parameters, analyse_hover, PARAMETER_RENAMES
from spearhead.gui.sitl_page import sitl_figures


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
    path = tmp_path / "shifted.yaml"; path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="does not balance"):
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
    home = {"latitude": 473977432, "longitude": 85455940, "altitude": 500000}
    assert ground_home_matches(home)
    assert not ground_home_matches({**home, "altitude": 520000})
    assert not ground_home_matches({**home, "latitude": 483977420})
    assert not ground_home_matches({})
