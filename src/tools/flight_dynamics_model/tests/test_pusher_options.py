from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

from spearhead.aircraft import load_aircraft_config, pusher_options, selected_pusher_option
from spearhead.cg_card import cg_card_text
from spearhead.gui.flight_setup import FlightSetup
from spearhead.gui.flight_values import live_alerts
from spearhead.pusher import estimate_pusher
from spearhead.sitl.actuators import Channel
from spearhead.sitl.backend import CONFIG_ROOT, check_pusher_endpoints, load_setup

FDM_ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, FDM_ROOT / f"scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pt15_offers_the_hobbywing_pusher_by_default():
    assert set(pusher_options("pt15")) == {"hobbywing"}
    assert selected_pusher_option("pt15") == "hobbywing"
    system = load_aircraft_config("pt15").pusher_system
    assert system.option == "hobbywing"
    assert system.motor["model"] == "HOBBYWING PLATINUM 5230SL 380KV" and system.motor["kv_rpm_per_v"] == 380
    assert system.esc["model"] == "MAD AMPX 80A (5-14S)" and system.esc["pwm_range_us"] == [1050, 1940]
    assert system.propeller["diameter_in"] == 20 and system.propeller["pitch_in"] == 12
    assert system.output_limit == pytest.approx(400/890, abs=1e-4)
    with pytest.raises(ValueError, match="Unknown pusher option"):
        load_aircraft_config("pt15", pusher="mad")
    assert pusher_options("pt1_as_built") == {} and selected_pusher_option("pt1_as_built") is None


def test_hobbywing_fit_matches_the_profile_and_the_sheet():
    fit = _script("fit_hobbywing_pusher_motor")
    data, points = fit.load_points()
    used = points["current_a"] <= fit.MAX_FIT_CURRENT_A
    result = fit.fit(data, {key: value[used] for key, value in points.items()})
    assert used.sum() == 53
    assert result["duty_error_max"] < .05
    assert result["current_rms"] < .03 and result["current_max"] < .05
    stored = json.loads((FDM_ROOT / "data/propulsion/hobbywing_5230sl_380_pusher_fit.json").read_text())
    assert stored["equivalent_resistance_ohm"] == pytest.approx(result["resistance"])
    assert stored["source_sha256"] == hashlib.sha256(fit.SOURCE.read_bytes()).hexdigest()
    system = load_aircraft_config("pt15").pusher_system
    assert system.equivalent_resistance_ohm == pytest.approx(result["resistance"], rel=1e-4)
    assert system.motor["no_load_current_a"] == pytest.approx(result["idle"], rel=1e-4)


def test_kde_fit_script_still_reproduces_its_stored_fit():
    # Kept for note 0003 traceability; the KDE pusher is no longer a PT1.5 option.
    fit = _script("fit_kde_pusher_motor")
    points = fit.load_points()
    result = fit.fit(points)
    check = fit.check_against_test(points, result)
    assert check["rpm_rms"] < .01 and check["current_rms"] < .07
    stored = json.loads((FDM_ROOT / "data/propulsion/kde5215xf_220_pusher_fit.json").read_text())
    assert stored["equivalent_resistance_ohm"] == pytest.approx(result["resistance"])


def test_output_limit_scales_command_to_duty_and_matches_the_servo_endpoint():
    system = load_aircraft_config("pt15").pusher_system
    assert system.duty(0, 44.4) == 0
    assert system.duty(1, 44.4) == pytest.approx(system.output_limit)
    assert system.duty(.5, 44.4) == pytest.approx(.5*system.output_limit)
    assert system.duty(2, 44.4) == system.duty(1, 44.4)
    _, mapping, _ = load_setup(CONFIG_ROOT / "quad_transition.yaml")
    assert (mapping.pusher.minimum, mapping.pusher.maximum) == (1050, 1450)
    check_pusher_endpoints(mapping.pusher, system)
    with pytest.raises(ValueError, match="Pusher channel must span 1050-1450"):
        check_pusher_endpoints(Channel(3, minimum=1050, maximum=1940), system)
    unlimited = replace(system, output_limit=1.)
    check_pusher_endpoints(Channel(3, minimum=1050, maximum=1940), unlimited)
    with pytest.raises(ValueError, match="output limit"):
        replace(system, output_limit=0.)


def test_limited_hobbywing_pusher_stays_inside_esc_and_propeller_limits():
    params = load_aircraft_config("pt15")
    system = params.pusher_system
    static = estimate_pusher(system, 44.4, 1, 0, params.rho)
    assert 70 < static["thrust_n"] < 80 and static["effective_duty"] == pytest.approx(system.output_limit)
    assert params.Tmax == pytest.approx(static["thrust_n"])
    for voltage in (42., 44.4, 50.4):
        for speed in range(0, 30, 2):
            sample = estimate_pusher(system, voltage, 1, speed, params.rho)
            assert not sample["limited"] and not sample["above_sustainable_current"]
            assert sample["current_a"] < 60 and sample["rpm"] < system.propeller["max_rpm"]
    cruise = estimate_pusher(system, 44.4, 1, 26, params.rho)
    assert cruise["thrust_n"] > 40
    # Without the limit the 380KV motor would far exceed the ESC on 12S.
    unlimited = estimate_pusher(replace(system, output_limit=1.), 44.4, 1, 0, params.rho)
    assert unlimited["limited"] and unlimited["current_a"] == pytest.approx(system.max_current_a)
    assert system.sustainable_current_a == 80 and system.current_limit_label == "AMPX continuous rating (needs cooling)"


def test_flight_setup_and_sitl_config_carry_the_pusher_choice(tmp_path):
    assert FlightSetup().pusher == "hobbywing"
    data = FlightSetup(initial_point="hover20", maneuver="hold", pusher="hobbywing").config()
    assert data["pusher"] == "hobbywing"
    with pytest.raises(ValueError, match="pusher"):
        FlightSetup(pusher="kde").validate()
    path = tmp_path / "setup.yaml"
    path.write_text(yaml.safe_dump(data))
    physics, _, manifest = load_setup(path)
    assert physics.params.pusher_system.option == "hobbywing"
    default = yaml.safe_load((CONFIG_ROOT / "quad_hover.yaml").read_text())
    assert "pusher" not in default


def test_alert_and_cg_summary_name_the_pusher_limit():
    alerts = live_alerts({"pusher": {"above_sustainable_current": True, "sustainable_current_a": 44,
                                     "current_limit_label": "180 s motor rating"}})
    assert "Pusher above 44 A 180 s motor rating" in alerts
    card = {"static_margin_22": .064, "hover_max_rotor_fraction": .65,
            "pusher": {"option": "kde", "label": "KDE", "current_limit_a": 44, "current_limit_label": "180 s motor rating"},
            "trims": {20.: None, 22.: {"elevator_deg": -18.2, "pusher_command": .69, "beyond_measured": True,
                                        "pusher_current_a": 18.7, "above_current_limit": False}, 24.: None}}
    assert "Pusher at 22 m/s: 69 % · 18.7 A" in cg_card_text(card)
