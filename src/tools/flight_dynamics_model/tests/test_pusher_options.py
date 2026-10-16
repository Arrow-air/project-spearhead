import importlib.util
import json
from pathlib import Path

import pytest
import yaml

from spearhead.aircraft import load_aircraft_config, pusher_options
from spearhead.cg_card import cg_card_text
from spearhead.gui.flight_setup import FlightSetup
from spearhead.gui.flight_values import live_alerts
from spearhead.pusher import estimate_pusher
from spearhead.sitl.backend import CONFIG_ROOT, load_setup

FDM_ROOT = Path(__file__).resolve().parents[1]


def _fit_module():
    spec = importlib.util.spec_from_file_location("fit_kde", FDM_ROOT / "scripts/fit_kde_pusher_motor.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pt15_offers_mad_and_kde_pushers_with_mad_default():
    assert set(pusher_options("pt15")) == {"mad", "kde"}
    assert load_aircraft_config("pt15").pusher_system.option == "mad"
    kde = load_aircraft_config("pt15", pusher="kde").pusher_system
    assert kde.motor["model"] == "KDE5215XF-220" and kde.esc["model"] == "KDE-UAS75HVC+"
    assert kde.propeller["diameter_in"] == 20 and kde.propeller["pitch_in"] == 12
    with pytest.raises(ValueError, match="Unknown pusher option"):
        load_aircraft_config("pt15", pusher="other")
    assert pusher_options("pt1_as_built") == {}


def test_kde_fit_matches_the_profile_and_reproduces_the_static_test():
    fit = _fit_module()
    points = fit.load_points()
    result = fit.fit(points)
    check = fit.check_against_test(points, result)
    assert check["rpm_rms"] < .01 and check["rpm_max"] < .01
    assert check["current_rms"] < .07 and check["current_max"] < .2
    stored = json.loads((FDM_ROOT / "data/propulsion/kde5215xf_220_pusher_fit.json").read_text())
    assert stored["equivalent_resistance_ohm"] == pytest.approx(result["resistance"])
    kde = load_aircraft_config("pt15", pusher="kde").pusher_system
    assert kde.equivalent_resistance_ohm == pytest.approx(result["resistance"], rel=1e-5)
    assert kde.motor["no_load_current_a"] == pytest.approx(result["idle"], rel=1e-3)


def test_kde_throttle_maps_to_duty_and_mad_uses_throttle_directly():
    kde = load_aircraft_config("pt15", pusher="kde").pusher_system
    mad = load_aircraft_config("pt15", pusher="mad").pusher_system
    assert kde.duty(0, 44.4) == 0
    duties = [kde.duty(t, 44.4) for t in (.1, .25, .5, .75, 1.)]
    assert all(b > a for a, b in zip(duties, duties[1:]))
    assert 0.79 < kde.duty(1, 44.4) < 0.82
    assert kde.duty(1, 60) == kde.duty(1, 50.4)
    assert mad.duty(.7, 44.4) == pytest.approx(.7)


def test_kde_pusher_performance_and_current_limit():
    kde_params = load_aircraft_config("pt15", pusher="kde")
    mad_params = load_aircraft_config("pt15", pusher="mad")
    kde, mad = kde_params.pusher_system, mad_params.pusher_system
    static_kde = estimate_pusher(kde, 44.4, 1, 0, kde_params.rho)
    static_mad = estimate_pusher(mad, 44.4, 1, 0, mad_params.rho)
    assert static_kde["thrust_n"] > static_mad["thrust_n"]
    assert 0.79 < static_kde["effective_duty"] < 0.82
    cruise = estimate_pusher(kde, 44.4, 1, 26, kde_params.rho)
    assert cruise["thrust_n"] > 35 and cruise["current_a"] < kde.sustainable_current_a
    assert kde.sustainable_current_a == 44 and kde.current_limit_label == "180 s motor rating"
    assert kde_params.Tmax == pytest.approx(static_kde["thrust_n"])


def test_flight_setup_and_sitl_config_carry_the_pusher_choice(tmp_path):
    data = FlightSetup(initial_point="hover20", maneuver="hold", pusher="kde").config()
    assert data["pusher"] == "kde"
    with pytest.raises(ValueError, match="pusher"):
        FlightSetup(pusher="other").validate()
    path = tmp_path / "setup.yaml"
    path.write_text(yaml.safe_dump(data))
    physics, _, manifest = load_setup(path)
    assert physics.params.pusher_system.option == "kde"
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
