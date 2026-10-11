import json

import pytest
import yaml

from spearhead.gui.flight_setup import FlightSetup
from spearhead.gui import sitl_session
from spearhead.sitl.backend import load_setup


@pytest.mark.parametrize("setup,aircraft,speed,altitude,scenario", [
    (FlightSetup(), "pt1_as_built", 0, 20, "quad_interactive"),
    (FlightSetup(initial_point="custom_hover", altitude_m=32, duration_s=30, maneuver="hold"), "pt1_as_built", 0, 32, "quad_hover"),
    (FlightSetup("reference_design", "trim26", "steps", 118, 120, 26), "reference_design", 26, 120, "steps"),
    (FlightSetup("reference_design", "custom_trim", "elevator", 45, 100, 27), "reference_design", 27, 100, "elevator"),
])
def test_selection_becomes_authoritative_runner_configuration(setup, aircraft, speed, altitude, scenario):
    config = setup.config()
    assert config["aircraft"] == aircraft
    assert config["speed_m_s"] == speed and config["altitude_m"] == altitude
    assert config["aerodynamics"] == "tail_update_2026_10_07"
    assert setup.scenario == scenario


@pytest.mark.parametrize("changes", [
    {"aircraft": "unknown"}, {"initial_point": "trim28"}, {"maneuver": "elevator"},
    {"duration_s": float("nan")}, {"duration_s": None}, {"duration_s": 601},
    {"altitude_m": 0}, {"speed_m_s": 28}, {"altitude_m": 30},
    {"initial_point": "custom_hover", "altitude_m": 40, "maneuver": "hold"},
])
def test_invalid_or_unmodelled_choices_are_rejected(changes):
    with pytest.raises(ValueError):
        FlightSetup(**changes).validate()


def test_named_trim_and_maneuver_duration_are_not_silently_overridden():
    with pytest.raises(ValueError, match="118"):
        FlightSetup("reference_design", "trim28", "steps", 60, 100, 28).validate()
    with pytest.raises(ValueError, match="Custom forward"):
        FlightSetup("reference_design", "trim28", "hold", 60, 100, 26).validate()


def test_custom_hover_reaches_physics_and_manifest(tmp_path):
    config = FlightSetup(initial_point="custom_hover", altitude_m=32, duration_s=30).config()
    path = tmp_path / "setup.yaml"
    path.write_text(yaml.safe_dump(config))
    physics, _, manifest = load_setup(path)
    assert physics.state[2] == -32
    assert physics.params.mass == 21
    assert manifest["config"]["altitude_m"] == 32


def test_ground_recipe_starts_on_gear_and_preserves_takeoff_target(tmp_path):
    setup = FlightSetup(initial_point="ground", maneuver="takeoff_hover_land", altitude_m=5, duration_s=15)
    assert setup.scenario == "quad_ground"
    path = tmp_path / "ground.yaml"
    path.write_text(yaml.safe_dump(setup.config()))
    physics, _, manifest = load_setup(path)
    assert -.3 < physics.state[2] < -.29
    assert not physics.hold and physics.ground is not None
    assert manifest["config"]["altitude_m"] == 5
    assert manifest["ground_contact"]["provisional"]
    with pytest.raises(ValueError):
        FlightSetup(initial_point="ground").validate()
    with pytest.raises(ValueError):
        FlightSetup(maneuver="takeoff_hover_land").validate()


def test_manager_launches_selected_configuration_and_duration(tmp_path, monkeypatch):
    calls = []
    class Process:
        def poll(self):
            return None
    monkeypatch.setattr(sitl_session, "verify_firmware", lambda _: {})
    monkeypatch.setattr(sitl_session.subprocess, "Popen", lambda command, **kw: calls.append(command) or Process())
    monkeypatch.setattr(sitl_session.SITLManager, "_watch", lambda _: None)
    manager = sitl_session.SITLManager(root=tmp_path)
    setup = FlightSetup(initial_point="custom_hover", altitude_m=32, duration_s=30)
    directory = manager.start_configured(setup)
    command = calls[0]
    assert command[2] == "spearhead.sitl.quad" and "--interactive" in command
    assert command[command.index("--duration")+1] == "30"
    path = command[command.index("--config")+1]
    from pathlib import Path
    assert yaml.safe_load(Path(path).read_text())["altitude_m"] == 32
    assert not directory.exists()  # Runner creates its own output directory.
    session = json.loads((tmp_path / (directory.name+".session.json")).read_text())
    assert session["flight_setup"] == setup.as_dict()


def test_manager_rejects_mismatched_setup_before_launch(tmp_path):
    manager = sitl_session.SITLManager(root=tmp_path)
    with pytest.raises(ValueError, match="agree"):
        manager.start("fbwa", setup=FlightSetup())


def test_recorded_recipe_survives_reload_without_accepting_corrupt_labels(tmp_path):
    from spearhead.gui.sitl_session import RunObserver
    directory = tmp_path / "flight"
    directory.mkdir()
    recipe = FlightSetup(initial_point="custom_hover", altitude_m=32, duration_s=30)
    path = tmp_path / "flight.session.json"
    path.write_text(json.dumps({"flight_setup": recipe.as_dict()}))
    observer = RunObserver(directory)
    observer.poll()
    assert observer.snapshot()["flight_setup"] == recipe.as_dict()
    path.write_text(json.dumps({"flight_setup": {"aircraft": "missing"}}))
    observer = RunObserver(directory)
    observer.poll()
    assert observer.snapshot()["flight_setup"] is None
