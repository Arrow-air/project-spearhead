import json

import pytest
import yaml

from spearhead.gui.flight_setup import FlightSetup
from spearhead.gui import sitl_session
from spearhead.sitl.backend import load_setup


@pytest.mark.parametrize("speed", [22, 24, 26])
@pytest.mark.parametrize("cg_case", ["target40", "pt1_log", "aft60", "forward25", "forward0"])
def test_prototype_forward_start_solves_actual_trim_with_requested_actuators(tmp_path, speed, cg_case):
    import numpy as np
    from spearhead.controls import Control
    from spearhead.dynamics import dynamics
    from spearhead.sitl.flight import prototype_forward_parameters
    from spearhead.sitl.transition import transition_parameters
    setup = FlightSetup(initial_point=f"trim{speed}", speed_m_s=speed, altitude_m=100,
                        maneuver="hold", cg_case=cg_case)
    assert setup.scenario == "fbwa"
    path = tmp_path / "forward.yaml"
    path.write_text(yaml.safe_dump(setup.config()))
    physics, mapping, manifest = load_setup(path)
    trim = Control(**manifest["trim_control"])
    residual = dynamics(0, physics.state, physics.params, trim)
    np.testing.assert_allclose(residual[[2, 3, 4, 5, 9, 10, 11]], 0, atol=1e-6)
    assert np.linalg.norm(physics.state[3:6]) == pytest.approx(speed)
    assert physics.params.mass == 21 and mapping.lift == ()
    assert physics.params.control_surface_limits_deg == (30, 30, 30)
    assert (mapping.pusher.minimum, mapping.pusher.maximum) == (1050, 1940)
    params = prototype_forward_parameters(mapping, physics.state, trim)
    baseline = transition_parameters("first_transition_2026_10_09")
    for key in baseline:
        if key.startswith(("RLL_RATE_", "PTCH_RATE_")):
            assert params[key] == baseline[key]
    pwm = [params["SERVO1_TRIM"], params["SERVO2_TRIM"], 1050, params["SERVO4_TRIM"]]
    actual = mapping.control(pwm)
    assert actual.throttle == 0 and actual.lift_throttle == 0
    assert abs(actual.de-trim.de) <= np.deg2rad(.031)


def test_forward_and_hover_maneuvers_are_separate_and_old_hover_recording_survives(tmp_path):
    from spearhead.gui.flight_setup import POINTS, maneuver_options
    from spearhead.gui.sitl_session import RunObserver
    assert "hover30" not in POINTS["pt15"]
    assert maneuver_options("pt15", "trim24")["hold"] == "Straight & level hold"
    with pytest.raises(ValueError, match="maneuver"):
        FlightSetup(initial_point="trim24", speed_m_s=24, altitude_m=100, maneuver="transition").validate()
    with pytest.raises(ValueError, match="26"):
        FlightSetup(initial_point="custom_trim", speed_m_s=28, altitude_m=100).validate()
    directory = tmp_path / "recording"
    directory.mkdir()
    saved = FlightSetup(initial_point="custom_hover", altitude_m=30).as_dict()
    saved["initial_point"] = "hover30"
    (tmp_path / "recording.session.json").write_text(json.dumps({"flight_setup": saved}))
    observer = RunObserver(directory)
    observer.poll()
    assert observer.flight_setup["initial_point"] == "custom_hover"
    assert observer.flight_setup["altitude_m"] == 30


@pytest.mark.parametrize("setup,aircraft,speed,altitude,scenario", [
    (FlightSetup(), "pt15", 0, 20, "quad_interactive"),
    (FlightSetup(initial_point="custom_hover", altitude_m=32, duration_s=30, maneuver="hold"), "pt15", 0, 32, "quad_hover"),
    (FlightSetup("pt15", "trim26", "hold", 60, 120, 26), "pt15", 26, 120, "fbwa"),
    (FlightSetup("pt15", "custom_trim", "roll_doublet", 45, 100, 25), "pt15", 25, 100, "roll_doublet"),
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
    with pytest.raises(ValueError, match="45"):
        FlightSetup("pt15", "trim24", "pitch_doublet", 30, 100, 24).validate()
    with pytest.raises(ValueError, match="Custom forward"):
        FlightSetup("pt15", "trim24", "hold", 60, 100, 22).validate()
    with pytest.raises(ValueError, match="aircraft"):
        FlightSetup("reference_design", "trim28", "hold", 60, 100, 28).validate()


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


def test_manager_routes_prototype_forward_start_to_fixed_wing_runner(tmp_path, monkeypatch):
    from types import SimpleNamespace
    calls = []
    monkeypatch.setattr(sitl_session, "verify_firmware", lambda _: {})
    monkeypatch.setattr(sitl_session.subprocess, "Popen", lambda command, **kw: calls.append(command) or SimpleNamespace(poll=lambda: None))
    monkeypatch.setattr(sitl_session.threading.Thread, "start", lambda _: None)
    setup = FlightSetup(initial_point="trim24", speed_m_s=24, altitude_m=100, cg_case="forward25")
    sitl_session.SITLManager(root=tmp_path).start_configured(setup)
    command = calls[0]
    assert command[2] == "spearhead.sitl.flight" and "--interactive" in command
    assert command[command.index("--forward-parameters")+1] == "first_transition_2026_10_09"


def test_manager_rejects_mismatched_setup_before_launch(tmp_path):
    manager = sitl_session.SITLManager(root=tmp_path)
    with pytest.raises(ValueError, match="agree"):
        manager.start("fbwa", setup=FlightSetup())


@pytest.mark.parametrize("cg_case,x", [("pt1_log", -.0475), ("aft60", -.06), ("forward25", -.025), ("forward0", 0.)])
def test_transition_cg_reaches_physics_without_changing_source_origin(tmp_path, cg_case, x):
    setup = FlightSetup(initial_point="hover80", altitude_m=80, maneuver="transition", cg_case=cg_case)
    config = setup.config()
    assert config["cg_body_xyz_m"] == [x, 0, .065]
    assert config["hover_trim"] == "independent_rotors"
    path = tmp_path / "transition.yaml"
    path.write_text(yaml.safe_dump(config))
    physics, _, manifest = load_setup(path)
    assert physics.params.cg_body_xyz.tolist() == [x, 0, .065]
    assert manifest["aircraft"]["moment_center_body_xyz_m"] == [0, 0, 0]


@pytest.mark.parametrize("parameters", ["first_transition_2026_10_09"])
def test_manager_launches_exact_transition_parameter_choice(tmp_path, monkeypatch, parameters):
    from types import SimpleNamespace
    calls = []
    monkeypatch.setattr(sitl_session, "verify_firmware", lambda _: {})
    monkeypatch.setattr(sitl_session.subprocess, "Popen", lambda command, **kw: calls.append(command) or SimpleNamespace(poll=lambda: None))
    monkeypatch.setattr(sitl_session.threading.Thread, "start", lambda _: None)
    manager = sitl_session.SITLManager(root=tmp_path)
    setup = FlightSetup(initial_point="hover80", altitude_m=80, maneuver="transition",
                        cg_case="forward25", transition_parameters=parameters)
    directory = manager.start_configured(setup)
    command = calls[0]
    assert "--transition" in command and "--interactive" not in command
    assert command[command.index("--transition-parameters")+1] == parameters
    session = json.loads((tmp_path / (directory.name+".session.json")).read_text())
    assert session["flight_setup"] == setup.as_dict()


@pytest.mark.parametrize("changes", [{"cg_case": "unknown"}, {"transition_parameters": "custom_tune"}])
def test_invalid_transition_choices_rejected(changes):
    with pytest.raises(ValueError, match="supported transition"):
        FlightSetup(**changes).validate()


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


def test_older_transition_recording_retains_its_parameter_identity(tmp_path):
    from spearhead.gui.sitl_session import RunObserver
    directory = tmp_path / "flight"
    directory.mkdir()
    saved = FlightSetup(initial_point="hover80", altitude_m=80, maneuver="transition").as_dict()
    saved.pop("transition_parameters")
    saved.pop("cg_case")
    (tmp_path / "flight.session.json").write_text(json.dumps({"flight_setup": saved}))
    (directory / "manifest.json").write_text(json.dumps({"transition_plan": {}}))
    observer = RunObserver(directory)
    observer.poll()
    assert observer.snapshot()["flight_setup"]["transition_parameters"] == "first_transition_2026_10_09"


@pytest.mark.parametrize("point,maneuver,altitude,speed", [
    ("hover20", "hold", 20, 0),
    ("hover80", "transition", 80, 0),
    ("trim24", "hold", 100, 24),
    ("ground", "takeoff_hover_land", 5, 0),
])
def test_forward_default_cg_is_used_for_every_prototype_start(tmp_path, point, maneuver, altitude, speed):
    import numpy as np
    from spearhead.controls import Control
    from spearhead.dynamics import dynamics
    setup = FlightSetup(initial_point=point, maneuver=maneuver, altitude_m=altitude, speed_m_s=speed)
    assert setup.cg_case == "target40"
    path = tmp_path / "setup.yaml"
    path.write_text(yaml.safe_dump(setup.config()))
    physics, _, manifest = load_setup(path)
    assert physics.params.cg_body_xyz.tolist() == [-.04, 0, .065]
    assert manifest["aircraft"]["moment_center_body_xyz_m"] == [0, 0, 0]
    control = Control(**manifest["trim_control"])
    if point != "ground":
        np.testing.assert_allclose(dynamics(0, physics.state, physics.params, control)[[2, 3, 4, 5, 9, 10, 11]], 0, atol=1e-6)
    if not setup.is_forward:
        assert setup.config()["hover_trim"] == "independent_rotors"
        assert control.lift_motor_throttles is not None


def test_wind_gust_and_failure_reach_the_runner_configuration(tmp_path):
    setup = FlightSetup(initial_point="hover20", maneuver="hold", wind_speed_m_s=6, wind_from_deg=270,
                        gust_m_s=3, failure="lift_fr_out", failure_time_s=30)
    config = setup.config()
    assert config["wind"] == {"speed_m_s": 6., "from_deg": 270., "gusts": [
        {"start_s": 20., "duration_s": 3., "speed_m_s": 3., "from_deg": 270.}]}
    assert config["failure"] == {"kind": "lift_fr_out", "start_s": 30.}
    path = tmp_path / "conditions.yaml"; path.write_text(yaml.safe_dump(config))
    physics, _, manifest = load_setup(path)
    assert physics.wind.speed_m_s == 6 and physics.failure.kind == "lift_fr_out"
    calm = FlightSetup().config()
    assert "wind" not in calm and "failure" not in calm


@pytest.mark.parametrize("changes,message", [
    ({"wind_speed_m_s": 20}, "Wind"),
    ({"failure": "battery_low"}, "approval"),
    ({"failure": "pusher_out", "maneuver": "hold"}, "Pusher"),
    ({"failure": "lift_fr_out", "initial_point": "trim24", "speed_m_s": 24, "altitude_m": 100, "maneuver": "hold"}, "Lift-motor"),
    ({"failure": "lift_fr_out", "failure_time_s": 70}, "inside the flight"),
])
def test_invalid_conditions_rejected(changes, message):
    with pytest.raises(ValueError, match=message):
        FlightSetup(**changes).validate()


def test_recordings_with_old_geometry_are_labelled(tmp_path):
    from spearhead.gui.sitl_session import recording_label
    new = tmp_path / "new"; new.mkdir()
    (new / "manifest.json").write_text(json.dumps({"aircraft": {"id": "pt15", "lift_positions_body_m": [[.675, .55, -.043]]}}))
    old = tmp_path / "old"; old.mkdir()
    (old / "manifest.json").write_text(json.dumps({"aircraft": {"id": "pt1_as_built"}}))
    assert recording_label("new", new) == "new"
    assert recording_label("old", old) == "old · legacy geometry (before 14 Oct)"
    assert recording_label("missing", tmp_path / "none") == "missing"
