import csv
import hashlib
import json
from dataclasses import replace
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import numpy as np
import pytest
import yaml

from spearhead import AircraftParams, load_aircraft_config, run_simulation
from spearhead.aircraft import AIRCRAFT_CONFIGS, aircraft_summary
from spearhead.controls import Control
from spearhead.force_moment import compute_force_moment_breakdown
from spearhead.gui.app import scenario_options_for_aircraft
from spearhead.gui.flight_deck import FlightDeckSession, maneuver_options
from spearhead.propulsion import hover_collective, lift_propulsion, load_static_performance, propulsion_model
from spearhead.scenarios import load_scenario
from spearhead.trim import trim_fixed_wing_longitudinal
from spearhead.aero import aero_model

FDM_ROOT = AIRCRAFT_CONFIGS.parent


def test_manufacturer_csv_matches_all_28_workbook_samples():
    """Check units/column selection against the preserved manufacturer workbook."""
    root = FDM_ROOT / "data/propulsion"
    workbook = root / "sources/KDE5215-220kV_rev2.xlsx"
    source = json.loads((root / "kde5215xf_220_source.json").read_text())
    assert hashlib.sha256(workbook.read_bytes()).hexdigest() == source["workbook_sha256"]
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(workbook) as archive:
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    values = {cell.attrib["r"]: cell.find("x:v", ns).text for cell in sheet.findall(".//x:c", ns) if cell.find("x:v", ns) is not None}
    with (root / "kde5215xf_220_cf215_dp.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 28
    for row in rows:
        number = row["source_row"]
        for key, column in {"throttle": "D", "thrust_n": "I", "current_a": "E", "power_w": "F", "rpm": "N", "torque_nm": "O", "estimated_sea_level_thrust_n": "L"}.items():
            assert float(row[key]) == pytest.approx(float(values[f"{column}{number}"]))


def test_profiles_are_independent_and_preserve_original_parameters():
    original = AircraftParams()
    reference = load_aircraft_config("reference_design")
    pt1 = load_aircraft_config("pt15")
    assert pt1.mass == 21.0 and reference.mass == 25.0
    assert reference.lift_system.motor["model"] == "MAD V8013 PRO IPE"
    assert pt1.lift_system.motor["model"] == "KDE5215XF-220"
    assert pt1.lift_system.esc["model"] == "KDE-UAS75HVC+"
    assert pt1.lift_system.propeller["blades"] == 2
    assert reference.pusher_installed and pt1.pusher_installed
    assert reference.pusher_system is None and pt1.pusher_system.motor["kv_rpm_per_v"] == 150
    assert reference.Tmax == original.Tmax
    assert reference.aero_database_path == original.aero_database_path
    assert reference.geometry == original.geometry
    np.testing.assert_equal(reference.cg_body_xyz, original.cg_body_xyz)
    reference.cg_body_xyz[0] = -0.2
    assert pt1.cg_body_xyz.tolist() == [-.04, 0, .065]


def test_pt1_as_flown_and_pt15_share_drawing_rotors_but_not_mass_cg_or_pusher():
    pt1 = load_aircraft_config("pt1_as_built")
    pt15 = load_aircraft_config("pt15")
    assert pt1.mass == 17.9 and pt15.mass == 21.0
    assert pt1.cg_body_xyz.tolist() == [-.0475, 0, .065]
    assert not pt1.pusher_installed and pt15.pusher_installed
    assert pt15.pusher_system.position_body_m == (-.374, 0., .05)
    # 14 October drawing, X aft / Z up converted to body x forward / z down.
    expected = ((.675, .55, -.043), (-.75, .55, -.043), (-.75, -.55, -.043), (.675, -.55, -.043))
    assert pt1.lift_system.positions_body_m == pt15.lift_system.positions_body_m == expected
    centroid = np.mean(np.array(expected), axis=0)
    assert centroid[0] == pytest.approx(-.0375)
    # The logged 49.3 % front hover share maps to the PT1 CG with these stations.
    assert .493*.675 + .507*-.75 == pytest.approx(pt1.cg_body_xyz[0], abs=5e-4)
    assert pt1.aero.Cmq == -9.6 and pt15.aero.Cnr == -.04
    assert aircraft_summary(pt15)["rate_derivatives"]["Clp"] == -.58


def test_profile_rate_derivatives_are_validated(tmp_path):
    data = yaml.safe_load((AIRCRAFT_CONFIGS / "pt1_as_built.yaml").read_text())
    data["propulsion"]["lift"]["performance_path"] = str(load_aircraft_config("pt1_as_built").lift_system.performance_path)
    data["aerodynamics"]["rate_derivatives"] = {"Cmalpha": -1.}
    profile = tmp_path / "bad.yaml"
    profile.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError, match="rate derivatives"):
        load_aircraft_config(profile)


def test_zero_lift_reference_loads_equal_legacy_forward_propulsion():
    reference = load_aircraft_config("reference_design")
    control = Control(0.0, 0.0, 0.0, 0.45)
    force, moment = propulsion_model(control, reference)
    np.testing.assert_allclose(force, [54.0, 0.0, 0.0])
    np.testing.assert_allclose(moment, [0.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="No lift performance curve"):
        propulsion_model(replace(control, lift_throttle=0.5), reference)


def test_static_curve_exact_points_interpolation_and_voltage_bounds():
    pt1 = load_aircraft_config("pt15")
    curve = load_static_performance(pt1.lift_system.performance_path)
    point = curve.lookup(50.4, 0.5)
    assert point == {"thrust_n": 42.39, "current_a": 13.4, "power_w": 676.0, "rpm": 4460.0, "torque_nm": 1.11}
    assert curve.lookup(46.2, 0.5)["thrust_n"] == pytest.approx((31.9 + 42.39) / 2)
    assert curve.lookup(50.4, 0.5625)["thrust_n"] == pytest.approx((42.39 + 56.36) / 2)
    assert curve.lookup(50.4, -1)["thrust_n"] == 0
    assert curve.lookup(50.4, 2)["thrust_n"] == 95.88
    with pytest.raises(ValueError, match="outside the curve range"):
        curve.lookup(60.0, 0.5)


def test_lift_force_body_sign_moment_balance_and_cg_transfer():
    pt1 = load_aircraft_config("pt15")
    pt1 = replace(pt1, battery=replace(pt1.battery, voltage_v=50.4))
    control = Control(0, 0, 0, 0, lift_throttle=0.5)
    force, moment = propulsion_model(control, pt1)
    np.testing.assert_allclose(force, [0, 0, -4 * 42.39])
    # Rotor centroid x = -37.5 mm is 2.5 mm ahead of the -40 mm PT1.5 target CG.
    np.testing.assert_allclose(moment, [0, (-.0375 - -.04) * 4 * 42.39, 0], atol=1e-12)
    shifted = replace(pt1, cg_body_xyz=pt1.cg_body_xyz + [0.1, 0.0, 0.0])
    shifted_force, shifted_moment = propulsion_model(control, shifted)
    np.testing.assert_allclose(shifted_force, force)
    np.testing.assert_allclose(shifted_moment - moment, np.cross([-0.1, 0.0, 0.0], force), atol=1e-12)
    breakdown = compute_force_moment_breakdown(0, np.zeros(12), pt1, control_override=control)
    np.testing.assert_allclose(breakdown["propulsion"]["force_b"], force)


def test_single_rotor_has_expected_roll_pitch_and_reaction_yaw_moment():
    pt1 = load_aircraft_config("pt15")
    pt1 = replace(pt1, battery=replace(pt1.battery, voltage_v=50.4))
    control = Control(0, 0, 0, 0, lift_motor_throttles=(0.5, 0, 0, 0))
    force, moment, samples = lift_propulsion(control, pt1)
    np.testing.assert_allclose(force, [0, 0, -42.39])
    np.testing.assert_allclose(moment, [-0.55 * 42.39, (.675 + .04) * 42.39, 1.11])
    assert samples[0]["rpm"] == 4460
    with pytest.raises(ValueError, match="one finite lift throttle"):
        lift_propulsion(replace(control, lift_motor_throttles=(0.5,)), pt1)


def test_pt1_scenario_recomputes_collective_for_mass_and_preserves_control_commands():
    path = FDM_ROOT / "scenarios/pt1_hover.yaml"
    config = load_scenario(path)
    assert config.params.aircraft_id == "pt1_as_built"
    assert not config.start_from_trim
    force, _ = propulsion_model(config.control_override, config.params)
    assert -force[2] == pytest.approx(config.params.mass * config.params.g)
    lighter = replace(config.params, mass=config.params.mass - 1)
    assert hover_collective(lighter) < config.control_override.lift_throttle
    config = replace(config, duration=0.1, dt=0.02)
    result = run_simulation(config)
    assert result.success and np.all(np.isfinite(result.states))
    assert result.control_override(0).lift_throttle == config.control_override.lift_throttle
    assert result.control_override(2.1).lift_throttle == pytest.approx(config.control_override.lift_throttle + 0.03)
    initial_rotors = np.array(config.control_override.lift_motor_throttles)
    for time, delta in [(0, 0), (2.1, .03), (3.1, -.03), (4.1, 0)]:
        command = result.control_override(time)
        np.testing.assert_allclose(command.lift_motor_throttles, initial_rotors + delta)
        actual_force, _, _ = lift_propulsion(command, config.params)
        if delta:
            assert (-actual_force[2] - config.params.mass * config.params.g) * delta > 0


@pytest.mark.parametrize("use_database", [False, True])
def test_aerodynamic_loads_vanish_at_rest_even_with_rate_regularization(use_database):
    params = replace(load_aircraft_config("pt15"), use_aero_database=use_database)
    force, moment, info = aero_model(np.zeros(3), np.array([0.1, 0.2, 0.3]), Control(0.1, 0.1, 0.1, 0), params)
    np.testing.assert_allclose(force, 0)
    np.testing.assert_allclose(moment, 0)
    assert info["qbar"] == 0


def test_pt1_realtime_holds_stationary_static_balance_before_the_input():
    from spearhead.realtime import RealtimeSimulation

    config = load_scenario(FDM_ROOT / "scenarios/pt1_hover.yaml")
    sim = RealtimeSimulation(config)
    for _ in range(10):
        sim.step()
    np.testing.assert_allclose(sim.state.vector, config.initial_state, atol=1e-10)
    assert sim.time == pytest.approx(0.2)


def test_gui_profile_switch_resets_simulation_and_filters_incompatible_scenarios():
    reference = load_aircraft_config("reference_design")
    reference_config = load_scenario(FDM_ROOT / "scenarios/pitch_doublet.yaml", params=reference)
    pt1_config = load_scenario(FDM_ROOT / "scenarios/pt1_hover.yaml")
    session = FlightDeckSession("reference", replace(reference_config, duration=0.1, dt=0.02))
    session.play()
    session.advance(0.04)
    session.set_scenario("pt1", pt1_config)
    assert session.current_state.time == 0
    assert session.base_config.params.aircraft_id == "pt1_as_built"
    assert session.controls()["pusher_thrust_n"] == 0
    assert session.controls()["lift_thrust_n"] == pytest.approx(17.9 * 9.80665)
    assert "lift_step" in maneuver_options(pt1_config)
    assert "pt1_hover.yaml" in scenario_options_for_aircraft("pt1_as_built")
    assert "pt1_hover.yaml" in scenario_options_for_aircraft("pt15")
    assert "pitch_doublet.yaml" in scenario_options_for_aircraft("reference_design")
    session.set_scenario("reference", reference_config)
    assert session.base_config.params.lift_system.motor["model"] == "MAD V8013 PRO IPE"
    assert session.controls()["lift_thrust_n"] == 0
    with pytest.raises(ValueError, match="requires an installed pusher"):
        trim_fixed_wing_longitudinal(replace(pt1_config.params, pusher_installed=False))


def test_new_geometry_and_aerodb_can_be_selected_without_overwriting_old_profile(tmp_path):
    data = yaml.safe_load((AIRCRAFT_CONFIGS / "pt15.yaml").read_text())
    data["geometry"]["tail_chord"] = 0.30
    data["aerodynamics"] = {"version": "future-test-db", "geometry_version": "test_tail_0.30m", "path": "new.csv"}
    data["propulsion"]["lift"]["performance_path"] = str(load_aircraft_config("pt15").lift_system.performance_path)
    data["propulsion"]["pusher"]["estimate"]["performance_path"] = str(load_aircraft_config("pt15").pusher_system.performance_path)
    (tmp_path / "new.csv").write_text("source supplied later\n")
    profile = tmp_path / "new.yaml"
    profile.write_text(yaml.safe_dump(data))
    updated = load_aircraft_config(profile)
    assert updated.geometry.tail_chord == 0.30
    assert updated.aero_database_path == tmp_path / "new.csv"
    assert aircraft_summary(updated)["geometry_version"] == "test_tail_0.30m"
    old = load_aircraft_config("pt15", aerodynamics="adb_v1_1")
    assert old.geometry.tail_chord == 0.28
    assert old.aero_database_version == "adb_v1.1"
