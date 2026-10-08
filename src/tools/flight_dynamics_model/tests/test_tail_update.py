import csv
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from spearhead import load_aircraft_config, select_aerodynamics
from spearhead.aero import aero_model, _placeholder_control_rate_increment
from spearhead.aerodb import dimensional_body_force_moment, load_aerodynamic_database
from spearhead.aircraft import FDM_ROOT, aircraft_summary
from spearhead.controls import Control
from spearhead.deflections import load_deflection_database
from spearhead.gui.app import scenario_options_for_aircraft
from spearhead.gui.flight_deck import FlightDeckSession
from spearhead.scenarios import load_scenario
from spearhead.simulation import run_simulation
from spearhead.trim import trim_fixed_wing_longitudinal

ROOT = FDM_ROOT / "data/aerodynamics/tail_update_2026_10_07"


def source_rows():
    with (ROOT / "sources/deflection.csv").open() as stream:
        reader = csv.DictReader(stream)
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        return [{key: float(value) for key, value in row.items() if key != "case_id"} for row in reader]


def test_preserves_original_deliveries_and_overrides_only_the_conflicting_header():
    source = json.loads((ROOT / "source.json").read_text())
    for name, info in source["sources"].items():
        assert hashlib.sha256((ROOT / info["path"]).read_bytes()).hexdigest() == info["sha256"]
    raw = (ROOT / "sources/database.csv").read_text()
    assert (ROOT / "database.csv").read_text() == raw.replace(
        "moment_center_body_xyz=(-0.1,0,0.15)", "moment_center_body_xyz=(0,0,0)", 1,
    )
    p = load_aircraft_config("pt1_as_built")
    database = load_aerodynamic_database(p.aero_database_path)
    assert len(database.alpha_grid_deg) * len(database.beta_grid_deg) == 32761
    assert p.rho == database.reference.rho == 1.021
    np.testing.assert_equal(database.reference.moment_center_body_xyz, [0, 0, 0])
    assert p.cref == 0.457  # Wing reference chord, not the updated tail chord.
    assert p.geometry.tail_chord is None


@pytest.mark.parametrize("row", source_rows())
def test_all_delivered_cases_reproduce_measured_load_increments(row):
    p = replace(load_aircraft_config("pt1_as_built"), cg_body_xyz=np.zeros(3))
    alpha = np.deg2rad(row["alpha"])
    velocity = row["velocity"] * np.array([np.cos(alpha), 0, np.sin(alpha)])
    zero = Control(0, 0, 0, 0)
    control = replace(zero, de=np.deg2rad(row["elevator"]),
                      da=np.deg2rad(row["aileron"]), dr=np.deg2rad(row["rudder"]))
    base = next(r for r in source_rows() if r["alpha"] == row["alpha"] and
                all(r[key] == 0 for key in ("aileron", "elevator", "rudder")))
    force0, moment0, _ = aero_model(velocity, np.zeros(3), zero, p)
    force, moment, info = aero_model(velocity, np.zeros(3), control, p)
    np.testing.assert_allclose(force - force0, [row[k] - base[k] for k in ("Fx", "Fy", "Fz")], atol=1e-10)
    np.testing.assert_allclose(moment - moment0, [row[k] - base[k] for k in ("Mx", "My", "Mz")], atol=1e-10)
    assert info["control_source"] == "measured_increment_csv"
    assert info["control_data_coverage"]["within_measured_coverage"]


def test_interpolates_sparse_controls_without_double_counting_and_scales_density_speed_cg():
    p = load_aircraft_config("pt1_as_built")
    zero = Control(0, 0, 0, 0)
    half = replace(zero, da=np.deg2rad(2.5))
    full = replace(zero, da=np.deg2rad(5))
    velocity = 25 * np.array([np.cos(np.deg2rad(-2.5)), 0, np.sin(np.deg2rad(-2.5))])
    f0, m0, _ = aero_model(velocity, np.zeros(3), zero, p)
    f1, m1, _ = aero_model(velocity, np.zeros(3), full, p)
    fh, mh, info = aero_model(velocity, np.zeros(3), half, p)
    np.testing.assert_allclose(fh - f0, (f1 - f0) / 2, atol=1e-12)
    np.testing.assert_allclose(mh - m0, (m1 - m0) / 2, atol=1e-12)
    p2 = replace(p, rho=p.rho * 0.8, cg_body_xyz=np.zeros(3))
    fs, ms, _ = aero_model(velocity * 2, np.zeros(3), full, p2)
    fs0, ms0, _ = aero_model(velocity * 2, np.zeros(3), zero, p2)
    np.testing.assert_allclose(fs - fs0, (f1 - f0) * 3.2, atol=1e-10)
    np.testing.assert_allclose(ms - ms0, ((m1 - m0) - np.cross(-p.cg_body_xyz, f1 - f0)) * 3.2, atol=1e-10)
    assert all(value == 0 for value in info["provisional_increment_coefficients"].values())


def test_outside_coverage_holds_endpoints_and_exposes_assumptions():
    p = load_aircraft_config("pt1_as_built")
    database = load_deflection_database(p.deflection_database_path, p.aerodb_reference)
    query = Control(np.deg2rad(20), np.deg2rad(5), 0, 0)
    values, coverage = database.increments_at(20, 2, query)
    bounded, _ = database.increments_at(5, 0, replace(query, de=np.deg2rad(5), da=0))
    aileron, _ = database.increments_at(0, 0, replace(query, de=0))
    assert values == pytest.approx({key: bounded[key] + aileron[key] for key in values})
    assert not coverage["within_measured_coverage"]
    assert coverage["held_surfaces"] == ["aileron", "elevator"]
    assert coverage["beta_held"] and coverage["combined_surfaces_assumed_additive"]
    zero, coverage = database.increments_at(90, 90, Control(0, 0, 0, 0))
    assert all(value == 0 for value in zero.values())


def test_rate_derivatives_remain_provisional_with_measured_controls():
    p = load_aircraft_config("pt1_as_built")
    control = Control(0, np.deg2rad(3), 0, 0)
    velocity = np.array([25, 0, 0])
    rates = np.array([0.1, -0.2, 0.3])
    f0, m0, _ = aero_model(velocity, np.zeros(3), control, p)
    f, m, info = aero_model(velocity, rates, control, p)
    df, dm, _ = _placeholder_control_rate_increment(0, 25, 0.5*p.rho*25**2, rates, Control(0, 0, 0, 0), p)
    np.testing.assert_allclose(f-f0, df, atol=1e-10)
    np.testing.assert_allclose(m-m0, dm, atol=1e-10)
    assert info["rate_source"] == "provisional_derivatives"


def test_hardware_and_databases_can_be_switched_independently():
    reference = load_aircraft_config("reference_design")
    new = select_aerodynamics(reference, "tail_update_2026_10_07")
    old = select_aerodynamics(new, "adb_v1_1")
    assert new.lift_system == reference.lift_system
    assert new.pusher_installed
    assert new.rho == 1.021 and old.rho == 1.09
    assert old.aero_database_path == reference.aero_database_path
    assert old.geometry == reference.geometry
    assert old.deflection_database_path is None
    assert aircraft_summary(new)["control_data"]["coverage"]["elevator"]["alpha_deg"] == [-5, 5]
    assert "pitch_doublet.yaml" not in scenario_options_for_aircraft("reference_design", "tail_update_2026_10_07")


def test_new_trim_does_not_accept_solver_termination_without_balance():
    p = load_aircraft_config("reference_design", aerodynamics="tail_update_2026_10_07")
    trim = trim_fixed_wing_longitudinal(p, target_speed=25)
    assert not trim["success"]
    assert -5 <= trim["de_deg"] <= 5
    assert "balance" in trim["message"]
    config = load_scenario(FDM_ROOT / "scenarios/pitch_doublet.yaml", aerodynamics="tail_update_2026_10_07")
    with pytest.raises(ValueError, match="Trim failed"):
        run_simulation(config)


def test_open_loop_new_database_runs_and_gui_switch_resets_history():
    config = load_scenario(FDM_ROOT / "scenarios/aero_database_check.yaml")
    result = run_simulation(config)
    assert result.success and np.all(np.isfinite(result.states))
    assert config.params.aero_dataset_id == "tail_update_2026_10_07"
    session = FlightDeckSession("new", config)
    session.step_once()
    assert session.aerodynamic_coverage() is not None
    old = load_scenario(FDM_ROOT / "scenarios/pitch_doublet.yaml")
    session.set_scenario("old", old)
    assert session.current_state.time == 0
    assert session.aerodynamic_coverage() is None


def test_missing_zero_baseline_is_rejected(tmp_path):
    text = (ROOT / "sources/deflection.csv").read_text()
    path = tmp_path / "bad.csv"
    path.write_text('\n'.join(line for line in text.splitlines() if not line.startswith('Case008,')))
    with pytest.raises(ValueError, match="zero-deflection baseline"):
        load_deflection_database(path, load_aircraft_config("pt1_as_built").aerodb_reference)
