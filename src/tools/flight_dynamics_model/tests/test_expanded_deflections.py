"""Check delivered dimensional increments, sparse coverage and frozen gains."""
import csv
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

from spearhead.aircraft import FDM_ROOT, load_aircraft_config
from spearhead.aero import aero_model
from spearhead.controls import Control
from spearhead.deflections import load_deflection_database
from spearhead.sitl.backend import load_setup, CONFIG_ROOT
from spearhead.sitl.transition import transition_parameters, forward_trim_parameters

SOURCE = FDM_ROOT / "data/aerodynamics/tail_update_2026_10_07/sources/deflection_2026_10_11.csv"
ROWS = list(csv.DictReader(SOURCE.open()))
CASES = [r for r in ROWS if float(r["beta"]) == 0]


@pytest.mark.parametrize("row", CASES, ids=[r["case_id"] for r in CASES])
def test_every_new_measured_case_reproduces_all_six_load_increments(row):
    params = replace(load_aircraft_config("pt1_as_built"), cg_body_xyz=np.zeros(3),
                     control_deflection_policy="linear", control_surface_limits_deg=(30,30,30))
    r = {k: float(v) for k,v in row.items() if k != "case_id"}
    base = next(b for b in ROWS if all(float(b[k]) == r[k] for k in ("alpha", "beta", "velocity"))
                and all(float(b[k]) == 0 for k in ("elevator", "aileron", "rudder")))
    a = np.deg2rad(r["alpha"])
    velocity = r["velocity"] * np.array([np.cos(a),0,np.sin(a)])
    control = Control(*np.deg2rad([r["elevator"],r["aileron"],r["rudder"]]),0)
    f0,m0,_ = aero_model(velocity,np.zeros(3),Control(0,0,0,0),params)
    f,m,info = aero_model(velocity,np.zeros(3),control,params)
    expected = [r[k]-float(base[k]) for k in ("Fx","Fy","Fz","Mx","My","Mz")]
    np.testing.assert_allclose(np.r_[f-f0,m-m0],expected,atol=1e-10)
    assert info["control_data_coverage"]["within_measured_coverage"]


def test_speed_planes_interpolate_coefficients_only_inside_shared_coverage():
    p = load_aircraft_config("pt1_as_built")
    db = load_deflection_database(p.deflection_database_path,p.aerodb_reference)
    command = Control(np.deg2rad(5),0,0,0)
    low,_ = db.increments_at(10,0,command,velocity_m_s=10)
    high,_ = db.increments_at(10,0,command,velocity_m_s=15)
    middle,coverage = db.increments_at(10,0,command,velocity_m_s=12.5)
    assert middle == pytest.approx({k:(low[k]+high[k])/2 for k in low})
    assert coverage["within_measured_coverage"]
    # Large deflection has no low-speed measurement; extrapolate only when the
    # configured linear policy requests it, flagging the participating planes.
    large = replace(command,de=np.deg2rad(12))
    extrapolated,coverage = db.increments_at(10,0,large,velocity_m_s=12.5,
                                           deflection_policy="linear",limits_deg=(20,20,20))
    low,_ = db.increments_at(10,0,large,velocity_m_s=10,deflection_policy="linear",limits_deg=(20,20,20))
    high,_ = db.increments_at(10,0,large,velocity_m_s=15,deflection_policy="linear",limits_deg=(20,20,20))
    assert extrapolated == pytest.approx({k:(low[k]+high[k])/2 for k in low})
    assert coverage["extrapolated_surfaces"] == ["elevator"]
    assert not coverage["within_measured_coverage"]


@pytest.mark.parametrize("alpha,angle", [(0.,5.),(0.,12.),(10.,5.)])
def test_low_speed_sparse_lookup_has_no_force_jump_at_coverage_boundary(alpha,angle):
    p = load_aircraft_config("pt1_as_built")
    db = load_deflection_database(p.deflection_database_path,p.aerodb_reference)
    def query(a,d):
        values,_ = db.increments_at(a,0,Control(np.deg2rad(d),0,0,0),velocity_m_s=20,
                                    deflection_policy="linear",limits_deg=(20,20,20))
        return np.array(list(values.values()))
    np.testing.assert_allclose(query(alpha-1e-7,angle),query(alpha+1e-7,angle),atol=1e-7)
    np.testing.assert_allclose(query(alpha,angle-1e-7),query(alpha,angle+1e-7),atol=1e-7)


def test_sparse_hull_extrapolation_and_assumed_travel_are_visible():
    p = load_aircraft_config("pt1_as_built")
    db = load_deflection_database(p.deflection_database_path,p.aerodb_reference)
    outside,coverage = db.increments_at(-5,0,Control(np.deg2rad(15),0,0,0),
                                       deflection_policy="linear",limits_deg=(20,20,20),velocity_m_s=25)
    assert coverage["extrapolated_surfaces"] == ["elevator"]
    assert not coverage["within_measured_coverage"]
    capped,cov = db.increments_at(0,0,Control(np.deg2rad(40),0,0,0),
                                 deflection_policy="linear",limits_deg=(20,20,20))
    limit,_ = db.increments_at(0,0,Control(np.deg2rad(20),0,0,0),
                              deflection_policy="linear",limits_deg=(20,20,20))
    assert capped == pytest.approx(limit)
    assert cov["travel_limited_surfaces"] == ["elevator"]


def test_sideslip_baselines_are_preserved_without_replacing_base_database():
    p = load_aircraft_config("pt1_as_built")
    db = load_deflection_database(p.deflection_database_path,p.aerodb_reference)
    assert db.case_count == 103 and db.audit_cases == [f"Case{n:03d}" for n in range(50,56)]
    value,cov = db.increments_at(20,5,Control(0,0,0,0),velocity_m_s=25)
    assert all(v == 0 for v in value.values())
    controlled,cov = db.increments_at(20,5,Control(np.deg2rad(5),0,0,0),velocity_m_s=25)
    assert cov["beta_held"] and not cov["within_measured_coverage"]


def test_interaction_residual_vanishes_when_either_control_is_neutral():
    p = load_aircraft_config("pt1_as_built")
    db = load_deflection_database(p.deflection_database_path,p.aerodb_reference)
    for control in (Control(0,0,np.deg2rad(15),0),Control(np.deg2rad(-10),0,0,0)):
        _,coverage = db.increments_at(0,0,control,velocity_m_s=25)
        assert not coverage["interaction_residual_applied"]
    _,coverage = db.increments_at(9,0,Control(np.deg2rad(-10),0,np.deg2rad(7),0),velocity_m_s=25)
    assert coverage["interaction_residual_applied"] and not coverage["interaction_at_measured_case"]
    assert not coverage["within_measured_coverage"]


def test_earlier_proposal_matches_all_24_values_and_does_not_copy_simulation_gains():
    earlier = transition_parameters("first_transition_2026_10_09")
    candidate = json.loads((Path(__file__).parent/"fixtures/first_transition_2026_10_09.json").read_text())
    for name,entry in candidate["proposed"].items():
        assert earlier[name] == entry["value"]
    physics,mapping,_ = load_setup(CONFIG_ROOT/"quad_transition.yaml")
    trims,_,_ = forward_trim_parameters(physics.params,mapping,80,include_simulation_gains=False)
    assert not any(k.startswith(("PTCH_RATE_","RLL_RATE_")) or k == "PTCH2SRV_TCONST" for k in trims)
    assert earlier["PTCH_RATE_FF"] == .345 and earlier["PTCH_RATE_P"] == .04
    assert earlier["RLL_RATE_I"] == .15 and earlier["PTCH2SRV_TCONST"] == .75
