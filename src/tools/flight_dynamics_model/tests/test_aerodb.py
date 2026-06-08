from pathlib import Path

import numpy as np
import pytest

from spearhead import aero as aero_module
from spearhead.controls import Control
from spearhead.aerodb import (
    dimensional_body_force_moment,
    load_aerodynamic_database,
)
from spearhead.params import AeroDBReference, AircraftParams


PROJECT_ROOT = Path(__file__).resolve().parents[4]
PR7_ADB_PATH = (
    PROJECT_ROOT
    / "docs/information-note/phase-1/Aerodynamics/"
    "0002-Preliminary-Design-Aerodynamic-Database/assets/"
    "adb_v1_1_drag_correction.csv"
)


def _write_sample_adb(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# Nondimit export,rho=1.09,S=1.805,V=25,q=340.625,span=3.95,chord=0.457,output_accuracy=0.0001",
                '# Moment center,"moment_center_body_xyz=(-0.15,0,0.15)",body_axes=x_forward_y_right_z_down',
                "# Source and options,alpha_column=Alpha(deg),beta_column=Beta(deg)",
                "alpha_deg,beta_deg,cd,cs,cl,cfx,cfy,cfz,cmx,cmy,cmz,cr,cm,cn",
                "0,0,0.01,0.00,0.10,1.0,2.0,3.0,0.10,0.20,0.30,0.0,0.0,0.0",
                "0,1,0.02,0.01,0.20,2.0,3.0,4.0,0.20,0.30,0.40,0.0,0.0,0.0",
                "1,0,0.03,0.02,0.30,3.0,4.0,5.0,0.30,0.40,0.50,0.0,0.0,0.0",
                "1,1,0.04,0.03,0.40,4.0,5.0,6.0,0.40,0.50,0.60,0.0,0.0,0.0",
            ]
        )
    )


def test_loads_adb_with_metadata_comment_rows(tmp_path):
    adb_path = tmp_path / "adb.csv"
    _write_sample_adb(adb_path)

    database = load_aerodynamic_database(adb_path)

    assert database.metadata["Nondimit export"].startswith("rho=1.09")
    assert database.alpha_grid_deg.tolist() == [0.0, 1.0]
    assert database.beta_grid_deg.tolist() == [0.0, 1.0]
    assert database.reference.rho == 1.09
    np.testing.assert_allclose(database.reference.moment_center_body_xyz, [-0.15, 0.0, 0.15])


def test_interpolates_at_exact_grid_point(tmp_path):
    adb_path = tmp_path / "adb.csv"
    _write_sample_adb(adb_path)
    database = load_aerodynamic_database(adb_path)

    coefficients = database.coefficients_at(1.0, 1.0)

    assert coefficients["cfx"] == 4.0
    assert coefficients["cfy"] == 5.0
    assert coefficients["cfz"] == 6.0
    assert coefficients["cmy"] == 0.50


def test_body_axis_force_dimensionalization_signs():
    reference = AeroDBReference()
    coefficients = {"cfx": 1.0, "cfy": 2.0, "cfz": 3.0, "cmx": 0.0, "cmy": 0.0, "cmz": 0.0}

    force_b, _ = dimensional_body_force_moment(
        coefficients,
        reference.q,
        reference,
        shift_to_cg=False,
    )

    np.testing.assert_allclose(force_b, reference.qS * np.array([1.0, 2.0, 3.0]), rtol=1e-6)


def test_default_aircraft_params_match_current_adb_reference():
    params = AircraftParams()

    assert params.rho == 1.09
    assert params.aerodb_reference.rho == 1.09
    assert params.aerodb_reference.q == 340.625
    np.testing.assert_allclose(params.cg_body_xyz, [-0.15, 0.0, 0.15])
    np.testing.assert_allclose(params.aerodb_reference.moment_center_body_xyz, [-0.15, 0.0, 0.15])


def test_moment_center_shift_from_adb_center_to_cg():
    reference = AeroDBReference(moment_center_body_xyz=np.array([-0.15, 0.0, 0.15]))
    coefficients = {"cfx": 1.0, "cfy": 2.0, "cfz": 3.0, "cmx": 0.1, "cmy": 0.2, "cmz": 0.3}
    cg_body_xyz = np.array([-0.05, 0.0, 0.10])

    force_b, moment_b = dimensional_body_force_moment(
        coefficients,
        reference.q,
        reference,
        shift_to_cg=True,
        cg_body_xyz=cg_body_xyz,
    )

    moment_before_shift = np.array([0.1 * reference.qSb, 0.2 * reference.qSc, 0.3 * reference.qSb])
    expected = moment_before_shift + np.cross(reference.moment_center_body_xyz - cg_body_xyz, force_b)
    np.testing.assert_allclose(moment_b, expected, rtol=1e-6)


def test_zero_cg_offset_reproduces_unshifted_adb_moment():
    reference = AeroDBReference()
    coefficients = {"cfx": 0.0, "cfy": 0.0, "cfz": -1.0, "cmx": 0.1, "cmy": -0.2, "cmz": 0.3}

    _, unshifted = dimensional_body_force_moment(
        coefficients,
        reference.q,
        reference,
        shift_to_cg=False,
    )
    _, shifted = dimensional_body_force_moment(
        coefficients,
        reference.q,
        reference,
        shift_to_cg=True,
        cg_body_xyz=reference.moment_center_body_xyz,
    )

    np.testing.assert_allclose(shifted, unshifted)


def test_forward_and_aft_cg_shift_pitching_moment_sign():
    reference = AeroDBReference(moment_center_body_xyz=np.array([0.0, 0.0, 0.0]))
    coefficients = {"cfx": 0.0, "cfy": 0.0, "cfz": -1.0, "cmx": 0.0, "cmy": 0.0, "cmz": 0.0}

    _, moment_forward_cg = dimensional_body_force_moment(
        coefficients,
        reference.q,
        reference,
        cg_body_xyz=np.array([0.1, 0.0, 0.0]),
    )
    _, moment_aft_cg = dimensional_body_force_moment(
        coefficients,
        reference.q,
        reference,
        cg_body_xyz=np.array([-0.1, 0.0, 0.0]),
    )

    assert moment_forward_cg[1] < 0.0
    assert moment_aft_cg[1] > 0.0


def test_adb_metadata_overrides_stale_reference_fallback(tmp_path):
    adb_path = tmp_path / "adb.csv"
    _write_sample_adb(adb_path)
    stale_reference = AeroDBReference(
        rho=9.9,
        V_ref=9.9,
        S=9.9,
        b=9.9,
        c=9.9,
        moment_center_body_xyz=np.array([9.9, 9.9, 9.9]),
    )

    database = load_aerodynamic_database(adb_path, stale_reference)

    assert database.reference.rho == 1.09
    assert database.reference.V_ref == 25.0
    assert database.reference.S == 1.805
    assert database.reference.b == 3.95
    assert database.reference.c == 0.457
    np.testing.assert_allclose(database.reference.moment_center_body_xyz, [-0.15, 0.0, 0.15])


def test_aero_model_adb_passes_cg_body_xyz_to_moment_shift(monkeypatch):
    reference = AeroDBReference(moment_center_body_xyz=np.array([-0.15, 0.0, 0.15]))
    expected_cg = np.array([-0.20, 0.0, 0.12])
    captured = {}

    class FakeDatabase:
        def __init__(self, adb_reference):
            self.reference = adb_reference

        def coefficients_at(self, alpha_deg, beta_deg, *, clamp=False):
            return {
                "cfx": 0.0,
                "cfy": 0.0,
                "cfz": -1.0,
                "cmx": 0.0,
                "cmy": 0.0,
                "cmz": 0.0,
            }

    def fake_load_cached(path, fallback_reference=None):
        return FakeDatabase(reference)

    def fake_dimensional_body_force_moment(
        coefficients,
        qbar,
        adb_reference,
        *,
        shift_to_cg=True,
        cg_body_xyz=None,
    ):
        captured["reference"] = adb_reference
        captured["cg_body_xyz"] = np.asarray(cg_body_xyz)
        return np.zeros(3), np.zeros(3)

    monkeypatch.setattr(aero_module, "load_aerodynamic_database_cached", fake_load_cached)
    monkeypatch.setattr(aero_module, "dimensional_body_force_moment", fake_dimensional_body_force_moment)

    params = AircraftParams(cg_body_xyz=expected_cg)
    aero_module.aero_model_adb(
        np.array([25.0, 0.0, 0.0]),
        np.zeros(3),
        Control(de=0.0, da=0.0, dr=0.0, throttle=0.0),
        params,
    )

    assert captured["reference"] is reference
    np.testing.assert_allclose(captured["cg_body_xyz"], expected_cg)


def test_no_extrapolation_outside_alpha_beta_bounds(tmp_path):
    adb_path = tmp_path / "adb.csv"
    _write_sample_adb(adb_path)
    database = load_aerodynamic_database(adb_path)

    with pytest.raises(ValueError, match="outside alpha/beta bounds"):
        database.coefficients_at(2.0, 0.0)

    clamped = database.coefficients_at(2.0, 0.0, clamp=True)
    assert clamped["cfx"] == 3.0


def test_loads_real_pr7_adb_grid_and_exact_corner():
    database = load_aerodynamic_database(PR7_ADB_PATH)

    assert database.alpha_grid_deg.size == 181
    assert database.beta_grid_deg.size == 181
    assert database.metadata["Nondimit export"].startswith("rho=1.09")
    assert database.reference.rho == 1.09
    assert database.reference.V_ref == 25.0
    assert database.reference.q == 340.625
    assert database.reference.S == 1.805
    assert database.reference.b == 3.95
    assert database.reference.c == 0.457
    np.testing.assert_allclose(database.reference.moment_center_body_xyz, [-0.15, 0.0, 0.15])

    coefficients = database.coefficients_at(-90.0, -90.0)
    np.testing.assert_allclose(
        [coefficients["cfx"], coefficients["cfy"], coefficients["cfz"]],
        [0.0144, 0.2224, -0.0624],
        atol=5e-5,
    )
