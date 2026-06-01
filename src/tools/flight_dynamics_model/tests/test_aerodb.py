from pathlib import Path

import numpy as np
import pytest

from spearhead.aerodb import (
    dimensional_body_force_moment,
    load_aerodynamic_database,
)
from spearhead.params import AeroDBReference


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
                "# adb_version,v1.1",
                "# moment_center_body_xyz,-0.15,0.0,0.0",
                "# notes,synthetic test table",
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

    assert database.metadata["adb_version"] == "v1.1"
    assert database.alpha_grid_deg.tolist() == [0.0, 1.0]
    assert database.beta_grid_deg.tolist() == [0.0, 1.0]


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


def test_moment_center_shift_from_adb_center_to_cg():
    reference = AeroDBReference()
    coefficients = {"cfx": 1.0, "cfy": 2.0, "cfz": 3.0, "cmx": 0.1, "cmy": 0.2, "cmz": 0.3}

    force_b, moment_b = dimensional_body_force_moment(
        coefficients,
        reference.q,
        reference,
        shift_to_cg=True,
    )

    moment_before_shift = np.array([0.1 * reference.qSb, 0.2 * reference.qSc, 0.3 * reference.qSb])
    expected = moment_before_shift + np.cross(reference.moment_center_body_xyz, force_b)
    np.testing.assert_allclose(moment_b, expected, rtol=1e-6)


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
    assert database.metadata["Nondimit export"].startswith("rho=1.225")

    coefficients = database.coefficients_at(-90.0, -90.0)
    np.testing.assert_allclose(
        [coefficients["cfx"], coefficients["cfy"], coefficients["cfz"]],
        [0.0128, 0.1979, -0.0555],
        atol=5e-5,
    )
