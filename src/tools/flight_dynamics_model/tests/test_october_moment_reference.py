"""A recentered table must preserve the wrench at every selected CG."""

import csv
from pathlib import Path

import numpy as np
import pytest

from spearhead.aerodb import dimensional_body_force_moment, load_aerodynamic_database_cached

ROOT = Path(__file__).resolve().parents[1] / "data/aerodynamics/tail_update_2026_10_07"


@pytest.mark.parametrize("cg", [(-.15, 0, .15), (.01, 0, .15), (.025, 0, .15), (.05, .02, -.03)])
@pytest.mark.parametrize("speed,alpha,beta", [(10, 0, 0), (25, 5, 0), (32, 30, -5)])
def test_origin_recovery_and_runtime_cg_shift_equal_direct_export_to_cg(cg, speed, alpha, beta):
    with (ROOT / "sources/database.csv").open() as stream:
        reader = csv.DictReader(line for line in stream if not line.startswith("#"))
        row = next(r for r in reader if float(r["Alpha(deg)"]) == alpha and float(r["Beta(deg)"]) == beta)
    force_export = np.array([float(row[k]) for k in ("Fx", "Fy", "Fz")])
    moment_export = np.array([float(row[k]) for k in ("Mx", "My", "Mz")])
    speed_scale = (speed / 25) ** 2
    expected_force = speed_scale * force_export
    expected_moment = speed_scale * (
        moment_export + np.cross(np.array([-.1, 0, .15]) - cg, force_export)
    )
    database = load_aerodynamic_database_cached(ROOT / "database.csv")
    actual_force, actual_moment = dimensional_body_force_moment(
        database.coefficients_at(alpha, beta), .5 * 1.021 * speed ** 2,
        database.reference, cg_body_xyz=np.array(cg),
    )
    np.testing.assert_allclose(actual_force, expected_force, atol=1e-6, rtol=0)
    np.testing.assert_allclose(actual_moment, expected_moment, atol=1e-6, rtol=0)
