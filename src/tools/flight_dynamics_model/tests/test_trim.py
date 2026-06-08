import numpy as np

from spearhead.controls import Control
from spearhead.force_moment import compute_force_moment_breakdown
from spearhead.params import AircraftParams
from spearhead.simulation import default_initial_state, run_open_loop
from spearhead.trim import trim_fixed_wing_longitudinal


def _longitudinal_residual_norm(breakdown, params):
    force_b = breakdown["total"]["force_b"]
    moment_b = breakdown["total"]["moment_b"]
    residuals = np.array(
        [
            force_b[0] / params.mass,
            force_b[2] / params.mass,
            moment_b[1] / params.inertia[1, 1],
        ]
    )
    return np.linalg.norm(residuals)


def test_trim_returns_finite_values_in_reasonable_bounds():
    params = AircraftParams()
    trim = trim_fixed_wing_longitudinal(params)

    assert trim["success"]
    for key in ("alpha_rad", "theta_rad", "de_rad", "throttle"):
        assert np.isfinite(trim[key])

    assert np.deg2rad(-10.0) <= trim["alpha_rad"] <= np.deg2rad(20.0)
    assert np.deg2rad(-10.0) <= trim["theta_rad"] <= np.deg2rad(20.0)
    assert np.deg2rad(-25.0) <= trim["de_rad"] <= np.deg2rad(25.0)
    assert 0.0 <= trim["throttle"] <= 1.0


def test_trim_residuals_are_smaller_than_default_condition():
    params = AircraftParams()
    trim = trim_fixed_wing_longitudinal(params)

    before = compute_force_moment_breakdown(
        0.0,
        default_initial_state(),
        params,
        control_override=Control(de=0.0, da=0.0, dr=0.0, throttle=0.45),
    )
    before_norm = _longitudinal_residual_norm(before, params)
    after_norm = np.linalg.norm(trim["residuals"][:3])

    assert after_norm < before_norm


def test_trimmed_simulation_runs_without_nans():
    params = AircraftParams()
    trim = trim_fixed_wing_longitudinal(params)

    sol = run_open_loop(
        t_final=1.0,
        n_points=11,
        params=params,
        initial_state=trim["x_trim"],
        control_override=trim["control_trim"],
    )

    assert sol.success
    assert np.all(np.isfinite(sol.y))

