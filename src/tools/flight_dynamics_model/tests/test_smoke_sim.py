import numpy as np

from spearhead import dynamics as dynamics_module
from spearhead.force_moment import compute_force_moment_breakdown
from spearhead.params import AircraftParams
from spearhead.simulation import run_open_loop


def test_open_loop_simulation_runs_and_returns_finite_states():
    sol = run_open_loop(t_final=1.0, n_points=11)
    assert sol.success
    assert sol.y.shape == (12, 11)
    assert np.all(np.isfinite(sol.y))


def test_open_loop_defaults_to_trimmed_equilibrium():
    params = AircraftParams()
    sol = run_open_loop(t_final=0.1, n_points=2, params=params)
    breakdown = compute_force_moment_breakdown(
        0.0,
        sol.y[:, 0],
        params,
        control_override=sol.control_override,
    )

    np.testing.assert_allclose(sol.y[:, 0], sol.trim_result["x_trim"])
    assert sol.control_override == sol.trim_result["control_trim"]
    np.testing.assert_allclose(breakdown["total"]["force_b"][[0, 2]], 0.0, atol=1e-7)
    np.testing.assert_allclose(breakdown["total"]["moment_b"][1], 0.0, atol=1e-7)


def test_dynamics_uses_shared_total_force_and_moment(monkeypatch):
    params = AircraftParams()
    x = np.zeros(12)
    expected_force_b = np.array([params.mass, 2.0 * params.mass, 3.0 * params.mass])
    expected_moment_b = np.array([params.inertia[0, 0], 2.0 * params.inertia[1, 1], 3.0 * params.inertia[2, 2]])

    def fake_force_moment_breakdown(t, state, model_params, control_override=None):
        return {
            "total": {
                "force_b": expected_force_b,
                "moment_b": expected_moment_b,
            }
        }

    monkeypatch.setattr(
        dynamics_module,
        "compute_force_moment_breakdown",
        fake_force_moment_breakdown,
    )

    x_dot = dynamics_module.dynamics(0.0, x, params)

    np.testing.assert_allclose(x_dot[3:6], expected_force_b / params.mass)
    np.testing.assert_allclose(
        x_dot[9:12],
        np.linalg.solve(params.inertia, expected_moment_b),
    )
