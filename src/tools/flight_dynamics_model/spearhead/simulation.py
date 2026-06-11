"""Simulation entry points."""

import numpy as np
from scipy.integrate import solve_ivp

from .controls import Control
from .dynamics import dynamics
from .params import AircraftParams
from .trim import trim_fixed_wing_longitudinal


def default_initial_state() -> np.ndarray:
    """Return the default 12-state initial condition."""
    return np.array(
        [
            0.0,
            0.0,
            -100.0,
            22.0,
            0.0,
            0.0,
            np.deg2rad(0.0),
            np.deg2rad(2.0),
            np.deg2rad(0.0),
            0.0,
            0.0,
            0.0,
        ]
    )


def run_open_loop(
    t_final: float = 20.0,
    n_points: int = 1001,
    params: AircraftParams | None = None,
    initial_state: np.ndarray | None = None,
    control_override: Control | dict[str, float] | None = None,
    start_from_trim: bool = True,
    trim_target_speed: float = 22.0,
    trim_altitude: float = 100.0,
):
    """Run an open-loop nonlinear simulation.

    By default this starts from the preliminary longitudinal trim condition and
    holds the trimmed control input constant. Pass ``start_from_trim=False`` to
    use the original arbitrary default initial state and scheduled controls.
    """
    if params is None:
        params = AircraftParams()

    trim_result = None
    if initial_state is None and start_from_trim:
        trim_result = trim_fixed_wing_longitudinal(
            params,
            target_speed=trim_target_speed,
            altitude=trim_altitude,
        )
        initial_state = trim_result["x_trim"]
        if control_override is None:
            control_override = trim_result["control_trim"]

    y0 = default_initial_state() if initial_state is None else initial_state
    t_eval = np.linspace(0.0, t_final, n_points)
    sol = solve_ivp(
        fun=lambda t, x: dynamics(t, x, params, control_override=control_override),
        t_span=(0.0, t_final),
        y0=y0,
        t_eval=t_eval,
        rtol=1e-7,
        atol=1e-9,
    )
    sol.initial_state = y0
    sol.control_override = control_override
    sol.trim_result = trim_result
    return sol
