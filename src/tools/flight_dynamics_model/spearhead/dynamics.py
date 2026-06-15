"""Nonlinear 6-DOF rigid-body dynamics."""

import numpy as np

from .controls import ControlInput
from .force_moment import compute_force_moment_breakdown
from .params import AircraftParams
from .rotations import body_rates_to_euler_rates, euler_body_to_ned


def dynamics(
    t: float,
    x: np.ndarray,
    params: AircraftParams,
    control_override: ControlInput | None = None,
) -> np.ndarray:
    """Return the state derivative for the 12-state nonlinear aircraft model."""
    v_b = x[3:6]
    phi, theta, psi = x[6:9]
    omega_b = x[9:12]

    R_body_to_ned = euler_body_to_ned(phi, theta, psi)
    breakdown = compute_force_moment_breakdown(t, x, params, control_override=control_override)
    F_total_b = breakdown["total"]["force_b"]
    M_total_b = breakdown["total"]["moment_b"]

    v_dot_b = F_total_b / params.mass - np.cross(omega_b, v_b)
    omega_dot_b = np.linalg.solve(
        params.inertia,
        M_total_b - np.cross(omega_b, params.inertia @ omega_b),
    )
    pos_dot_ned = R_body_to_ned @ v_b
    euler_dot = body_rates_to_euler_rates(phi, theta, omega_b)

    return np.concatenate([pos_dot_ned, v_dot_b, euler_dot, omega_dot_b])
