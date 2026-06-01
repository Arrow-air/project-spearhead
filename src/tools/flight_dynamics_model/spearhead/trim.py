"""Simple preliminary fixed-wing trim utilities."""

import numpy as np
from scipy.optimize import least_squares

from .controls import Control
from .force_moment import compute_force_moment_breakdown
from .params import AircraftParams


def _state_from_trim_variables(
    alpha: float,
    theta: float,
    target_speed: float,
    altitude: float,
) -> np.ndarray:
    return np.array(
        [
            0.0,
            0.0,
            -altitude,
            target_speed * np.cos(alpha),
            0.0,
            target_speed * np.sin(alpha),
            0.0,
            theta,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )


def _parse_initial_guess(initial_guess) -> np.ndarray:
    if initial_guess is None:
        return np.array([np.deg2rad(2.0), np.deg2rad(2.0), 0.0, 0.45])
    if isinstance(initial_guess, dict):
        return np.array(
            [
                float(initial_guess["alpha"]),
                float(initial_guess["theta"]),
                float(initial_guess["de"]),
                float(initial_guess["throttle"]),
            ]
        )
    return np.asarray(initial_guess, dtype=float)


def trim_fixed_wing_longitudinal(
    params: AircraftParams,
    target_speed: float = 22.0,
    altitude: float = 100.0,
    initial_guess=None,
) -> dict:
    """Trim the placeholder model for simple wings-level longitudinal flight.

    This is a preliminary solver for the current placeholder aero/propulsion
    model. It is not a validated Spearhead trim result.
    """
    x0 = _parse_initial_guess(initial_guess)
    lower_bounds = np.array([np.deg2rad(-10.0), np.deg2rad(-10.0), np.deg2rad(-25.0), 0.0])
    upper_bounds = np.array([np.deg2rad(20.0), np.deg2rad(20.0), np.deg2rad(25.0), 1.0])
    x0 = np.clip(x0, lower_bounds, upper_bounds)

    def residual(unknowns: np.ndarray) -> np.ndarray:
        alpha, theta, de, throttle = unknowns
        x_trim = _state_from_trim_variables(alpha, theta, target_speed, altitude)
        control_trim = Control(de=de, da=0.0, dr=0.0, throttle=throttle)
        breakdown = compute_force_moment_breakdown(
            0.0,
            x_trim,
            params,
            control_override=control_trim,
        )
        F_total_b = breakdown["total"]["force_b"]
        M_total_b = breakdown["total"]["moment_b"]
        return np.array(
            [
                F_total_b[0] / params.mass,
                F_total_b[2] / params.mass,
                M_total_b[1] / params.inertia[1, 1],
                theta - alpha,
            ]
        )

    result = least_squares(
        residual,
        x0,
        bounds=(lower_bounds, upper_bounds),
        xtol=1e-10,
        ftol=1e-10,
        gtol=1e-10,
    )

    alpha, theta, de, throttle = result.x
    x_trim = _state_from_trim_variables(alpha, theta, target_speed, altitude)
    control_trim = Control(de=de, da=0.0, dr=0.0, throttle=throttle)
    breakdown = compute_force_moment_breakdown(
        0.0,
        x_trim,
        params,
        control_override=control_trim,
    )

    return {
        "success": bool(result.success),
        "message": result.message,
        "alpha_rad": alpha,
        "alpha_deg": np.rad2deg(alpha),
        "theta_rad": theta,
        "theta_deg": np.rad2deg(theta),
        "de_rad": de,
        "de_deg": np.rad2deg(de),
        "throttle": throttle,
        "state_vector": x_trim,
        "x_trim": x_trim,
        "control_trim": control_trim,
        "residuals": residual(result.x),
        "force_moment_breakdown": breakdown,
    }

