"""High-level stability analysis built on the shared simulation backend."""

from __future__ import annotations

import numpy as np

from ..config import SimulationConfig
from ..force_moment import compute_force_moment_breakdown
from ..simulation import run_simulation
from .linearize import linearize
from .modes import classify_lateral_directional_modes, classify_longitudinal_modes, modal_analysis
from .types import LinearizationConfig, StabilityResult, StaticDerivativeResult


def _state_with_alpha_beta(trim_state: np.ndarray, alpha: float, beta: float, speed: float) -> np.ndarray:
    state = trim_state.copy()
    state[3] = speed * np.cos(alpha) * np.cos(beta)
    state[4] = speed * np.sin(beta)
    state[5] = speed * np.sin(alpha) * np.cos(beta)
    return state


def _aero_coefficients(config: SimulationConfig, state: np.ndarray, control) -> dict[str, float]:
    params = config.params
    breakdown = compute_force_moment_breakdown(0.0, state, params, control_override=control)
    force_b = (
        breakdown["wing"]["force_b"]
        + breakdown["tail"]["force_b"]
        + breakdown["fuselage"]["force_b"]
    )
    moment_b = (
        breakdown["wing"]["moment_b"]
        + breakdown["tail"]["moment_b"]
        + breakdown["fuselage"]["moment_b"]
    )
    speed = max(float(np.linalg.norm(state[3:6])), params.Vmin)
    qbar = 0.5 * params.rho * speed**2
    qS = qbar * params.Sref
    return {
        "Cfx": float(force_b[0] / qS),
        "Cfy": float(force_b[1] / qS),
        "Cfz": float(force_b[2] / qS),
        "CL": float(-force_b[2] / qS),
        "Cl": float(moment_b[0] / (qS * params.bref)),
        "Cm": float(moment_b[1] / (qS * params.cref)),
        "Cn": float(moment_b[2] / (qS * params.bref)),
    }


def extract_static_derivatives(
    config: SimulationConfig,
    trim_state: np.ndarray,
    trim_control,
    *,
    alpha_step_rad: float = np.deg2rad(1.0),
    beta_step_rad: float = np.deg2rad(1.0),
) -> StaticDerivativeResult:
    """Estimate static derivatives with central finite differences."""
    if alpha_step_rad <= 0.0 or beta_step_rad <= 0.0:
        raise ValueError("Static derivative perturbations must be positive")

    velocity = trim_state[3:6]
    speed = max(float(np.linalg.norm(velocity)), config.params.Vmin)
    alpha = float(np.arctan2(velocity[2], velocity[0]))
    beta = float(np.arcsin(np.clip(velocity[1] / speed, -1.0, 1.0)))

    alpha_plus = _aero_coefficients(
        config,
        _state_with_alpha_beta(trim_state, alpha + alpha_step_rad, beta, speed),
        trim_control,
    )
    alpha_minus = _aero_coefficients(
        config,
        _state_with_alpha_beta(trim_state, alpha - alpha_step_rad, beta, speed),
        trim_control,
    )
    beta_plus = _aero_coefficients(
        config,
        _state_with_alpha_beta(trim_state, alpha, beta + beta_step_rad, speed),
        trim_control,
    )
    beta_minus = _aero_coefficients(
        config,
        _state_with_alpha_beta(trim_state, alpha, beta - beta_step_rad, speed),
        trim_control,
    )

    alpha_scale = 2.0 * alpha_step_rad
    beta_scale = 2.0 * beta_step_rad
    derivatives = {
        "CL_alpha_per_rad": (alpha_plus["CL"] - alpha_minus["CL"]) / alpha_scale,
        "Cm_alpha_per_rad": (alpha_plus["Cm"] - alpha_minus["Cm"]) / alpha_scale,
        "Cfx_alpha_per_rad": (alpha_plus["Cfx"] - alpha_minus["Cfx"]) / alpha_scale,
        "Cfz_alpha_per_rad": (alpha_plus["Cfz"] - alpha_minus["Cfz"]) / alpha_scale,
        "CY_beta_per_rad": (beta_plus["Cfy"] - beta_minus["Cfy"]) / beta_scale,
        "Cl_beta_per_rad": (beta_plus["Cl"] - beta_minus["Cl"]) / beta_scale,
        "Cn_beta_per_rad": (beta_plus["Cn"] - beta_minus["Cn"]) / beta_scale,
    }
    return StaticDerivativeResult(
        derivatives={key: float(value) for key, value in derivatives.items()},
        source="force_moment_breakdown aerodynamic components, central finite differences",
        alpha_deg=float(np.rad2deg(alpha)),
        beta_deg=float(np.rad2deg(beta)),
    )


def analyze_stability(
    config: SimulationConfig,
    linearization_config: LinearizationConfig | None = None,
) -> StabilityResult:
    """Run trim, linearization, static derivatives, and modal analysis."""
    sim_result = run_simulation(config)
    if sim_result.trim_result is None:
        raise ValueError("Stability analysis requires config.start_from_trim=True")
    if not bool(sim_result.trim_result["success"]):
        raise ValueError(f"Trim failed: {sim_result.trim_result['message']}")

    linearization = linearize(config, settings=linearization_config)
    longitudinal_labels = classify_longitudinal_modes(np.linalg.eigvals(linearization.A_longitudinal))
    lateral_labels = classify_lateral_directional_modes(
        np.linalg.eigvals(linearization.A_lateral_directional)
    )
    longitudinal = modal_analysis(
        linearization.A_longitudinal,
        subsystem="longitudinal",
        labels=longitudinal_labels,
    )
    lateral_directional = modal_analysis(
        linearization.A_lateral_directional,
        subsystem="lateral_directional",
        labels=lateral_labels,
    )
    static_derivatives = extract_static_derivatives(
        config,
        linearization.trim_state,
        linearization.trim_control,
    )
    return StabilityResult(
        config=config,
        linearization=linearization,
        longitudinal=longitudinal,
        lateral_directional=lateral_directional,
        static_derivatives=static_derivatives,
        trim_result=sim_result.trim_result,
    )
