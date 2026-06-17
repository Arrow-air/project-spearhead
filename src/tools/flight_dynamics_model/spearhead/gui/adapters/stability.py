"""Stability-analysis adapters for GUI code."""

from __future__ import annotations

from typing import Any

import numpy as np

from ...config import SimulationConfig
from ...stability import analyze_stability, run_cg_sweep
from ...stability.types import CGSweepResult, LinearizationConfig, StabilityResult


def build_linearization_config(
    *,
    state_perturbation: float = 1e-5,
    control_perturbation: float = 1e-5,
) -> LinearizationConfig:
    """Build backend finite-difference settings from GUI form values."""
    return LinearizationConfig(
        state_perturbation=state_perturbation,
        control_perturbation=control_perturbation,
    )


def run_stability_for_gui(
    config: SimulationConfig,
    linearization_config: LinearizationConfig | None = None,
) -> StabilityResult:
    """Run the shared stability backend for a GUI request."""
    return analyze_stability(config, linearization_config=linearization_config)


def run_cg_sweep_for_gui(
    config: SimulationConfig,
    cg_positions: list[np.ndarray] | tuple[np.ndarray, ...],
    linearization_config: LinearizationConfig | None = None,
) -> CGSweepResult:
    """Run the shared CG sweep backend for a GUI request."""
    return run_cg_sweep(config, cg_positions, linearization_config=linearization_config)


def stability_mode_table_rows(result: StabilityResult, case: str | None = None) -> list[dict[str, Any]]:
    """Return modal metrics as GUI table rows."""
    rows: list[dict[str, Any]] = []
    for analysis in (result.longitudinal, result.lateral_directional):
        for mode in analysis.modes:
            row = mode.as_dict()
            row["subsystem"] = analysis.subsystem
            row["eigenvalue_label"] = f"{mode.eigenvalue.real:+.6g}{mode.eigenvalue.imag:+.6g}j"
            if case is not None:
                row["case"] = case
            rows.append(row)
    return rows


def static_derivative_table_rows(result: StabilityResult) -> list[dict[str, Any]]:
    """Return static aerodynamic derivatives as GUI table rows."""
    return [
        {
            "name": name,
            "value": float(value),
            "source": result.static_derivatives.source,
            "alpha_deg": result.static_derivatives.alpha_deg,
            "beta_deg": result.static_derivatives.beta_deg,
        }
        for name, value in sorted(result.static_derivatives.derivatives.items())
    ]


def cg_sweep_summary_rows(result: CGSweepResult) -> list[dict[str, Any]]:
    """Return one compact GUI table row per CG sweep point."""
    rows: list[dict[str, Any]] = []
    for index, point in enumerate(result.points):
        row: dict[str, Any] = {
            "case": f"cg_{index}",
            "cg_x": float(point.cg_body_xyz[0]),
            "cg_y": float(point.cg_body_xyz[1]),
            "cg_z": float(point.cg_body_xyz[2]),
            "success": point.success,
            "error": point.error,
            "longitudinal_modes": 0,
            "lateral_directional_modes": 0,
            "trim_alpha_deg": None,
            "trim_theta_deg": None,
        }
        if point.result is not None:
            trim = point.result.trim_result
            row.update(
                {
                    "longitudinal_modes": len(point.result.longitudinal.modes),
                    "lateral_directional_modes": len(point.result.lateral_directional.modes),
                    "trim_alpha_deg": float(trim["alpha_deg"]),
                    "trim_theta_deg": float(trim["theta_deg"]),
                }
            )
        rows.append(row)
    return rows
