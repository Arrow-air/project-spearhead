"""JSON and CSV exports for stability results."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, TextIO

import numpy as np

from .types import CGSweepResult, ModalAnalysisResult, StabilityResult


def _json_default(value: Any):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _mode_rows(result: StabilityResult, case: str | None = None) -> list[dict[str, Any]]:
    rows = []
    for analysis in (result.longitudinal, result.lateral_directional):
        for mode in analysis.modes:
            row = mode.as_dict()
            row["subsystem"] = analysis.subsystem
            if case is not None:
                row["case"] = case
            rows.append(row)
    return rows


def stability_to_dict(result: StabilityResult) -> dict[str, Any]:
    """Return a JSON-serializable stability result summary."""
    return {
        "name": result.config.name,
        "warning": result.warning,
        "cg_body_xyz": result.config.params.cg_body_xyz,
        "trim": {
            "success": bool(result.trim_result["success"]),
            "message": str(result.trim_result["message"]),
            "alpha_deg": float(result.trim_result["alpha_deg"]),
            "theta_deg": float(result.trim_result["theta_deg"]),
            "de_deg": float(result.trim_result["de_deg"]),
            "throttle": float(result.trim_result["throttle"]),
            "residuals": np.asarray(result.trim_result["residuals"], dtype=float),
        },
        "static_derivatives": result.static_derivatives.derivatives,
        "static_derivative_source": result.static_derivatives.source,
        "linearization": {
            "A": result.linearization.A,
            "B": result.linearization.B,
            "A_longitudinal": result.linearization.A_longitudinal,
            "B_longitudinal": result.linearization.B_longitudinal,
            "A_lateral_directional": result.linearization.A_lateral_directional,
            "B_lateral_directional": result.linearization.B_lateral_directional,
            "state_labels": result.linearization.state_labels,
            "control_labels": result.linearization.control_labels,
        },
        "modes": _mode_rows(result),
    }


def cg_sweep_to_dict(result: CGSweepResult) -> dict[str, Any]:
    """Return a JSON-serializable CG sweep summary."""
    points = []
    for index, point in enumerate(result.points):
        payload: dict[str, Any] = {
            "case": f"cg_{index}",
            "cg_body_xyz": point.cg_body_xyz,
            "success": point.success,
            "error": point.error,
        }
        if point.result is not None:
            payload["stability"] = stability_to_dict(point.result)
        points.append(payload)
    return {"warning": result.warning, "points": points}


def export_json(payload: StabilityResult | CGSweepResult, destination: str | Path | TextIO) -> None:
    """Write a stability result or CG sweep to JSON."""
    data = stability_to_dict(payload) if isinstance(payload, StabilityResult) else cg_sweep_to_dict(payload)
    should_close = False
    if hasattr(destination, "write"):
        stream = destination
    else:
        stream = Path(destination).open("w")
        should_close = True
    try:
        json.dump(data, stream, indent=2, default=_json_default)
        stream.write("\n")
    finally:
        if should_close:
            stream.close()


def export_modes_csv(payload: StabilityResult | CGSweepResult, destination: str | Path | TextIO) -> None:
    """Write modal eigenvalue rows to CSV."""
    rows: list[dict[str, Any]] = []
    if isinstance(payload, StabilityResult):
        rows.extend(_mode_rows(payload))
    else:
        for index, point in enumerate(payload.points):
            if point.result is not None:
                rows.extend(_mode_rows(point.result, case=f"cg_{index}"))

    fieldnames = [
        "case",
        "subsystem",
        "mode",
        "eigenvalue_real",
        "eigenvalue_imag",
        "natural_frequency_rad_s",
        "damping_ratio",
        "oscillation_period_s",
        "time_to_half_s",
        "time_to_double_s",
        "classification",
    ]
    should_close = False
    if hasattr(destination, "write"):
        stream = destination
    else:
        stream = Path(destination).open("w", newline="")
        should_close = True
    try:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    finally:
        if should_close:
            stream.close()
