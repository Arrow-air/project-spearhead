"""JSON export for structured simulation results."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TextIO

import numpy as np

from ..result import SimulationResult


def _trim_summary(result: SimulationResult) -> dict | None:
    trim = result.trim_result
    if trim is None:
        return None
    return {
        "success": bool(trim["success"]),
        "message": str(trim["message"]),
        "alpha_deg": float(trim["alpha_deg"]),
        "theta_deg": float(trim["theta_deg"]),
        "de_deg": float(trim["de_deg"]),
        "throttle": float(trim["throttle"]),
        "residuals": np.asarray(trim["residuals"], dtype=float).tolist(),
    }


def result_to_json_dict(result: SimulationResult) -> dict:
    """Return a JSON-serializable representation of a simulation result."""
    return {
        "name": result.config.name,
        "success": result.success,
        "message": result.message,
        "duration": result.config.duration,
        "initial_state": result.initial_state.tolist(),
        "trim": _trim_summary(result),
        "time": result.time.tolist(),
        "states": [result.state_at_index(index).as_dict() for index in range(result.time.size)],
    }


def export_json(result: SimulationResult, destination: str | Path | TextIO, *, indent: int = 2) -> None:
    """Write a simulation result to JSON."""
    payload = result_to_json_dict(result)
    should_close = False
    if hasattr(destination, "write"):
        stream = destination
    else:
        stream = Path(destination).open("w")
        should_close = True

    try:
        json.dump(payload, stream, indent=indent)
        stream.write("\n")
    finally:
        if should_close:
            stream.close()
