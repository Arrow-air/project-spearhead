"""Simulation adapters for GUI code."""

from __future__ import annotations

from typing import Any

from ...config import SimulationConfig
from ...result import SimulationResult
from ...simulation import run_simulation


def build_simulation_config(**overrides: Any) -> SimulationConfig:
    """Build a backend ``SimulationConfig`` from GUI form values."""
    return SimulationConfig(**overrides)


def run_simulation_for_gui(config: SimulationConfig) -> SimulationResult:
    """Run the shared simulation backend for a GUI request."""
    return run_simulation(config)


def _trim_summary(result: SimulationResult) -> dict[str, Any] | None:
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
        "residuals": [float(value) for value in trim["residuals"]],
    }


def simulation_result_summary(result: SimulationResult) -> dict[str, Any]:
    """Return a compact summary suitable for GUI status panels."""
    final_state = result.final_state.as_dict()
    return {
        "name": result.config.name,
        "success": result.success,
        "message": result.message,
        "duration": result.config.duration,
        "samples": int(result.time.size),
        "start_time": float(result.time[0]),
        "end_time": float(result.time[-1]),
        "trim": _trim_summary(result),
        "final_state": final_state,
    }


def simulation_time_history_rows(result: SimulationResult) -> list[dict[str, float]]:
    """Return one row per simulation sample using existing state labels."""
    return [result.state_at_index(index).as_dict() for index in range(result.time.size)]
