"""Scenario-loading adapters for GUI code."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ...params import AircraftParams
from ...scenarios import load_scenario
from ...config import SimulationConfig


def load_scenario_config(path: str | Path, params: AircraftParams | None = None) -> SimulationConfig:
    """Load a YAML scenario with the shared scenario loader."""
    return load_scenario(path, params=params)


def scenario_summary(config: SimulationConfig) -> dict[str, Any]:
    """Return a compact, GUI-friendly summary of a simulation config."""
    return {
        "name": config.name,
        "duration": config.duration,
        "dt": config.dt,
        "n_points": config.n_points,
        "start_from_trim": config.start_from_trim,
        "trim_target_speed": config.trim_target_speed,
        "trim_altitude": config.trim_altitude,
        "input_commands": [
            {
                "surface": command.surface,
                "kind": command.kind,
                "amplitude": command.amplitude,
                "start_time": command.start_time,
                "duration": command.duration,
            }
            for command in config.input_commands
        ],
        "cg_body_xyz": config.params.cg_body_xyz.tolist(),
    }

