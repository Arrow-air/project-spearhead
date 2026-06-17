"""YAML scenario loading for shared simulation configs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from ..config import ControlInputCommand, SimulationConfig
from ..params import AircraftParams


_ANGLE_SURFACES = {"elevator", "de", "aileron", "da", "rudder", "ruddervator", "dr"}


def _command_from_yaml(surface: str, data: dict) -> ControlInputCommand:
    kind = str(data.get("type", "step"))
    if surface in _ANGLE_SURFACES:
        if "amplitude_deg" in data:
            amplitude = np.deg2rad(float(data["amplitude_deg"]))
        else:
            amplitude = float(data.get("amplitude_rad", 0.0))
    else:
        amplitude = float(data.get("amplitude", data.get("amplitude_fraction", 0.0)))

    return ControlInputCommand(
        surface=surface,
        kind=kind,
        amplitude=amplitude,
        start_time=float(data.get("start_time", 0.0)),
        duration=float(data.get("duration", data.get("duration_s", 2.0))),
    )


def load_scenario(path: str | Path, params: AircraftParams | None = None) -> SimulationConfig:
    """Load a YAML scenario into a ``SimulationConfig``."""
    scenario_path = Path(path)
    with scenario_path.open() as stream:
        data = yaml.safe_load(stream) or {}

    aircraft_case = data.get("aircraft_case", "nominal")
    if aircraft_case != "nominal":
        raise ValueError(f"Unsupported aircraft_case: {aircraft_case}")

    trim = data.get("trim") or {}
    simulation = data.get("simulation") or {}
    inputs = data.get("inputs") or {}
    commands = tuple(
        _command_from_yaml(surface, command_data or {})
        for surface, command_data in inputs.items()
    )

    dt = simulation.get("dt")
    n_points = simulation.get("n_points")
    return SimulationConfig(
        name=str(data.get("name", scenario_path.stem)),
        duration=float(simulation.get("duration", 20.0)),
        dt=float(dt) if dt is not None else None,
        n_points=int(n_points) if n_points is not None else None if dt is not None else 1001,
        params=AircraftParams() if params is None else params,
        input_commands=commands,
        start_from_trim=bool(trim.get("enabled", True)),
        trim_target_speed=float(trim.get("airspeed", 22.0)),
        trim_altitude=float(trim.get("altitude", 100.0)),
    )
