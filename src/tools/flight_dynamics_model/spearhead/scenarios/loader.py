"""YAML scenario loading for shared simulation configs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from ..config import ControlInputCommand, SimulationConfig
from ..params import AircraftParams
from ..aircraft import load_aircraft_config, select_aerodynamics
from ..controls import as_control
from ..propulsion import hover_collective, hover_trim_control


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


def load_scenario(path: str | Path, params: AircraftParams | None = None, *, aerodynamics: str | None = None) -> SimulationConfig:
    """Load a YAML scenario into a ``SimulationConfig``."""
    scenario_path = Path(path)
    with scenario_path.open() as stream:
        data = yaml.safe_load(stream) or {}

    aircraft_case = data.get("aircraft_config", data.get("aircraft_case", "nominal"))
    if params is None:
        if aircraft_case == "nominal":
            params = load_aircraft_config("reference_design")
        else:
            profile = Path(aircraft_case)
            if profile.suffix in {".yaml", ".yml"} and not profile.is_absolute():
                profile = scenario_path.parent / profile
            params = load_aircraft_config(profile)
        if data.get("aerodynamics") is not None:
            params = select_aerodynamics(params, str(data["aerodynamics"]))
    if aerodynamics is not None:
        params = select_aerodynamics(params, aerodynamics)

    trim = data.get("trim") or {}
    simulation = data.get("simulation") or {}
    inputs = data.get("inputs") or {}
    commands = tuple(
        _command_from_yaml(surface, command_data or {})
        for surface, command_data in inputs.items()
    )

    dt = simulation.get("dt")
    n_points = simulation.get("n_points")
    control = data.get("control")
    if control is not None:
        control = dict(control)
        if control.get("lift_throttle") == "balanced_hover":
            balanced = hover_trim_control(params)
            control["lift_throttle"] = balanced.lift_throttle
            control["lift_motor_throttles"] = balanced.lift_motor_throttles
        elif control.get("lift_throttle") == "hover":
            control["lift_throttle"] = hover_collective(params)
        control = as_control({"de": 0.0, "da": 0.0, "dr": 0.0, "throttle": 0.0, **control})
    return SimulationConfig(
        name=str(data.get("name", scenario_path.stem)),
        duration=float(simulation.get("duration", 20.0)),
        dt=float(dt) if dt is not None else None,
        n_points=int(n_points) if n_points is not None else None if dt is not None else 1001,
        params=params,
        initial_state=data.get("initial_state"),
        control_override=control,
        input_commands=commands,
        start_from_trim=bool(trim.get("enabled", True)),
        trim_target_speed=float(trim.get("airspeed", 22.0)),
        trim_altitude=float(trim.get("altitude", 100.0)),
    )
