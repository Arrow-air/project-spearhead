"""Typed configuration for Spearhead simulations."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .controls import ControlInput
from .params import AircraftParams


@dataclass(frozen=True)
class ControlInputCommand:
    """A simple time-domain perturbation applied to one control channel."""

    surface: str
    kind: str
    amplitude: float
    start_time: float
    duration: float = 2.0

    def value_at(self, t: float) -> float:
        """Return the command increment at time ``t``."""
        if self.kind == "doublet":
            half_duration = 0.5 * self.duration
            if self.start_time <= t < self.start_time + half_duration:
                return self.amplitude
            if self.start_time + half_duration <= t < self.start_time + self.duration:
                return -self.amplitude
            return 0.0
        if self.kind == "step":
            return self.amplitude if t >= self.start_time else 0.0
        raise ValueError(f"Unsupported input command type: {self.kind}")


@dataclass(frozen=True)
class SimulationConfig:
    """Configuration consumed by the shared simulation backend."""

    name: str = "open_loop"
    duration: float = 20.0
    dt: float | None = None
    n_points: int | None = 1001
    params: AircraftParams = field(default_factory=AircraftParams)
    initial_state: np.ndarray | None = None
    control_override: ControlInput | None = None
    input_commands: tuple[ControlInputCommand, ...] = ()
    start_from_trim: bool = True
    trim_target_speed: float = 22.0
    trim_altitude: float = 100.0
    rtol: float = 1e-7
    atol: float = 1e-9

    def __post_init__(self) -> None:
        if self.duration <= 0.0:
            raise ValueError("Simulation duration must be positive")
        if self.dt is not None and self.dt <= 0.0:
            raise ValueError("Simulation dt must be positive when provided")
        if self.n_points is not None and self.n_points < 2:
            raise ValueError("Simulation n_points must be at least 2 when provided")
        if self.initial_state is not None:
            state = np.asarray(self.initial_state, dtype=float)
            if state.shape != (12,):
                raise ValueError("Simulation initial_state must have shape (12,)")
            object.__setattr__(self, "initial_state", state)
