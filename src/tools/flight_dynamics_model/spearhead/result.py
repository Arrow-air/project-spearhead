"""Simulation result containers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .config import SimulationConfig
from .state import SimulationState


@dataclass(frozen=True)
class SimulationResult:
    """Structured result returned by the shared simulation API."""

    config: SimulationConfig
    time: np.ndarray
    states: np.ndarray
    success: bool
    message: str
    initial_state: np.ndarray
    control_override: Any = None
    trim_result: dict[str, Any] | None = None
    raw_solution: Any = None

    def __post_init__(self) -> None:
        time = np.asarray(self.time, dtype=float)
        states = np.asarray(self.states, dtype=float)
        initial_state = np.asarray(self.initial_state, dtype=float)
        if time.ndim != 1:
            raise ValueError("SimulationResult time must be one-dimensional")
        if states.shape != (time.size, 12):
            raise ValueError("SimulationResult states must have shape (N, 12)")
        if initial_state.shape != (12,):
            raise ValueError("SimulationResult initial_state must have shape (12,)")
        object.__setattr__(self, "time", time)
        object.__setattr__(self, "states", states)
        object.__setattr__(self, "initial_state", initial_state)

    @property
    def final_state(self) -> SimulationState:
        """Return the final timestamped state."""
        return SimulationState(float(self.time[-1]), self.states[-1])

    def state_at_index(self, index: int) -> SimulationState:
        """Return a timestamped state by sample index."""
        return SimulationState(float(self.time[index]), self.states[index])

    @classmethod
    def from_solve_ivp(
        cls,
        *,
        config: SimulationConfig,
        solution,
        initial_state: np.ndarray,
        control_override=None,
        trim_result: dict[str, Any] | None = None,
    ) -> "SimulationResult":
        """Create a structured result from a SciPy ODE solution."""
        return cls(
            config=config,
            time=solution.t,
            states=solution.y.T,
            success=bool(solution.success),
            message=str(solution.message),
            initial_state=initial_state,
            control_override=control_override,
            trim_result=trim_result,
            raw_solution=solution,
        )
