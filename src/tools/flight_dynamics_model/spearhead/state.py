"""State containers for the 12-state nonlinear aircraft model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


STATE_COLUMNS = ("pn", "pe", "pd", "u", "v", "w", "phi", "theta", "psi", "p", "q", "r")


@dataclass(frozen=True)
class SimulationState:
    """A timestamped 12-state aircraft state vector."""

    time: float
    vector: np.ndarray

    def __post_init__(self) -> None:
        vector = np.asarray(self.vector, dtype=float)
        if vector.shape != (12,):
            raise ValueError("SimulationState vector must have shape (12,)")
        object.__setattr__(self, "vector", vector)

    def as_dict(self) -> dict[str, float]:
        """Return a JSON/CSV-friendly state dictionary."""
        data = {"time": float(self.time)}
        data.update({name: float(value) for name, value in zip(STATE_COLUMNS, self.vector)})
        return data
