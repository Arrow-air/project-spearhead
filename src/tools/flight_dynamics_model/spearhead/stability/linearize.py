"""Finite-difference linearization around a shared-backend trim point."""

from __future__ import annotations

import numpy as np

from ..config import SimulationConfig
from ..controls import Control
from ..dynamics import dynamics
from ..simulation import run_simulation
from .types import (
    LATERAL_CONTROL_INDICES,
    LATERAL_STATE_INDICES,
    LONGITUDINAL_CONTROL_INDICES,
    LONGITUDINAL_STATE_INDICES,
    LinearizationConfig,
    LinearizationResult,
)


def control_to_vector(control: Control) -> np.ndarray:
    """Return the canonical ``[de, da, dr, throttle]`` vector."""
    return np.array([control.de, control.da, control.dr, control.throttle], dtype=float)


def vector_to_control(vector: np.ndarray) -> Control:
    """Return a ``Control`` object from the canonical control vector."""
    values = np.asarray(vector, dtype=float)
    if values.shape != (4,):
        raise ValueError("Expected control vector with shape (4,)")
    return Control(de=values[0], da=values[1], dr=values[2], throttle=values[3])


def _central_difference_state(
    trim_state: np.ndarray,
    trim_control: Control,
    config: SimulationConfig,
    index: int,
    perturbation: float,
) -> np.ndarray:
    plus = trim_state.copy()
    minus = trim_state.copy()
    plus[index] += perturbation
    minus[index] -= perturbation
    f_plus = dynamics(0.0, plus, config.params, control_override=trim_control)
    f_minus = dynamics(0.0, minus, config.params, control_override=trim_control)
    return (f_plus - f_minus) / (2.0 * perturbation)


def _central_difference_control(
    trim_state: np.ndarray,
    trim_control_vector: np.ndarray,
    config: SimulationConfig,
    index: int,
    perturbation: float,
) -> np.ndarray:
    plus = trim_control_vector.copy()
    minus = trim_control_vector.copy()
    plus[index] += perturbation
    minus[index] -= perturbation
    f_plus = dynamics(0.0, trim_state, config.params, control_override=vector_to_control(plus))
    f_minus = dynamics(0.0, trim_state, config.params, control_override=vector_to_control(minus))
    return (f_plus - f_minus) / (2.0 * perturbation)


def _reduced(A: np.ndarray, B: np.ndarray, state_indices: tuple[int, ...], control_indices: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    return A[np.ix_(state_indices, state_indices)], B[np.ix_(state_indices, control_indices)]


def linearize(config: SimulationConfig, settings: LinearizationConfig | None = None) -> LinearizationResult:
    """Linearize the nonlinear FDM around the trim implied by ``config``."""
    settings = LinearizationConfig() if settings is None else settings
    trimmed_result = run_simulation(config)
    trim = trimmed_result.trim_result
    if trim is None:
        raise ValueError("Stability linearization requires config.start_from_trim=True")
    if not bool(trim["success"]):
        raise ValueError(f"Trim failed: {trim['message']}")

    trim_state = np.asarray(trim["x_trim"], dtype=float)
    trim_control = trim["control_trim"]
    trim_control_vector = control_to_vector(trim_control)

    A = np.zeros((12, 12))
    B = np.zeros((12, 4))
    for index in range(12):
        A[:, index] = _central_difference_state(
            trim_state,
            trim_control,
            config,
            index,
            settings.state_perturbation,
        )
    for index in range(4):
        B[:, index] = _central_difference_control(
            trim_state,
            trim_control_vector,
            config,
            index,
            settings.control_perturbation,
        )

    A_long, B_long = _reduced(A, B, LONGITUDINAL_STATE_INDICES, LONGITUDINAL_CONTROL_INDICES)
    A_lat, B_lat = _reduced(A, B, LATERAL_STATE_INDICES, LATERAL_CONTROL_INDICES)
    return LinearizationResult(
        A=A,
        B=B,
        A_longitudinal=A_long,
        B_longitudinal=B_long,
        A_lateral_directional=A_lat,
        B_lateral_directional=B_lat,
        trim_state=trim_state,
        trim_control=trim_control,
    )
