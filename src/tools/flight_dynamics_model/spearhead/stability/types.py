"""Typed result structures for stability analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..config import SimulationConfig
from ..controls import Control


STATE_LABELS = ("pn", "pe", "pd", "u", "v", "w", "phi", "theta", "psi", "p", "q", "r")
CONTROL_LABELS = ("de", "da", "dr", "throttle")
LONGITUDINAL_STATE_LABELS = ("u", "w", "q", "theta")
LONGITUDINAL_CONTROL_LABELS = ("de", "throttle")
LATERAL_STATE_LABELS = ("v", "p", "r", "phi", "psi")
LATERAL_CONTROL_LABELS = ("da", "dr")
LONGITUDINAL_STATE_INDICES = (3, 5, 10, 7)
LONGITUDINAL_CONTROL_INDICES = (0, 3)
LATERAL_STATE_INDICES = (4, 9, 11, 6, 8)
LATERAL_CONTROL_INDICES = (1, 2)

VALIDATION_WARNING = (
    "PRELIMINARY / NOT VALIDATED: inertia, control derivatives, rate damping, "
    "and propulsion are placeholders; ADB v1.1 contains alpha-beta data only."
)


@dataclass(frozen=True)
class LinearizationConfig:
    """Numerical perturbation settings for finite-difference linearization."""

    state_perturbation: float = 1e-5
    control_perturbation: float = 1e-5

    def __post_init__(self) -> None:
        if self.state_perturbation <= 0.0:
            raise ValueError("state_perturbation must be positive")
        if self.control_perturbation <= 0.0:
            raise ValueError("control_perturbation must be positive")


@dataclass(frozen=True)
class LinearizationResult:
    """Full and reduced linear systems around a trimmed operating point."""

    A: np.ndarray
    B: np.ndarray
    A_longitudinal: np.ndarray
    B_longitudinal: np.ndarray
    A_lateral_directional: np.ndarray
    B_lateral_directional: np.ndarray
    trim_state: np.ndarray
    trim_control: Control
    state_labels: tuple[str, ...] = STATE_LABELS
    control_labels: tuple[str, ...] = CONTROL_LABELS
    longitudinal_state_labels: tuple[str, ...] = LONGITUDINAL_STATE_LABELS
    longitudinal_control_labels: tuple[str, ...] = LONGITUDINAL_CONTROL_LABELS
    lateral_state_labels: tuple[str, ...] = LATERAL_STATE_LABELS
    lateral_control_labels: tuple[str, ...] = LATERAL_CONTROL_LABELS


@dataclass(frozen=True)
class ModeResult:
    """Eigenvalue-derived modal metrics."""

    eigenvalue: complex
    mode: str
    natural_frequency_rad_s: float
    damping_ratio: float | None
    oscillation_period_s: float | None
    time_to_half_s: float | None
    time_to_double_s: float | None
    classification: str

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON/CSV-friendly representation."""
        return {
            "mode": self.mode,
            "eigenvalue_real": float(self.eigenvalue.real),
            "eigenvalue_imag": float(self.eigenvalue.imag),
            "natural_frequency_rad_s": self.natural_frequency_rad_s,
            "damping_ratio": self.damping_ratio,
            "oscillation_period_s": self.oscillation_period_s,
            "time_to_half_s": self.time_to_half_s,
            "time_to_double_s": self.time_to_double_s,
            "classification": self.classification,
        }


@dataclass(frozen=True)
class ModalAnalysisResult:
    """Eigenvalues and modes for one subsystem."""

    subsystem: str
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    modes: tuple[ModeResult, ...]


@dataclass(frozen=True)
class StaticDerivativeResult:
    """Static aerodynamic derivatives estimated around trim."""

    derivatives: dict[str, float]
    source: str
    alpha_deg: float
    beta_deg: float


@dataclass(frozen=True)
class StabilityResult:
    """Complete stability analysis for one simulation configuration."""

    config: SimulationConfig
    linearization: LinearizationResult
    longitudinal: ModalAnalysisResult
    lateral_directional: ModalAnalysisResult
    static_derivatives: StaticDerivativeResult
    trim_result: dict[str, Any]
    warning: str = VALIDATION_WARNING


@dataclass(frozen=True)
class CGSweepPointResult:
    """Stability result or failure for one CG point."""

    cg_body_xyz: np.ndarray
    success: bool
    result: StabilityResult | None = None
    error: str | None = None


@dataclass(frozen=True)
class CGSweepResult:
    """CG sweep over runtime moment-shift positions."""

    base_config: SimulationConfig
    points: tuple[CGSweepPointResult, ...]
    warning: str = VALIDATION_WARNING

    @property
    def successful_points(self) -> tuple[CGSweepPointResult, ...]:
        """Return sweep points with successful stability results."""
        return tuple(point for point in self.points if point.success)
