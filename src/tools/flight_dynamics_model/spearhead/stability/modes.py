"""Eigenvalue and modal-metric helpers for stability analysis."""

from __future__ import annotations

import math

import numpy as np

from .types import ModalAnalysisResult, ModeResult


def _modal_metrics(eigenvalue: complex) -> tuple[float, float | None, float | None, float | None, float | None]:
    sigma = float(eigenvalue.real)
    omega = float(eigenvalue.imag)
    natural_frequency = abs(eigenvalue)
    damping_ratio = None if natural_frequency == 0.0 else -sigma / natural_frequency
    period = None if abs(omega) == 0.0 else 2.0 * math.pi / abs(omega)
    time_to_half = math.log(0.5) / sigma if sigma < 0.0 else None
    time_to_double = math.log(2.0) / sigma if sigma > 0.0 else None
    return natural_frequency, damping_ratio, period, time_to_half, time_to_double


def _classification(eigenvalue: complex) -> str:
    if eigenvalue.real < 0.0:
        return "stable"
    if eigenvalue.real > 0.0:
        return "unstable"
    return "neutral"


def classify_longitudinal_modes(eigenvalues: np.ndarray) -> tuple[str, ...]:
    """Classify longitudinal roots with preliminary frequency heuristics."""
    labels = []
    for eigenvalue in eigenvalues:
        if abs(eigenvalue.imag) < 1e-9:
            labels.append("non-oscillatory")
        elif abs(eigenvalue.imag) > 1.0:
            labels.append("short-period")
        else:
            labels.append("phugoid")
    return tuple(labels)


def classify_lateral_directional_modes(eigenvalues: np.ndarray) -> tuple[str, ...]:
    """Classify lateral-directional roots with simple conventional heuristics."""
    labels = ["yaw/directional-related"] * len(eigenvalues)
    real_indices = []
    for index, eigenvalue in enumerate(eigenvalues):
        if abs(eigenvalue.imag) > 1e-9:
            labels[index] = "dutch-roll-like"
        else:
            real_indices.append(index)

    near_zero = [index for index in real_indices if abs(eigenvalues[index].real) < 1e-9]
    if near_zero:
        labels[near_zero[0]] = "heading integrator"
        real_indices = [index for index in real_indices if index != near_zero[0]]

    stable_real = [index for index in real_indices if eigenvalues[index].real < 0.0]
    if stable_real:
        roll_index = min(stable_real, key=lambda index: eigenvalues[index].real)
        labels[roll_index] = "roll subsidence"
        remaining = [index for index in stable_real if index != roll_index]
        if remaining:
            spiral_index = max(remaining, key=lambda index: eigenvalues[index].real)
            labels[spiral_index] = "spiral"

    unstable_real = [index for index in real_indices if eigenvalues[index].real > 0.0]
    if unstable_real:
        spiral_index = min(unstable_real, key=lambda index: abs(eigenvalues[index].real))
        labels[spiral_index] = "spiral"

    return tuple(labels)


def modal_analysis(A: np.ndarray, subsystem: str, labels: tuple[str, ...] | None = None) -> ModalAnalysisResult:
    """Compute eigenvalues, eigenvectors, and modal metrics for ``A``."""
    eigenvalues, eigenvectors = np.linalg.eig(A)
    if labels is None:
        labels = tuple("unclassified" for _ in eigenvalues)
    if len(labels) != eigenvalues.size:
        raise ValueError("Expected one mode label per eigenvalue")

    modes = []
    for eigenvalue, label in zip(eigenvalues, labels):
        natural_frequency, damping_ratio, period, time_to_half, time_to_double = _modal_metrics(
            complex(eigenvalue)
        )
        modes.append(
            ModeResult(
                eigenvalue=complex(eigenvalue),
                mode=label,
                natural_frequency_rad_s=natural_frequency,
                damping_ratio=damping_ratio,
                oscillation_period_s=period,
                time_to_half_s=time_to_half,
                time_to_double_s=time_to_double,
                classification=_classification(complex(eigenvalue)),
            )
        )

    return ModalAnalysisResult(
        subsystem=subsystem,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        modes=tuple(modes),
    )
