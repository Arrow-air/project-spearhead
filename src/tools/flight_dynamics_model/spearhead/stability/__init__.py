"""Preliminary stability analysis built on the shared FDM backend."""

from .analysis import analyze_stability, extract_static_derivatives
from .export import export_json, export_modes_csv
from .linearize import linearize
from .modes import modal_analysis
from .report import export_markdown, stability_markdown
from .sweep import config_with_cg, run_cg_sweep
from .types import (
    CGSweepPointResult,
    CGSweepResult,
    LinearizationConfig,
    LinearizationResult,
    ModalAnalysisResult,
    ModeResult,
    StabilityResult,
    StaticDerivativeResult,
)

__all__ = [
    "CGSweepPointResult",
    "CGSweepResult",
    "LinearizationConfig",
    "LinearizationResult",
    "ModalAnalysisResult",
    "ModeResult",
    "StabilityResult",
    "StaticDerivativeResult",
    "analyze_stability",
    "config_with_cg",
    "export_json",
    "export_markdown",
    "export_modes_csv",
    "extract_static_derivatives",
    "linearize",
    "modal_analysis",
    "run_cg_sweep",
    "stability_markdown",
]
