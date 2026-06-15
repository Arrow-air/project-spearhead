"""CG sweep helpers for stability analysis."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from ..config import SimulationConfig
from .analysis import analyze_stability
from .types import CGSweepPointResult, CGSweepResult, LinearizationConfig


def config_with_cg(config: SimulationConfig, cg_body_xyz: np.ndarray) -> SimulationConfig:
    """Return a config copy with runtime CG/moment-shift location updated."""
    cg = np.asarray(cg_body_xyz, dtype=float)
    if cg.shape != (3,):
        raise ValueError("cg_body_xyz must have shape (3,)")
    params = replace(config.params, cg_body_xyz=cg.copy())
    return replace(config, params=params)


def run_cg_sweep(
    config: SimulationConfig,
    cg_positions: list[np.ndarray] | tuple[np.ndarray, ...],
    linearization_config: LinearizationConfig | None = None,
) -> CGSweepResult:
    """Run stability analysis at multiple runtime CG positions."""
    points = []
    for cg_position in cg_positions:
        cg = np.asarray(cg_position, dtype=float)
        try:
            point_config = config_with_cg(config, cg)
            result = analyze_stability(point_config, linearization_config=linearization_config)
        except Exception as exc:  # noqa: BLE001 - sweep should record point failures.
            points.append(CGSweepPointResult(cg_body_xyz=cg, success=False, error=str(exc)))
        else:
            points.append(CGSweepPointResult(cg_body_xyz=cg, success=True, result=result))
    return CGSweepResult(base_config=config, points=tuple(points))
