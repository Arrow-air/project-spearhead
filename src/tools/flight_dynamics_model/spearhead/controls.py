"""Open-loop control schedules."""

from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np


@dataclass(frozen=True)
class Control:
    de: float
    da: float
    dr: float
    throttle: float


def as_control(control: Control | Mapping[str, float]) -> Control:
    """Return a Control object from a Control instance or control-like dict."""
    if isinstance(control, Control):
        return control
    return Control(
        de=float(control["de"]),
        da=float(control["da"]),
        dr=float(control["dr"]),
        throttle=float(control["throttle"]),
    )


def control_schedule(t: float) -> Control:
    """Return the default open-loop control inputs at time ``t``."""
    step = np.deg2rad(5.0)
    return Control(
        de=step if 4.0 <= t <= 8.0 else 0.0,
        da=step if 10.0 <= t <= 13.0 else 0.0,
        dr=step if 15.0 <= t <= 17.0 else 0.0,
        throttle=0.45,
    )
