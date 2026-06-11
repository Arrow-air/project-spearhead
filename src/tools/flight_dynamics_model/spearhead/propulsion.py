"""Placeholder propulsion model."""

import numpy as np

from .controls import Control
from .params import AircraftParams


def propulsion_model_placeholder(
    control: Control,
    params: AircraftParams,
) -> tuple[np.ndarray, np.ndarray]:
    """Return placeholder propulsion force and moment in body axes."""
    throttle = np.clip(control.throttle, 0.0, 1.0)
    F_prop_b = np.array([params.Tmax * throttle, 0.0, 0.0])
    M_prop_b = np.zeros(3)
    return F_prop_b, M_prop_b

