"""Centralized body-frame force and moment computation."""

import numpy as np

from .aero import aero_model
from .controls import ControlInput, control_at
from .params import AircraftParams
from .propulsion import propulsion_model_placeholder
from .rotations import euler_body_to_ned

COMPONENTS = ("wing", "tail", "fuselage", "propulsion", "gravity", "total")


def _entry(force_b: np.ndarray, moment_b: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "force_b": np.asarray(force_b, dtype=float),
        "moment_b": np.asarray(moment_b, dtype=float),
    }


def _placeholder_aero_component_split(
    F_aero_b: np.ndarray,
    M_aero_b: np.ndarray,
) -> dict[str, dict[str, np.ndarray]]:
    """Split total placeholder aero loads into bookkeeping components.

    This is not validated component-level physics. It only preserves a useful
    debugging interface until the Spearhead aerodynamic database provides real
    wing, tail, and fuselage loads.
    """
    wing_fraction = 0.65
    tail_fraction = 0.25

    F_wing_b = wing_fraction * F_aero_b
    M_wing_b = wing_fraction * M_aero_b
    F_tail_b = tail_fraction * F_aero_b
    M_tail_b = tail_fraction * M_aero_b

    # Residual assignment keeps wing + tail + fuselage exactly equal to the
    # current total placeholder aerodynamic force and moment.
    F_fuselage_b = F_aero_b - F_wing_b - F_tail_b
    M_fuselage_b = M_aero_b - M_wing_b - M_tail_b

    return {
        "wing": _entry(F_wing_b, M_wing_b),
        "tail": _entry(F_tail_b, M_tail_b),
        "fuselage": _entry(F_fuselage_b, M_fuselage_b),
    }


def compute_force_moment_breakdown(
    t: float,
    x: np.ndarray,
    params: AircraftParams,
    control_override: ControlInput | None = None,
) -> dict[str, dict[str, np.ndarray]]:
    """Return body-frame force and moment contributions about the CG.

    Forces are in N, moments are in N*m, and all vectors are expressed in the
    aircraft body frame.
    """
    v_b = x[3:6]
    phi, theta, psi = x[6:9]
    omega_b = x[9:12]

    control = control_at(t, control_override)
    F_aero_b, M_aero_b, _ = aero_model(v_b, omega_b, control, params)
    F_prop_b, M_prop_b = propulsion_model_placeholder(control, params)

    R_body_to_ned = euler_body_to_ned(phi, theta, psi)
    F_gravity_ned = np.array([0.0, 0.0, params.mass * params.g])
    F_gravity_b = R_body_to_ned.T @ F_gravity_ned

    breakdown = _placeholder_aero_component_split(F_aero_b, M_aero_b)
    breakdown["propulsion"] = _entry(F_prop_b, M_prop_b)
    breakdown["gravity"] = _entry(F_gravity_b, np.zeros(3))

    non_total_components = [name for name in COMPONENTS if name != "total"]
    F_total_b = sum((breakdown[name]["force_b"] for name in non_total_components), np.zeros(3))
    M_total_b = sum((breakdown[name]["moment_b"] for name in non_total_components), np.zeros(3))
    breakdown["total"] = _entry(F_total_b, M_total_b)

    return breakdown


def compute_force_moment_history(
    t: np.ndarray,
    y: np.ndarray,
    params: AircraftParams,
    control_override: ControlInput | None = None,
) -> dict[str, np.ndarray | dict[str, np.ndarray]]:
    """Compute force and moment histories for each saved simulation point."""
    t = np.asarray(t)
    y = np.asarray(y)
    if y.shape[0] != 12:
        raise ValueError("Expected y with shape (12, N)")
    if y.shape[1] != t.size:
        raise ValueError("Expected t length to match y.shape[1]")

    history: dict[str, np.ndarray | dict[str, np.ndarray]] = {"time": t.copy()}
    for component in COMPONENTS:
        history[component] = {
            "force_b": np.zeros((t.size, 3)),
            "moment_b": np.zeros((t.size, 3)),
        }

    for i, ti in enumerate(t):
        breakdown = compute_force_moment_breakdown(
            float(ti),
            y[:, i],
            params,
            control_override=control_override,
        )
        for component in COMPONENTS:
            component_history = history[component]
            component_history["force_b"][i, :] = breakdown[component]["force_b"]
            component_history["moment_b"][i, :] = breakdown[component]["moment_b"]

    return history
