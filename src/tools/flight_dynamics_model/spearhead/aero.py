"""Aerodynamic models.

ADB v1.1 is preferred when the drag-corrected CSV is available. The original
derivative model remains available as a fallback and as provisional control/rate
increments because ADB v1.1 has no control-surface dimensions yet.
"""

from pathlib import Path

import numpy as np

from .aerodb import dimensional_body_force_moment, load_aerodynamic_database_cached
from .controls import Control
from .params import AircraftParams


def aero_model_placeholder(
    v_b: np.ndarray,
    omega_b: np.ndarray,
    control: Control,
    params: AircraftParams,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | dict[str, float]]]:
    """Compute placeholder aerodynamic forces and moments in body axes."""
    u, v, w = v_b
    p, q, r = omega_b

    V_raw = np.linalg.norm(v_b)
    V = max(V_raw, params.Vmin)
    alpha = np.arctan2(w, u)
    beta = np.arcsin(np.clip(v / V, -1.0, 1.0))
    alpha_eff = np.clip(alpha, -params.alpha_limit, params.alpha_limit)
    beta_eff = np.clip(beta, -params.beta_limit, params.beta_limit)
    qbar = 0.5 * params.rho * V**2

    p_hat = p * params.bref / (2.0 * V)
    q_hat = q * params.cref / (2.0 * V)
    r_hat = r * params.bref / (2.0 * V)

    d = params.aero
    CL = d.CL0 + d.CLalpha * alpha_eff + d.CLde * control.de
    CD = d.CD0 + d.k * CL**2
    CY = d.CYbeta * beta_eff + d.CYdr * control.dr
    Cl = d.Clbeta * beta_eff + d.Clda * control.da + d.Clp * p_hat + d.Clr * r_hat
    Cm = d.Cm0 + d.Cmalpha * alpha_eff + d.Cmde * control.de + d.Cmq * q_hat
    Cn = d.Cnbeta * beta_eff + d.Cndr * control.dr + d.Cnp * p_hat + d.Cnr * r_hat

    CX = -CD * np.cos(alpha_eff) + CL * np.sin(alpha_eff)
    CZ = -CD * np.sin(alpha_eff) - CL * np.cos(alpha_eff)

    F_aero_b = qbar * params.Sref * np.array([CX, CY, CZ])
    M_aero_b = qbar * params.Sref * np.array(
        [params.bref * Cl, params.cref * Cm, params.bref * Cn]
    )

    info = {
        "alpha": alpha,
        "beta": beta,
        "V": V_raw,
        "qbar": qbar,
        "coefficients": {
            "CL": CL,
            "CD": CD,
            "CY": CY,
            "Cl": Cl,
            "Cm": Cm,
            "Cn": Cn,
            "CX": CX,
            "CZ": CZ,
        },
    }
    return F_aero_b, M_aero_b, info


def _placeholder_control_rate_increment(
    alpha: float,
    V: float,
    qbar: float,
    omega_b: np.ndarray,
    control: Control,
    params: AircraftParams,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Return provisional increments not present in ADB v1.1.

    These are not validated control-surface or damping data. They preserve the
    old model's ability to respond to controls/rates without baking those terms
    into the alpha-beta ADB interpolation.
    """
    p, q, r = omega_b
    p_hat = p * params.bref / (2.0 * V)
    q_hat = q * params.cref / (2.0 * V)
    r_hat = r * params.bref / (2.0 * V)

    d = params.aero
    CL = d.CLde * control.de
    CY = d.CYdr * control.dr
    Cl = d.Clda * control.da + d.Clp * p_hat + d.Clr * r_hat
    Cm = d.Cmde * control.de + d.Cmq * q_hat
    Cn = d.Cndr * control.dr + d.Cnp * p_hat + d.Cnr * r_hat

    CX = CL * np.sin(alpha)
    CZ = -CL * np.cos(alpha)

    F_increment_b = qbar * params.Sref * np.array([CX, CY, CZ])
    M_increment_b = qbar * params.Sref * np.array(
        [params.bref * Cl, params.cref * Cm, params.bref * Cn]
    )
    coefficients = {"CX": CX, "CY": CY, "CZ": CZ, "Cl": Cl, "Cm": Cm, "Cn": Cn}
    return F_increment_b, M_increment_b, coefficients


def aero_model_adb(
    v_b: np.ndarray,
    omega_b: np.ndarray,
    control: Control,
    params: AircraftParams,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | str | dict[str, float]]]:
    """Compute aerodynamic loads from ADB v1.1 plus provisional increments."""
    u, v, w = v_b
    V_raw = np.linalg.norm(v_b)
    V = max(V_raw, params.Vmin)
    alpha = np.arctan2(w, u)
    beta = np.arcsin(np.clip(v / V, -1.0, 1.0))
    qbar = 0.5 * params.rho * V**2

    adb_path = Path(params.aero_database_path)
    database = load_aerodynamic_database_cached(adb_path, params.aerodb_reference)
    coefficients = database.coefficients_at(
        np.rad2deg(alpha),
        np.rad2deg(beta),
        clamp=params.aero_database_clamp,
    )
    F_adb_b, M_adb_b = dimensional_body_force_moment(
        coefficients,
        qbar,
        params.aerodb_reference,
        shift_to_cg=params.shift_adb_moments_to_cg,
    )
    F_increment_b, M_increment_b, increment_coefficients = _placeholder_control_rate_increment(
        alpha,
        V,
        qbar,
        omega_b,
        control,
        params,
    )

    info = {
        "source": "adb_v1_1_drag_correction",
        "alpha": alpha,
        "beta": beta,
        "V": V_raw,
        "qbar": qbar,
        "database_path": str(adb_path),
        "coefficients": coefficients,
        "provisional_increment_coefficients": increment_coefficients,
    }
    return F_adb_b + F_increment_b, M_adb_b + M_increment_b, info


def aero_model(
    v_b: np.ndarray,
    omega_b: np.ndarray,
    control: Control,
    params: AircraftParams,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | str | dict[str, float]]]:
    """Compute aerodynamic loads, preferring ADB v1.1 when available."""
    adb_path = Path(params.aero_database_path)
    if params.use_aero_database and adb_path.exists():
        return aero_model_adb(v_b, omega_b, control, params)

    F_aero_b, M_aero_b, info = aero_model_placeholder(v_b, omega_b, control, params)
    info["source"] = "placeholder_derivative_fallback"
    return F_aero_b, M_aero_b, info
