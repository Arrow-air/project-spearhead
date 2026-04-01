"""
aircraft_sizing.py
------------------
Initial aerodynamic sizing for a conventional UAV.
Supports V-tail / inverse V-tail (config 1) and conventional tail (config 2).

Computes:
  - Main wing area from stall speed requirement
  - Tail incidence for cruise trim (lift + moment balance) with fixed tail geometry
  - Induced drag on both surfaces
  - L/D at cruise

Tail geometry (span, chord, area, AR) is fixed by inputs. The trim solver
finds [alpha, tail_incidence] that balance lift and pitching moment at cruise.
This gives zero-elevator trim — the minimum-drag cruise condition.

Sign convention (consistently applied throughout):
  - x positions are relative to CG, positive = aft
    e.g. wing_x_cg = x_ac_wing - x_cg
  - Positive lift = upward
  - Positive pitching moment = nose-up
  - A lift force at aft location (x > 0) creates a nose-DOWN moment (negative)
  - cm_ac follows NACA convention: negative for nose-down (cambered airfoils)

Requires: numpy, scipy
"""

import math
import os
import csv
import numpy as np
from scipy.optimize import fsolve


# ---------------------------------------------------------------------------
# Airfoil database lookup
# ---------------------------------------------------------------------------
_DB_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_RE = 500000


def lookup_airfoil(name: str, re_val: int = _DEFAULT_RE) -> dict:
    """
    Look up an airfoil by name in the CSV database.

    Returns a dict with keys: cl_0, cl_alpha_2d, cl_max, cd_0, cm_ac,
    cd_min, alpha_minD_deg, k_drag, thickness_pct.
    Raises ValueError if not found.
    """
    db_path = os.path.join(_DB_DIR, f"airfoil_db_re{re_val}.csv")
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Airfoil database not found: {db_path}")

    name_lower = name.strip().lower()
    with open(db_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["name"].strip().lower() == name_lower:
                cd_min = float(row["cd_min"]) if row.get("cd_min") else None
                alpha_minD = float(row["alpha_minD_deg"]) if row.get("alpha_minD_deg") else None
                k_drag = float(row["k_drag"]) if row.get("k_drag") else None
                return {
                    "cl_0":           float(row["cl_0"]),
                    "cl_alpha_2d":    float(row["cl_alpha_2d"]),
                    "cl_max":         float(row["cl_max"]),
                    "cd_0":           float(row["cd_0"]),
                    "cm_ac":          float(row["cm_ac"]),
                    "cd_min":         cd_min,
                    "alpha_minD_deg": alpha_minD,
                    "k_drag":         k_drag,
                    "thickness_pct":  float(row["thickness_pct"]) if row.get("thickness_pct") else None,
                }
    raise ValueError(f"Airfoil '{name}' not found in {db_path}")

# ---------------------------------------------------------------------------
# ISA constants
# ---------------------------------------------------------------------------
_G    = 9.80665   # m/s^2
_R    = 287.058   # J/(kg·K)
_P0   = 101325.0  # Pa
_T0   = 288.15    # K
_L    = -0.0065   # K/m  (lapse rate, negative = decreases with altitude)
_T11  = 216.65    # K    (stratosphere isothermal temperature)
_H11  = 11000.0   # m    (tropopause altitude)

# Sutherland viscosity constants
_MU_REF   = 1.716e-5   # Pa·s  (reference dynamic viscosity at T_MU_REF)
_T_MU_REF = 273.15     # K
_C_SUTH   = 110.4      # K     (Sutherland constant for air)


# ---------------------------------------------------------------------------
# Atmosphere
# ---------------------------------------------------------------------------
def isa_atmosphere(altitude_m: float) -> tuple[float, float]:
    """
    ISA standard atmosphere.

    Parameters
    ----------
    altitude_m : float
        Geometric altitude [m]. Valid range: 0 - 20 000 m.

    Returns
    -------
    rho : float  Air density [kg/m3]
    T   : float  Air temperature [K]
    """
    if not (0.0 <= altitude_m <= 20000.0):
        raise ValueError(f"Altitude {altitude_m} m is outside model range [0, 20000] m.")

    if altitude_m <= _H11:
        T = _T0 + _L * altitude_m
        P = _P0 * (T / _T0) ** (-_G / (_L * _R))
    else:
        # Compute conditions at tropopause first
        T11 = _T0 + _L * _H11
        P11 = _P0 * (T11 / _T0) ** (-_G / (_L * _R))
        T   = _T11
        P   = P11 * math.exp(-_G * (altitude_m - _H11) / (_R * _T11))

    rho = P / (_R * T)
    return rho, T


# ---------------------------------------------------------------------------
# Aerodynamic helpers
# ---------------------------------------------------------------------------
def oswald_efficiency(AR: float) -> float:
    """
    Raymer empirical Oswald span efficiency for straight, unswept wings.

        e = 1.78 * (1 - 0.045 * AR^0.68) - 0.64

    Result clamped to [0.5, 1.0].
    """
    if AR <= 0:
        raise ValueError(f"Aspect ratio must be positive, got {AR}.")
    e = 1.78 * (1.0 - 0.045 * AR ** 0.68) - 0.64
    if not (0.5 <= e <= 1.0):
        e = max(0.5, min(1.0, e))
    return e


def finite_cl_alpha(cl_alpha_2d: float, AR: float, e: float) -> float:
    """
    Helmbold equation: 2D → 3D finite-wing lift slope correction.

        CL_alpha = cl_alpha_2d / (sqrt(1 + kappa^2) + kappa)
        where kappa = cl_alpha_2d / (pi * e * AR)

    More accurate than Prandtl lifting-line for moderate/low AR.

    Parameters
    ----------
    cl_alpha_2d : float  2D (infinite-span) lift slope [rad^-1]
    AR          : float  Wing aspect ratio
    e           : float  Oswald efficiency factor

    Returns
    -------
    CL_alpha : float  3D lift slope [rad^-1]
    """
    kappa = cl_alpha_2d / (math.pi * e * AR)
    return cl_alpha_2d / (math.sqrt(1.0 + kappa ** 2) + kappa)


# ---------------------------------------------------------------------------
# Fuselage drag
# ---------------------------------------------------------------------------
def fuselage_drag(
    v: float,
    rho: float,
    T: float,
    diameter: float,
    fus_length: float,
    nose_length: float,
    base_drag_recovery: float = 1.0,
) -> dict:
    """
    Estimate cruise drag of a cylindrical fuselage with a conical nose cone.

    Assumes fully turbulent boundary layer (conservative for a UAV fuselage).
    Base drag is estimated using Hoerner's blunt-base correlation and reduced
    by a recovery factor representing pusher propeller wake ingestion.
    Engine cooling drag is not modelled here.

    Parameters
    ----------
    v          : Cruise speed [m/s]
    rho        : Air density [kg/m3]
    T          : Air temperature [K]
    diameter   : Fuselage outer diameter [m]
    fus_length : Total fuselage length, nose tip to tail [m]
    nose_length: Conical nose cone length [m]
    base_drag_recovery : Fraction of base drag recovered by pusher prop [0-1].
                         1.0 = full recovery (no base drag), 0.0 = no prop.

    Returns
    -------
    dict with drag breakdown and intermediate quantities
    """
    if nose_length >= fus_length:
        raise ValueError("nose_length must be less than fus_length.")

    q = 0.5 * rho * v ** 2

    # -- Dynamic viscosity (Sutherland's law) --
    mu = _MU_REF * (T / _T_MU_REF) ** 1.5 * (_T_MU_REF + _C_SUTH) / (T + _C_SUTH)
    nu = mu / rho                                    # kinematic viscosity [m2/s]

    # -- Reynolds number (based on total fuselage length) --
    Re_fus = v * fus_length / nu

    # -- Turbulent skin friction coefficient (Prandtl-Schlichting) --
    Cf = 0.455 / math.log10(Re_fus) ** 2.58

    # -- Wetted areas --
    r            = diameter / 2.0
    cyl_length   = fus_length - nose_length
    S_wet_nose   = math.pi * r * math.sqrt(r ** 2 + nose_length ** 2)  # cone slant area
    S_wet_cyl    = math.pi * diameter * cyl_length
    S_wet_total  = S_wet_nose + S_wet_cyl

    # -- Form factor (Raymer, body of revolution) --
    f  = fus_length / diameter                       # fineness ratio
    FF = 1.0 + 60.0 / f ** 3 + f / 400.0

    # -- Friction + pressure drag --
    D_friction = q * FF * Cf * S_wet_total

    # -- Base drag (Hoerner blunt-base correlation, 3D corrected) --
    # 2D: Cd_base = 0.029 / sqrt(Cf) referenced to base area (Hoerner 1965).
    # 3D correction: axisymmetric bodies have ~60% of the 2D base drag
    # because the wake can close from all directions (360°), recovering
    # more base pressure than the 2D case (Hoerner, Chapter 3).
    S_base = math.pi * r ** 2
    Cd_base_2d = 0.029 / math.sqrt(Cf)
    Cd_base = 0.60 * Cd_base_2d
    D_base_full = q * Cd_base * S_base
    D_base = D_base_full * (1.0 - base_drag_recovery)

    D_fus = D_friction + D_base

    return {
        "Re_fus":       Re_fus,
        "Cf":           Cf,
        "fineness":     f,
        "FF":           FF,
        "S_wet_nose":   S_wet_nose,
        "S_wet_cyl":    S_wet_cyl,
        "S_wet_total":  S_wet_total,
        "D_friction":   D_friction,
        "D_base_full":  D_base_full,
        "D_base":       D_base,
        "D_fus":        D_fus,
    }


# ---------------------------------------------------------------------------
# Main sizing function
# ---------------------------------------------------------------------------
def size_aircraft(
    # Atmosphere
    altitude_m: float,
    # Flight conditions
    v_cruise: float,
    v_stall: float,
    # Aircraft
    MTOW_kg: float,
    fuselage_width: float,
    fus_length: float,          # m, nose tip to tail end
    nose_length: float,         # m, conical nose cone length
    # Main wing — 2D airfoil data + geometry
    # (these can be omitted when wing_airfoil is set)
    wing_cl_0: float = None,
    wing_cl_alpha_2d: float = None,
    wing_cl_max: float = None,
    wing_cd_0: float = None,
    wing_cm_ac: float = None,
    wing_AR: float = None,
    wing_incidence: float = None,      # rad
    wing_x_cg: float = None,           # m, positive = AC is aft of CG
    # Tail configuration: 1 = V-tail or inverse V-tail, 2 = conventional
    tail_config: int = 1,
    # Tail — 2D airfoil data + geometry
    # (aero coefficients can be omitted when tail_airfoil is set)
    tail_cl_0: float = None,
    tail_cl_alpha_2d: float = None,
    tail_cd_0: float = None,
    tail_cm_ac: float = None,
    tail_proj_span: float = None,      # m, horizontal projected span (V-tail) or full span (conventional)
    tail_chord: float = None,          # m, tail chord (fixed geometry)
    tail_x_cg: float = None,           # m, positive = AC is aft of CG (must be > 0)
    tail_dihedral: float = None,       # rad, V-tail panel dihedral (config 1 only, ignored for config 2)
    # Vertical fin (config 2 only, ignored for config 1)
    vert_chord: float = None,          # m, vertical fin chord
    vert_span: float = None,           # m, vertical fin span (height)
    # Downwash factor: 0.0 = tail sees no downwash (above wing plane),
    #                  1.0 = full downwash (tail in wing wake)
    downwash_factor: float = 1.0,
    # Miscellaneous drag multiplier on total clean drag (gaps, booms, surface roughness)
    # 1.0 = no extra drag, 1.25 = 25% additional drag
    misc_drag_factor: float = 1.0,
    # Base drag recovery by pusher propeller [0-1]
    # 1.0 = full recovery (prop ingests all wake), 0.0 = no prop / no recovery
    base_drag_recovery: float = 1.0,
    # Pusher propeller thrust moment about CG
    # thrust_N     : cruise thrust [N] (0 = ignore thrust moment)
    # z_thrust_m   : vertical offset of thrust line from CG [m]
    #                positive = thrust line BELOW CG (nose-up moment)
    thrust_N: float = 0.0,
    z_thrust_m: float = 0.0,
    # Optional drag polar fit: cd = cd_min + k*(alpha_deg - alpha_minD_deg)^2
    wing_cd_min: float = None,
    wing_alpha_minD_deg: float = None,
    wing_k_drag: float = None,
    tail_cd_min: float = None,
    tail_alpha_minD_deg: float = None,
    tail_k_drag: float = None,
    # Airfoil input mode: 1 = lookup by name, 2 = use manual coefficients above
    wing_airfoil_mode: int = 2,
    tail_airfoil_mode: int = 2,
    wing_airfoil: str = None,
    tail_airfoil: str = None,
    re_db: int = _DEFAULT_RE,
    # Tail incidence override [deg] — when set, skips the trim solver for
    # tail incidence and uses this fixed value. Reports moment imbalance.
    tail_incidence_override_deg: float = None,
) -> dict:
    """
    Size main wing and tail for a conventional UAV.
    Supports V-tail/inverse V-tail (config 1) and conventional tail (config 2).

    Wing area is determined from the stall speed requirement.
    Tail incidence is determined by solving lift + pitch moment balance at cruise
    with fixed tail geometry (zero-elevator trim = minimum drag at cruise).

    Parameters
    ----------
    altitude_m       : Cruise/sizing altitude [m]
    v_cruise         : Cruise speed [m/s]
    v_stall          : Required stall speed [m/s]
    MTOW_kg          : Maximum take-off weight [kg]
    fuselage_width   : Fuselage outer diameter / width [m]; subtracted from wing span for effective area
    fus_length       : Total fuselage length, nose tip to tail [m]
    nose_length      : Conical nose cone length [m]
    wing_cl_0        : 2D lift coefficient at zero angle of attack
    wing_cl_alpha_2d : 2D lift slope [rad^-1] (~ 2pi for thin airfoils)
    wing_cl_max      : 2D maximum lift coefficient (determines stall)
    wing_cd_0        : 2D profile drag coefficient
    wing_cm_ac       : 2D pitching moment coefficient about aerodynamic centre
    wing_AR          : Wing aspect ratio
    wing_incidence   : Wing geometric incidence relative to fuselage datum [rad]
    wing_x_cg        : Wing AC position relative to CG [m], + = aft
    tail_cl_0        : Tail 2D CL at zero AoA
    tail_cl_alpha_2d : Tail 2D lift slope [rad^-1]
    tail_cd_0        : Tail 2D profile drag coefficient
    tail_cm_ac       : Tail 2D pitching moment coefficient about AC
    tail_proj_span   : Horizontal projected span of V-tail [m]
    tail_chord       : Tail chord [m] (fixed geometry; tail area = chord × actual span)
    tail_x_cg        : Tail AC position relative to CG [m], + = aft (must be > 0)
    tail_dihedral    : V-tail panel dihedral angle from horizontal [rad]

    Returns
    -------
    dict  — see result dictionary keys below
    """
    # ------------------------------------------------------------------
    # Airfoil database lookup (mode 1 = DB, mode 2 = manual coefficients)
    # ------------------------------------------------------------------
    if wing_airfoil_mode == 1:
        if wing_airfoil is None:
            raise ValueError("wing_airfoil_mode=1 but wing_airfoil name not set.")
        db = lookup_airfoil(wing_airfoil, re_db)
        wing_cl_0          = db["cl_0"]
        wing_cl_alpha_2d   = db["cl_alpha_2d"]
        wing_cl_max        = db["cl_max"]
        wing_cd_0          = db["cd_0"]
        wing_cm_ac         = db["cm_ac"]
        if db["cd_min"] is not None:
            wing_cd_min         = db["cd_min"]
            wing_alpha_minD_deg = db["alpha_minD_deg"]
            wing_k_drag         = db["k_drag"]

    if tail_airfoil_mode == 1:
        if tail_airfoil is None:
            raise ValueError("tail_airfoil_mode=1 but tail_airfoil name not set.")
        db = lookup_airfoil(tail_airfoil, re_db)
        tail_cl_0          = db["cl_0"]
        tail_cl_alpha_2d   = db["cl_alpha_2d"]
        tail_cd_0          = db["cd_0"]
        tail_cm_ac         = db["cm_ac"]
        if db["cd_min"] is not None:
            tail_cd_min         = db["cd_min"]
            tail_alpha_minD_deg = db["alpha_minD_deg"]
            tail_k_drag         = db["k_drag"]

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------
    required = [
        ("wing_cl_0", wing_cl_0), ("wing_cl_alpha_2d", wing_cl_alpha_2d),
        ("wing_cl_max", wing_cl_max), ("wing_cd_0", wing_cd_0),
        ("wing_cm_ac", wing_cm_ac), ("wing_AR", wing_AR),
        ("wing_incidence", wing_incidence), ("wing_x_cg", wing_x_cg),
        ("tail_cl_0", tail_cl_0), ("tail_cl_alpha_2d", tail_cl_alpha_2d),
        ("tail_cd_0", tail_cd_0), ("tail_cm_ac", tail_cm_ac),
        ("tail_proj_span", tail_proj_span), ("tail_chord", tail_chord),
        ("tail_x_cg", tail_x_cg),
    ]
    if tail_config == 1:
        required.append(("tail_dihedral", tail_dihedral))
    for param_name, param_val in required:
        if param_val is None:
            raise ValueError(f"'{param_name}' is required (provide it directly or via wing_airfoil/tail_airfoil).")

    if v_stall >= v_cruise:
        raise ValueError(f"Stall speed ({v_stall} m/s) must be less than cruise speed ({v_cruise} m/s).")
    if tail_x_cg <= 0:
        raise ValueError(f"tail_x_cg must be positive (tail aft of CG), got {tail_x_cg}.")

    # ------------------------------------------------------------------
    # Atmosphere and dynamic pressures
    # ------------------------------------------------------------------
    rho, T_atm = isa_atmosphere(altitude_m)
    q_cruise = 0.5 * rho * v_cruise ** 2
    q_stall  = 0.5 * rho * v_stall  ** 2
    W        = MTOW_kg * _G             # weight [N]

    # ------------------------------------------------------------------
    # Oswald efficiency and 3D lift slopes
    # ------------------------------------------------------------------
    e_wing = oswald_efficiency(wing_AR)
    # Tail geometry — depends on configuration
    if tail_config == 1:
        # V-tail or inverse V-tail: projected span given, actual span derived
        tail_span_actual = tail_proj_span / math.cos(tail_dihedral)
    elif tail_config == 2:
        # Conventional: span is the full horizontal tail span, no dihedral
        tail_span_actual = tail_proj_span
        tail_dihedral = 0.0
    else:
        raise ValueError(f"tail_config must be 1 (V-tail) or 2 (conventional), got {tail_config}")
    S_tail  = tail_chord * tail_span_actual
    tail_AR = tail_span_actual / tail_chord
    e_tail  = oswald_efficiency(tail_AR)

    CL_alpha_wing = finite_cl_alpha(wing_cl_alpha_2d, wing_AR, e_wing)
    CL_alpha_tail = finite_cl_alpha(tail_cl_alpha_2d, tail_AR, e_tail)

    # 3D correction for cl_0 — same ratio as lift slope, because the
    # zero-lift angle is a geometric property that doesn't change in 3D
    wing_cl_0_3d = wing_cl_0 * (CL_alpha_wing / wing_cl_alpha_2d)
    tail_cl_0_3d = tail_cl_0 * (CL_alpha_tail / tail_cl_alpha_2d)

    # ------------------------------------------------------------------
    # Main wing area from stall speed
    # ------------------------------------------------------------------
    # Scale 2D cl_max to 3D using the ratio of lift slopes
    CL_max_3d = wing_cl_max * (CL_alpha_wing / wing_cl_alpha_2d)

    # Stall: W = q_stall * S_eff * CL_max_3d
    S_wing_eff = W / (q_stall * CL_max_3d)

    # Rectangular untapered wing geometry
    chord_wing  = math.sqrt(S_wing_eff / wing_AR)   # constant chord
    span_eff    = S_wing_eff / chord_wing             # effective (aerodynamic) span
    span_wing   = span_eff + fuselage_width           # physical span including fuselage
    S_wing_gross = span_wing * chord_wing             # gross planform area

    if S_wing_eff <= 0:
        raise ValueError("Computed effective wing area is non-positive. Check fuselage width vs wing span.")

    # ------------------------------------------------------------------
    # V-tail lift factor
    # ------------------------------------------------------------------
    # For an inverse V-tail, each panel's lift vector is tilted inward by the
    # dihedral angle. In symmetric flight the horizontal components cancel;
    # only the vertical (weight-supporting) component counts:
    #   L_tail_vertical = q * S_tail * CL_tail * cos(dihedral)
    # Full horizontal/vertical area split (yaw/roll contribution) to be added later.
    vtail_lift_factor = math.cos(tail_dihedral)

    # ------------------------------------------------------------------
    # Trim solver — two unknowns: [alpha_cruise, tail_incidence]
    # ------------------------------------------------------------------
    # Tail geometry is FIXED. Solver finds the fuselage AoA and tail
    # setting angle that balance lift and moment at cruise (zero elevator).
    #
    # Equations:
    #   (1) Lift balance:   L_wing + L_tail_vertical = W
    #   (2) Moment balance about CG = 0:
    #       q*S_eff*c_wing*cm_ac_wing  - L_wing*wing_x_cg
    #     + q*S_tail*c_tail*cm_ac_tail - L_tail_vert*tail_x_cg = 0
    #
    # Downwash: epsilon = 2*CL_wing / (pi * AR_wing)
    # Effective tail AoA: alpha_tail = alpha - epsilon + i_tail

    # Thrust pitching moment (positive = nose-up)
    M_thrust = thrust_N * z_thrust_m

    alpha0  = math.radians(3.0)
    i_tail0 = math.radians(-2.0)

    if tail_incidence_override_deg is not None:
        # Fixed tail incidence — solve only for alpha (lift balance)
        i_tail_fixed = math.radians(tail_incidence_override_deg)

        def lift_balance(alpha_guess):
            alpha = alpha_guess[0]
            alpha_wing_geo = alpha + wing_incidence
            CL_wing = wing_cl_0_3d + CL_alpha_wing * alpha_wing_geo
            L_wing = q_cruise * S_wing_eff * CL_wing
            epsilon = downwash_factor * 2.0 * CL_wing / (math.pi * wing_AR)
            alpha_tail_eff = alpha - epsilon + i_tail_fixed
            CL_tail = tail_cl_0_3d + CL_alpha_tail * alpha_tail_eff
            L_tail_vert = q_cruise * S_tail * CL_tail * vtail_lift_factor
            return [L_wing + L_tail_vert - W]

        sol, _info, ier, msg = fsolve(lift_balance, [alpha0], full_output=True)
        if ier != 1:
            raise RuntimeError(f"Lift balance solver did not converge: {msg}")
        alpha_cruise = sol[0]
        tail_incidence = i_tail_fixed

    else:
        # Normal trim — solve for both alpha and tail incidence
        def trim_equations(unknowns):
            """Lift + moment balance with fixed tail geometry."""
            alpha, i_tail = unknowns

            # Wing
            alpha_wing_geo = alpha + wing_incidence
            CL_wing = wing_cl_0_3d + CL_alpha_wing * alpha_wing_geo
            L_wing = q_cruise * S_wing_eff * CL_wing

            # Downwash
            epsilon = downwash_factor * 2.0 * CL_wing / (math.pi * wing_AR)

            # Tail
            alpha_tail_eff = alpha - epsilon + i_tail
            CL_tail = tail_cl_0_3d + CL_alpha_tail * alpha_tail_eff
            L_tail_vert = q_cruise * S_tail * CL_tail * vtail_lift_factor

            # Moments
            M_wing_ac = q_cruise * S_wing_eff * chord_wing * wing_cm_ac
            M_tail_ac = q_cruise * S_tail * tail_chord * tail_cm_ac

            F1 = L_wing + L_tail_vert - W
            F2 = (M_wing_ac - L_wing * wing_x_cg
                + M_tail_ac - L_tail_vert * tail_x_cg
                + M_thrust)
            return [F1, F2]

        solution, _info, ier, msg = fsolve(trim_equations, [alpha0, i_tail0], full_output=True)

        if ier != 1:
            raise RuntimeError(f"Trim solver did not converge: {msg}")

        alpha_cruise, tail_incidence = solution

    # ------------------------------------------------------------------
    # Post-solution quantities
    # ------------------------------------------------------------------
    alpha_wing_geo = alpha_cruise + wing_incidence
    CL_wing_cruise = wing_cl_0_3d + CL_alpha_wing * alpha_wing_geo

    epsilon_cruise    = downwash_factor * 2.0 * CL_wing_cruise / (math.pi * wing_AR)
    alpha_tail_eff    = alpha_cruise - epsilon_cruise + tail_incidence
    CL_tail_cruise    = tail_cl_0_3d + CL_alpha_tail * alpha_tail_eff

    # Cruise pitching moment (for override mode; zero when trim-solved)
    L_wing_cruise      = q_cruise * S_wing_eff * CL_wing_cruise
    L_tail_cruise_vert = q_cruise * S_tail * CL_tail_cruise * vtail_lift_factor
    M_cruise = (q_cruise * S_wing_eff * chord_wing * wing_cm_ac
               - L_wing_cruise * wing_x_cg
               + q_cruise * S_tail * tail_chord * tail_cm_ac
               - L_tail_cruise_vert * tail_x_cg
               + M_thrust)

    # Sanity checks
    if abs(math.degrees(alpha_cruise)) > 20.0:
        print(f"  [Warning] Cruise AoA = {math.degrees(alpha_cruise):.1f} deg seems unrealistic.")
    if abs(math.degrees(tail_incidence)) > 15.0:
        print(f"  [Warning] Tail incidence = {math.degrees(tail_incidence):.1f} deg — "
              "consider resizing tail geometry.")

    # ------------------------------------------------------------------
    # Stall condition
    # ------------------------------------------------------------------
    # At stall, wing CL = CL_max_3d (fixed). Tail incidence and area are
    # unchanged from the cruise trim solution. This is NOT a trimmed state —
    # it shows the aerodynamic forces and pitching moment at the onset of stall,
    # which determines whether the aircraft naturally pitches nose-down (safe)
    # or nose-up (dangerous) at stall.

    # Fuselage AoA when wing reaches CL_max_3d
    alpha_stall        = (CL_max_3d - wing_cl_0_3d) / CL_alpha_wing - wing_incidence
    alpha_wing_geo_stall = alpha_stall + wing_incidence

    epsilon_stall      = downwash_factor * 2.0 * CL_max_3d / (math.pi * wing_AR)
    alpha_tail_eff_stall = alpha_stall - epsilon_stall + tail_incidence
    CL_tail_stall      = tail_cl_0_3d + CL_alpha_tail * alpha_tail_eff_stall

    L_wing_stall       = q_stall * S_wing_eff * CL_max_3d        # = W by construction
    L_tail_stall_vert  = q_stall * S_tail * CL_tail_stall * vtail_lift_factor

    # Pitching moment about CG at stall (positive = nose-up)
    M_stall = (q_stall * S_wing_eff * chord_wing  * wing_cm_ac
             - L_wing_stall               * wing_x_cg
             + q_stall * S_tail * tail_chord * tail_cm_ac
             - L_tail_stall_vert          * tail_x_cg
             + M_thrust)

    # Tail incidence required for neutral stability at stall (M_stall = 0)
    # Solve: M_wing_ac - L_wing*x_w + M_tail_ac - q*S_t*(cl0+CLa*(a-e+i))*vf*x_t + M_thrust = 0
    M_rest = (q_stall * S_wing_eff * chord_wing * wing_cm_ac
             - L_wing_stall * wing_x_cg
             + q_stall * S_tail * tail_chord * tail_cm_ac
             + M_thrust)
    tail_lift_per_rad = q_stall * S_tail * CL_alpha_tail * vtail_lift_factor
    CL_tail_base = tail_cl_0_3d + CL_alpha_tail * (alpha_stall - epsilon_stall)
    L_tail_base_vert = q_stall * S_tail * CL_tail_base * vtail_lift_factor
    # M_rest - L_tail_base_vert*x_t - tail_lift_per_rad*i*x_t = 0
    i_tail_neutral_stall = (M_rest - L_tail_base_vert * tail_x_cg) / (tail_lift_per_rad * tail_x_cg)

    if M_stall > 0:
        print(f"  [Warning] Pitching moment at stall is nose-UP ({M_stall:.2f} N.m). "
              "Aircraft may deepen the stall rather than recover.")

    # ------------------------------------------------------------------
    # Drag
    # ------------------------------------------------------------------
    CDi_wing = CL_wing_cruise ** 2 / (math.pi * wing_AR * e_wing)
    CDi_tail = CL_tail_cruise ** 2 / (math.pi * tail_AR * e_tail)

    # Use drag polar fit if available, otherwise fall back to constant cd_0
    # Polar fit is cd = cd_min + k_drag * (alpha_deg - alpha_minD_deg)^2
    if wing_cd_min is not None and wing_alpha_minD_deg is not None and wing_k_drag is not None:
        alpha_w_deg = math.degrees(alpha_wing_geo)
        cd_profile_wing = wing_cd_min + wing_k_drag * (alpha_w_deg - wing_alpha_minD_deg) ** 2
    else:
        cd_profile_wing = wing_cd_0
    D_wing_profile = q_cruise * S_wing_eff * cd_profile_wing
    D_wing_induced = q_cruise * S_wing_eff * CDi_wing
    if tail_cd_min is not None and tail_alpha_minD_deg is not None and tail_k_drag is not None:
        alpha_t_deg = math.degrees(alpha_tail_eff)
        cd_profile_tail = tail_cd_min + tail_k_drag * (alpha_t_deg - tail_alpha_minD_deg) ** 2
    else:
        cd_profile_tail = tail_cd_0
    D_tail_profile = q_cruise * S_tail * cd_profile_tail
    D_tail_induced = q_cruise * S_tail * CDi_tail

    fus = fuselage_drag(v_cruise, rho, T_atm, fuselage_width, fus_length, nose_length, base_drag_recovery)
    D_clean = D_wing_profile + D_wing_induced + D_tail_profile + D_tail_induced + fus["D_fus"]
    D_misc = D_clean * (misc_drag_factor - 1.0)  # additional drag from gaps, booms, surface roughness
    D_total = D_clean + D_misc

    # ------------------------------------------------------------------
    # Tail volume coefficients
    # ------------------------------------------------------------------
    # Effective V_h uses the lift-producing component of tail area:
    # V-tail: S_tail * cos(dihedral) (horizontal projection)
    # Conventional: S_tail (full area, cos(0) = 1)
    S_tail_h_eff = S_tail * vtail_lift_factor
    V_h = (S_tail_h_eff * tail_x_cg) / (S_wing_eff * chord_wing)
    if tail_config == 1:
        # V-tail: vertical area is the vertical projection of the V panels
        S_tail_vert = S_tail * math.sin(tail_dihedral)
    elif tail_config == 2:
        # Conventional: separate vertical fin
        S_tail_vert = (vert_chord * vert_span) if (vert_chord and vert_span) else 0.0
    V_v = (S_tail_vert * tail_x_cg) / (S_wing_eff * span_wing)

    if V_h < 0.40:
        print(f"  [Warning] V_h = {V_h:.3f} is low (typical 0.50-0.70). Consider increasing tail chord or moment arm.")
    if V_v > 0 and V_v < 0.025:
        print(f"  [Warning] V_v = {V_v:.4f} is low (typical 0.03-0.07). Consider increasing tail area or dihedral.")

    # ------------------------------------------------------------------
    # Alpha sweep (L and D vs AoA, 0-14 deg)
    # ------------------------------------------------------------------
    sweep_alpha_deg = list(range(0, 15))
    sweep_L = []
    sweep_D = []
    for _a_deg in sweep_alpha_deg:
        _a_rad = math.radians(_a_deg)
        _a_wing = _a_rad + wing_incidence
        _CL_w = wing_cl_0_3d + CL_alpha_wing * _a_wing
        _L_w = q_cruise * S_wing_eff * _CL_w
        _eps = downwash_factor * 2.0 * _CL_w / (math.pi * wing_AR)
        _a_tail = _a_rad - _eps + tail_incidence
        _CL_t = tail_cl_0_3d + CL_alpha_tail * _a_tail
        _L_t = q_cruise * S_tail * _CL_t * vtail_lift_factor
        sweep_L.append(_L_w + _L_t)
        # Drag
        _a_w_deg = math.degrees(_a_wing)
        if wing_cd_min is not None and wing_k_drag is not None:
            _cd_pw = wing_cd_min + wing_k_drag * (_a_w_deg - wing_alpha_minD_deg) ** 2
        else:
            _cd_pw = wing_cd_0
        _D_pw = q_cruise * S_wing_eff * _cd_pw
        _CDi_w = _CL_w ** 2 / (math.pi * wing_AR * e_wing)
        _D_iw = q_cruise * S_wing_eff * _CDi_w
        _a_t_deg = math.degrees(_a_tail)
        if tail_cd_min is not None and tail_k_drag is not None:
            _cd_pt = tail_cd_min + tail_k_drag * (_a_t_deg - tail_alpha_minD_deg) ** 2
        else:
            _cd_pt = tail_cd_0
        _D_pt = q_cruise * S_tail * _cd_pt
        _CDi_t = _CL_t ** 2 / (math.pi * tail_AR * e_tail)
        _D_it = q_cruise * S_tail * _CDi_t
        _D_clean = _D_pw + _D_iw + _D_pt + _D_it + fus["D_fus"]
        _D_total = _D_clean * misc_drag_factor
        sweep_D.append(_D_total)

    # ------------------------------------------------------------------
    # Result dictionary
    # ------------------------------------------------------------------
    return {
        # Atmosphere
        "altitude_m":           altitude_m,
        "rho":                  rho,
        "T_atm":                T_atm,
        "q_cruise":             q_cruise,
        # Main wing geometry
        "S_wing_eff":           S_wing_eff,
        "S_wing_gross":         S_wing_gross,
        "span_wing":            span_wing,
        "span_wing_eff":        span_eff,
        "chord_wing":           chord_wing,
        "AR_wing_eff":          span_eff ** 2 / S_wing_eff,
        "CL_max_3d":            CL_max_3d,
        "CL_alpha_wing":        CL_alpha_wing,
        "e_wing":               e_wing,
        # Tail geometry
        "S_tail":               S_tail,
        "span_tail":            tail_span_actual,
        "tail_proj_span":       tail_proj_span,
        "tail_chord":           tail_chord,
        "tail_AR":              tail_AR,
        "CL_alpha_tail":        CL_alpha_tail,
        "e_tail":               e_tail,
        "tail_config":          tail_config,
        "V_h":                  V_h,
        "V_v":                  V_v,
        "vert_chord":           vert_chord,
        "vert_span":            vert_span,
        # Trim state (solved)
        "tail_incidence_rad":   tail_incidence,
        "tail_incidence_deg":   math.degrees(tail_incidence),
        "alpha_cruise_rad":     alpha_cruise,
        "alpha_cruise_deg":     math.degrees(alpha_cruise),
        "alpha_wing_geo_deg":   math.degrees(alpha_wing_geo),
        "epsilon_cruise_deg":   math.degrees(epsilon_cruise),
        "alpha_tail_eff_deg":   math.degrees(alpha_tail_eff),
        "CL_wing_cruise":       CL_wing_cruise,
        "CL_tail_cruise":       CL_tail_cruise,
        # Stall condition
        "alpha_stall_deg":          math.degrees(alpha_stall),
        "alpha_wing_geo_stall_deg": math.degrees(alpha_wing_geo_stall),
        "epsilon_stall_deg":        math.degrees(epsilon_stall),
        "alpha_tail_eff_stall_deg": math.degrees(alpha_tail_eff_stall),
        "CL_wing_stall":            CL_max_3d,
        "CL_tail_stall":            CL_tail_stall,
        "L_tail_stall_vert":        L_tail_stall_vert,
        "M_stall":                  M_stall,
        "i_tail_neutral_stall_deg": math.degrees(i_tail_neutral_stall),
        # Drag and performance
        "cd_profile_wing":      cd_profile_wing,
        "cd_profile_tail":      cd_profile_tail,
        "D_wing_profile":       D_wing_profile,
        "D_wing_induced":       D_wing_induced,
        "D_tail_profile":       D_tail_profile,
        "D_tail_induced":       D_tail_induced,
        "D_fus":                fus["D_fus"],
        "D_fus_friction":       fus["D_friction"],
        "D_base_full":          fus["D_base_full"],
        "D_base":               fus["D_base"],
        "base_drag_recovery":   base_drag_recovery,
        "D_misc":               D_misc,
        "misc_drag_factor":     misc_drag_factor,
        "downwash_factor":      downwash_factor,
        "M_thrust":             M_thrust,
        "M_cruise":             M_cruise,
        "tail_incidence_override": tail_incidence_override_deg is not None,
        "D_total":              D_total,
        "L_total":              W,
        "LD":                   W / D_total,
        # Fuselage details
        "fus_Re":               fus["Re_fus"],
        "fus_Cf":               fus["Cf"],
        "fus_FF":               fus["FF"],
        "fus_fineness":         fus["fineness"],
        "fus_S_wet_nose":       fus["S_wet_nose"],
        "fus_S_wet_cyl":        fus["S_wet_cyl"],
        "fus_S_wet_total":      fus["S_wet_total"],
        # Alpha sweep
        "sweep_alpha_deg":      sweep_alpha_deg,
        "sweep_L":              sweep_L,
        "sweep_D":              sweep_D,
        # Solver diagnostics
        "solver_converged":     (ier == 1),
    }


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------
def print_results(r: dict) -> None:
    sep = "=" * 52
    print(f"\n{sep}")
    print("    Aircraft Initial Sizing Results")
    print(sep)

    print(f"\n  ISA Atmosphere")
    print(f"    Altitude          {r['altitude_m']:>10.1f}  m")
    print(f"    Density           {r['rho']:>10.4f}  kg/m3")
    print(f"    Temperature       {r['T_atm']:>10.2f}  K")
    print(f"    Dynamic pressure  {r['q_cruise']:>10.2f}  Pa")

    print(f"\n  Main Wing")
    print(f"    Effective area    {r['S_wing_eff']:>10.4f}  m2")
    print(f"    Gross area        {r['S_wing_gross']:>10.4f}  m2")
    print(f"    Physical span     {r['span_wing']:>10.3f}  m")
    print(f"    Effective span    {r['span_wing_eff']:>10.3f}  m")
    print(f"    Chord             {r['chord_wing']:>10.3f}  m")
    print(f"    Effective AR      {r['AR_wing_eff']:>10.2f}")
    print(f"    3D CL_alpha       {r['CL_alpha_wing']:>10.4f}  rad^-1")
    print(f"    3D CL_max         {r['CL_max_3d']:>10.4f}")
    print(f"    Oswald e          {r['e_wing']:>10.4f}")
    if "wing_thickness_pct" in r and r["wing_thickness_pct"] > 0:
        t_abs_m = r["wing_thickness_pct"] / 100.0 * r["chord_wing"]
        print(f"    Thickness         {r['wing_thickness_pct']:>10.1f}  %c")
        print(f"    Abs. thickness    {t_abs_m*100:>10.1f}  cm")

    tail_type = "V-Tail" if r.get("tail_config", 1) == 1 else "Conventional Tail"
    print(f"\n  {tail_type}  (fixed geometry)")
    print(f"    Panel area total  {r['S_tail']:>10.4f}  m2")
    print(f"    Projected span    {r['tail_proj_span']:>10.3f}  m")
    print(f"    Geometric span    {r['span_tail']:>10.3f}  m")
    print(f"    Chord             {r['tail_chord']:>10.3f}  m")
    print(f"    AR                {r['tail_AR']:>10.2f}")
    print(f"    3D CL_alpha       {r['CL_alpha_tail']:>10.4f}  rad^-1")
    print(f"    Oswald e          {r['e_tail']:>10.4f}")
    vh_note = "LOW" if r['V_h'] < 0.40 else "ok" if r['V_h'] < 0.70 else "high"
    print(f"    V_h               {r['V_h']:>10.3f}        ({vh_note})")
    if r.get("tail_config", 1) == 1:
        vv_note = "LOW" if r['V_v'] < 0.025 else "ok" if r['V_v'] < 0.07 else "high"
        print(f"    V_v               {r['V_v']:>10.4f}       ({vv_note})")
    elif r['V_v'] > 0:
        vv_note = "LOW" if r['V_v'] < 0.025 else "ok" if r['V_v'] < 0.07 else "high"
        print(f"    V_v               {r['V_v']:>10.4f}       ({vv_note}, vertical fin {r['vert_chord']*100:.0f}x{r['vert_span']*100:.0f} cm)")
    else:
        print(f"    V_v                      N/A        (set vert_chord and vert_span for config 2)")
    print(f"                       Sailplane  Homebuilt  GA-single  GA-twin")
    print(f"            V_h ref:     0.50       0.48       0.67      0.81")
    if r.get("tail_config", 1) == 1:
        print(f"            V_v ref:     0.019      0.038      0.044     0.066")

    if r.get("tail_incidence_override"):
        print(f"\n  Cruise State (tail incidence FIXED)")
        print(f"    Tail incidence    {r['tail_incidence_deg']:>10.2f}  deg  (OVERRIDE)")
    else:
        print(f"\n  Trim State (cruise, zero-elevator)")
        print(f"    Tail incidence    {r['tail_incidence_deg']:>10.2f}  deg  (SOLVED)")
    print(f"    Fuselage AoA      {r['alpha_cruise_deg']:>10.2f}  deg")
    print(f"    Wing geo. AoA     {r['alpha_wing_geo_deg']:>10.2f}  deg")
    print(f"    Downwash eps      {r['epsilon_cruise_deg']:>10.2f}  deg")
    print(f"    Tail eff. AoA     {r['alpha_tail_eff_deg']:>10.2f}  deg")
    print(f"    Wing CL           {r['CL_wing_cruise']:>10.4f}")
    print(f"    Tail CL           {r['CL_tail_cruise']:>10.4f}")

    if r.get("M_thrust", 0.0) != 0.0:
        print(f"    Thrust moment     {r['M_thrust']:>10.2f}  N.m  "
              f"({'nose-up' if r['M_thrust'] > 0 else 'nose-down'})")

    if r.get("tail_incidence_override"):
        m = r['M_cruise']
        note = "TRIMMED" if abs(m) < 0.01 else ("nose-up" if m > 0 else "nose-down")
        print(f"    Moment imbalance  {m:>10.2f}  N.m  ({note})")

    stall_moment_note = "NOSE-DOWN (safe)" if r['M_stall'] <= 0 else "NOSE-UP (CHECK!)"
    print(f"\n  Stall Condition  (wing at CL_max, tail fixed from cruise trim)")
    print(f"    Fuselage AoA      {r['alpha_stall_deg']:>10.2f}  deg")
    print(f"    Wing geo. AoA     {r['alpha_wing_geo_stall_deg']:>10.2f}  deg")
    print(f"    Downwash eps      {r['epsilon_stall_deg']:>10.2f}  deg")
    print(f"    Tail eff. AoA     {r['alpha_tail_eff_stall_deg']:>10.2f}  deg")
    print(f"    Wing CL           {r['CL_wing_stall']:>10.4f}  (= CL_max_3d)")
    print(f"    Tail CL           {r['CL_tail_stall']:>10.4f}")
    print(f"    Tail lift (vert)  {r['L_tail_stall_vert']:>10.2f}  N")
    print(f"    Pitch moment      {r['M_stall']:>10.2f}  N.m  {stall_moment_note}")
    print(f"    Neutral i_tail    {r['i_tail_neutral_stall_deg']:>10.2f}  deg  (M=0 at stall)")

    recovery_pct = r['base_drag_recovery'] * 100
    print(f"\n  Fuselage  (conical nose + cylinder, base recovery {recovery_pct:.0f}%)")
    print(f"    Fineness ratio    {r['fus_fineness']:>10.2f}")
    print(f"    Form factor FF    {r['fus_FF']:>10.4f}")
    print(f"    Re (fuselage)     {r['fus_Re']:>10.3e}")
    print(f"    Cf (turbulent)    {r['fus_Cf']:>10.5f}")
    print(f"    Swet nose         {r['fus_S_wet_nose']:>10.4f}  m2")
    print(f"    Swet cylinder     {r['fus_S_wet_cyl']:>10.4f}  m2")
    print(f"    Swet total        {r['fus_S_wet_total']:>10.4f}  m2")
    print(f"    Friction drag     {r['D_fus_friction']:>10.2f}  N")
    print(f"    Base drag (full)  {r['D_base_full']:>10.2f}  N")
    print(f"    Base drag (resid) {r['D_base']:>10.2f}  N")
    print(f"    Fuselage total    {r['D_fus']:>10.2f}  N")

    print(f"\n  Drag (cruise)")
    print(f"    Wing profile      {r['D_wing_profile']:>10.2f}  N")
    print(f"    Wing induced      {r['D_wing_induced']:>10.2f}  N")
    print(f"    Tail profile      {r['D_tail_profile']:>10.2f}  N")
    print(f"    Tail induced      {r['D_tail_induced']:>10.2f}  N")
    print(f"    Fuselage          {r['D_fus']:>10.2f}  N")
    print(f"    Misc (x{r['misc_drag_factor']:.2f})     {r['D_misc']:>10.2f}  N")
    print(f"    Total drag        {r['D_total']:>10.2f}  N")
    print(f"    Total lift        {r['L_total']:>10.2f}  N")
    print(f"    L/D               {r['LD']:>10.2f}")

    print(f"\n  Lift & Drag vs AoA")
    print(f"    {'AoA':>4}  {'Lift':>8}  {'Drag':>8}  {'L/D':>7}")
    for a, L, D in zip(r['sweep_alpha_deg'], r['sweep_L'], r['sweep_D']):
        ld = L / D if D > 0 else 0
        print(f"    {a:>3d}°  {L:>8.1f}  {D:>8.1f}  {ld:>7.1f}")

    print(f"\n{sep}\n")


# ---------------------------------------------------------------------------
# Run directly with: python aircraft_sizing.py
# All inputs come from config.py
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from config import AIRCRAFT, AIRFOIL

    results = size_aircraft(**AIRCRAFT, **AIRFOIL)

    print_results(results)

    # -- Save L & D vs AoA plots --
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    alpha = results["sweep_alpha_deg"]
    L = results["sweep_L"]
    D = results["sweep_D"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(alpha, L, "b-o", markersize=4, linewidth=2)
    ax1.set_xlabel("AoA (deg)", fontsize=12)
    ax1.set_ylabel("Lift (N)", fontsize=12)
    ax1.set_title("Lift vs AoA")
    ax1.xaxis.set_major_locator(MultipleLocator(1))
    ax1.grid(True, alpha=0.3)

    ax2.plot(alpha, D, "r-o", markersize=4, linewidth=2)
    ax2.set_xlabel("AoA (deg)", fontsize=12)
    ax2.set_ylabel("Drag (N)", fontsize=12)
    ax2.set_title("Drag vs AoA")
    ax2.xaxis.set_major_locator(MultipleLocator(1))
    ax2.yaxis.set_major_locator(MultipleLocator(10))
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out_dir = os.path.dirname(os.path.abspath(__file__))
    plot_path = os.path.join(out_dir, "LD_vs_alpha.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved: {plot_path}")
