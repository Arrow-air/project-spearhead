"""
Spearhead UAV - RCGF Stinger 35cc RE Engine Model
===================================================
Produces 8 subplots:
  Row 0: Power & Torque vs RPM  |  Static Thrust vs RPM (multiple props)
  Row 1: Fuel Consumption vs RPM  |  Endurance vs Tank Size
  Row 2: Force Balance vs Speed  |  Airspeed vs Time
  Row 3: Acceleration vs Speed   |  Distance vs Time

Physics notes
-------------
Power curve  : anchored to rated HP @ rated RPM, shaped by a volumetric-
               efficiency curve typical of piston-port 2-strokes.
               P = P_anchor * (VE/VE_anchor) * (RPM/RPM_anchor)
               which follows from P ~ BMEP * RPM ~ VE * RPM for fixed
               displacement.

Static thrust: combined engine+propeller model using momentum theory.
               T ~ P^(2/3) * (rho*A)^(1/3), and if P ~ RPM^2 (approx.
               for 2-strokes above torque peak), then T ~ RPM^(4/3).
               Anchored to a manufacturer test point.
               NOTE: the pitch scaling (p/p_ref)^(1/3) is empirical, not
               derived from actuator disk theory.

Thrust in flight: Ct(J) quadratic model,
               Ct = Ct0 * max(0, 1 - (J / J_zero)^2)
               where J = V / (n * D)  and  J_zero = eta_pitch * (p/D)
               eta_pitch ~ 1.15-1.25 for fixed-pitch props (zero-thrust
               advance ratio is beyond geometric pitch speed).

Drag         : D = D_cruise * (rho/rho_sl) * (V / V_cruise)^2
               Calibrated to total drag at cruise from the sizing tool.
               Simple quadratic model -- detailed aero breakdown (parasite
               + induced) belongs in aircraft_sizing.py.

Transition   : Euler integration, dt = 0.02 s, 1-D horizontal level flight.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

G      = 9.81    # m/s^2
RHO_SL = 1.225   # kg/m^3  sea-level ISA


# =============================================================================
# Input configuration
# =============================================================================
# All tunable parameters live here. Modify this dict or pass overrides to
# make_plots().

DEFAULT_CONFIG = dict(
    # --- Engine ---
    rated_hp        = 4.1,       # peak shaft HP
    rated_rpm       = 9000.0,    # RPM at rated HP
    ve_peak         = 0.78,      # peak volumetric efficiency (0-1)
    rpm_ve_peak     = 7500.0,    # RPM at which VE peaks
    bsfc            = 520.0,     # brake specific fuel consumption [g/kWh]

    # --- Manufacturer thrust anchor ---
    anchor_thrust_kg = 9.0,      # static thrust [kg-force]
    anchor_thrust_rpm = 7900.0,  # RPM for that thrust
    anchor_prop_d    = 20.0,     # prop diameter [inches] used in anchor test
    anchor_prop_p    = 10.0,      # prop pitch [inches] used in anchor test
    anchor_alt_m    = 100.0,       # altitude of manufacturer thrust data point [m]

    # --- Props to plot on static-thrust chart ---
    # (diameter_in, pitch_in, color)
    props = [
        (19, 8,  "#888780"),
        (20, 8,  "#BA7517"),
        (20, 10, "#1D9E75"),
        (22, 12, "#7F77DD"),
    ],

    # --- Aircraft ---
    mass_kg         = 25.0,      # MTOW [kg]
    drag_cruise_N   = 30.0,      # total drag at cruise speed [N] (from sizing tool)

    # --- Cruise ---
    v_cruise        = 25.0,      # cruise speed [m/s]
    endurance_target_hr = 4.0,   # target endurance [hr]
    battery_Wh      = 300.0,     # battery capacity [Wh] for charge time estimate
    eta_generator   = 0.70,      # generator efficiency (mechanical -> electrical)

    # --- Transition simulation ---
    sim_rpm         = 8000,      # engine RPM during transition
    sim_prop_d      = 20,        # prop diameter [inches]
    sim_prop_p      = 10,        # prop pitch [inches]
    eta_pitch       = 1.20,      # zero-thrust advance ratio factor
    eta_prop_max    = 0.70,      # peak propeller efficiency (typical 0.60-0.75 for fixed-pitch)
    cruise_alt_m    = 2000.0,       # mission / cruise altitude [m]
    v0              = 0.0,       # initial airspeed [m/s]
    v_target        = 25.0,      # target airspeed [m/s]
    dt              = 0.02,      # integration timestep [s]
    t_max           = 60.0,     # max sim time [s]

    # --- RPM sweep range ---
    rpm_min         = 1800.0,
    rpm_max         = 9000.0,
    rpm_points      = 400,
)


# =============================================================================
# Atmosphere
# =============================================================================
def air_density(alt_m: float) -> float:
    """ISA troposphere air density [kg/m^3] at altitude [m]."""
    return RHO_SL * (1 - 2.2557e-5 * alt_m) ** 5.2559


# =============================================================================
# Engine model
# =============================================================================
def ve_shape(rpm: np.ndarray, ve_peak: float, rpm_peak: float) -> np.ndarray:
    """
    Volumetric efficiency curve for a piston-port 2-stroke.
    Quadratic + cubic around peak gives asymmetric falloff (drops faster
    above peak than below, typical of port-timed engines). Clipped at 0.15
    to prevent negative values at extreme RPMs.
    """
    x = rpm / rpm_peak
    raw = ve_peak * (1 - 0.55 * (x - 1)**2 - 0.07 * (x - 1)**3)
    return np.maximum(0.15, raw)


def power_curve(rpms: np.ndarray, cfg: dict) -> np.ndarray:
    """
    Shaft power [kW] at each RPM.
    P = P_anchor * (VE(RPM) / VE(RPM_anchor)) * (RPM / RPM_anchor)
    """
    anchor_kw   = cfg["rated_hp"] * 0.7457
    anchor_rpm  = cfg["rated_rpm"]
    ve_peak     = cfg["ve_peak"]
    rpm_ve_peak = cfg["rpm_ve_peak"]

    ve      = ve_shape(rpms, ve_peak, rpm_ve_peak)
    ve_anch = ve_shape(np.array([anchor_rpm]), ve_peak, rpm_ve_peak)[0]
    ratio   = (ve / ve_anch) * (rpms / anchor_rpm)
    return np.maximum(0.0, anchor_kw * ratio)


def torque_curve(rpms: np.ndarray, powers_kw: np.ndarray) -> np.ndarray:
    """Torque [Nm] = Power [W] / omega [rad/s]."""
    omega = rpms * (2 * np.pi / 60)
    return (powers_kw * 1000) / omega


# =============================================================================
# Propeller / thrust model
# =============================================================================
def prop_scale(d_in: float, p_in: float, cfg: dict) -> float:
    """
    Scale static thrust relative to the anchor prop.
    Diameter: (D/D_ref)^(4/3) from momentum theory T ~ (P^2 * rho * A)^(1/3)
              with A ~ D^2 -> T ~ D^(2/3) at constant power, but since the
              engine RPM^2 also scales with disk area, net is ~D^(4/3).
    Pitch:    (p/p_ref)^(1/3) is empirical - different pitch changes blade
              loading and Ct0 but NOT via a simple power law. Treat with caution.
    """
    d_ref = cfg["anchor_prop_d"]
    p_ref = cfg["anchor_prop_p"]
    return (d_in / d_ref) ** (4/3) * (p_in / p_ref) ** (1/3)


def static_thrust(rpms: np.ndarray,
                  d_in: float, p_in: float,
                  rho: float, cfg: dict) -> np.ndarray:
    """
    Static thrust [N] via RPM^(4/3) combined engine+prop scaling.
    Assumes P ~ RPM^2 (approx. for 2-stroke above torque peak) and
    momentum theory T ~ P^(2/3), giving T ~ RPM^(4/3).
    """
    anchor_T_N   = cfg["anchor_thrust_kg"] * G
    anchor_T_rpm = cfg["anchor_thrust_rpm"]
    rho_anchor   = air_density(cfg["anchor_alt_m"])
    ps    = prop_scale(d_in, p_in, cfg)
    rho_f = rho / rho_anchor  # density ratio relative to anchor conditions
    return anchor_T_N * ps * rho_f * (rpms / anchor_T_rpm) ** (4/3)


def advance_ratio(v: float, rpm: float, d_in: float) -> float:
    """J = V / (n * D), where n = RPM/60 and D in meters."""
    n = rpm / 60.0
    D = d_in * 0.0254
    return v / (n * D + 1e-9)


def thrust_in_flight(v: float, rpm: float,
                     d_in: float, p_in: float,
                     rho: float, cfg: dict) -> float:
    """
    In-flight thrust [N] at airspeed v [m/s].
    Ct(J) = Ct0 * max(0, 1 - (J / J_zero)^2)
    J_zero = eta_pitch * (p/D) -- actual zero-thrust advance ratio.
    eta_pitch > 1 because real props produce thrust beyond geometric
    pitch speed. Typical 1.10-1.30 for fixed-pitch RC props.
    """
    eta_pitch = cfg["eta_pitch"]
    T0     = static_thrust(np.array([rpm]), d_in, p_in, rho, cfg)[0]
    J      = advance_ratio(v, rpm, d_in)
    J_geo  = p_in / d_in
    J_zero = eta_pitch * J_geo
    ct_ratio = max(0.0, 1.0 - (J / J_zero) ** 2)
    return T0 * ct_ratio


# =============================================================================
# Propeller efficiency & fuel consumption
# =============================================================================
def prop_efficiency(v: float, rpm: float, d_in: float, p_in: float,
                    cfg: dict) -> float:
    """
    Propeller efficiency eta = T*V / P_shaft.
    Model: eta peaks at ~70-80% of J_zero (the advance ratio where thrust
    goes to zero), and falls off on both sides.
    At static (J=0): eta = 0 (no useful work despite producing thrust).
    At J_zero: eta = 0 (no thrust).
    Shape: eta = eta_max * 4 * x * (1 - x)  where x = J/J_zero
    This parabola peaks at x=0.5 with value eta_max.
    """
    eta_max = cfg["eta_prop_max"]
    J = advance_ratio(v, rpm, d_in)
    J_geo = p_in / d_in
    J_zero = cfg["eta_pitch"] * J_geo
    if J_zero <= 0 or J <= 0:
        return 0.0
    x = J / J_zero
    if x >= 1.0:
        return 0.0
    return eta_max * 4.0 * x * (1.0 - x)


def shaft_power_from_thrust(T: float, v: float, rpm: float,
                            d_in: float, p_in: float,
                            cfg: dict) -> float:
    """
    Shaft power [kW] required to produce thrust T [N] at airspeed v [m/s].
    P_shaft = T * V / eta_prop.  At very low speed, falls back to engine
    power curve (WOT) as a ceiling.
    """
    eta = prop_efficiency(v, rpm, d_in, p_in, cfg)
    if eta < 0.01 or v < 1.0:
        # At static / very low speed, use engine power curve as upper bound
        return float(power_curve(np.array([rpm]), cfg)[0])
    return (T * v / eta) / 1000.0  # W -> kW


def fuel_flow(powers_kw: np.ndarray, bsfc: float) -> np.ndarray:
    """Fuel flow [ml/hr] from shaft power [kW] and BSFC [g/kWh]. Gasoline 0.74 kg/L."""
    return (powers_kw * bsfc) / 0.74


# =============================================================================
# Drag model
# =============================================================================
def aero_drag(v: float, cfg: dict) -> float:
    """
    Aerodynamic drag [N] = D_cruise * (V/V_cruise)^2.
    Quadratic scaling from a single calibration point (total drag at cruise
    speed and altitude, supplied by the sizing tool). No additional density
    correction — drag_cruise_N already includes the correct density.
    """
    D_cr  = cfg["drag_cruise_N"]
    V_cr  = cfg["v_cruise"]
    return D_cr * (v / V_cr) ** 2


# =============================================================================
# Cruise RPM solver
# =============================================================================
def find_cruise_rpm(d_in: float, p_in: float, rho: float, cfg: dict) -> float:
    """
    Find the RPM at which thrust equals drag at cruise speed.
    Bisection search over [rpm_min, rpm_max]. Returns NaN if no solution.
    """
    v_cr = cfg["v_cruise"]
    D_cr = aero_drag(v_cr, cfg)
    lo, hi = cfg["rpm_min"], cfg["rpm_max"]

    # Check that a solution exists (thrust at max RPM must exceed drag)
    T_hi = thrust_in_flight(v_cr, hi, d_in, p_in, rho, cfg)
    T_lo = thrust_in_flight(v_cr, lo, d_in, p_in, rho, cfg)
    if T_hi < D_cr:
        return float("nan")  # can't reach cruise even at max RPM
    if T_lo >= D_cr:
        return lo  # even idle RPM is enough

    for _ in range(60):  # ~18 digits of precision
        mid = (lo + hi) / 2
        T_mid = thrust_in_flight(v_cr, mid, d_in, p_in, rho, cfg)
        if T_mid < D_cr:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# =============================================================================
# Transition simulation
# =============================================================================
def simulate_transition(cfg: dict):
    """
    Euler integration of 1-D horizontal acceleration.
    Returns (time, velocity, acceleration, distance) arrays.

    Limitation: drag model is a simple v^2 scaling from cruise drag.
    Does not model ground roll friction or climb-out.
    """
    rho = air_density(cfg["cruise_alt_m"])
    rpm = cfg["sim_rpm"]
    d_in = cfg["sim_prop_d"]
    p_in = cfg["sim_prop_p"]
    v_target = cfg["v_target"]
    mass = cfg["mass_kg"]
    dt = cfg["dt"]
    t_max = cfg["t_max"]

    t, v, x = 0.0, cfg["v0"], 0.0
    ts, vs, accs, xs = [], [], [], []

    while t < t_max:
        T    = thrust_in_flight(v, rpm, d_in, p_in, rho, cfg)
        D    = aero_drag(v, cfg)
        Fnet = T - D
        a    = Fnet / mass

        ts.append(t)
        vs.append(v)
        accs.append(a)
        xs.append(x)

        v = max(0.0, v + a * dt)
        x += v * dt
        t += dt

        if v >= v_target:
            ts.append(t)
            vs.append(v)
            accs.append(0.0)
            xs.append(x)
            break

    return (np.array(ts), np.array(vs),
            np.array(accs), np.array(xs))


# =============================================================================
# Plotting
# =============================================================================
def make_plots(overrides: dict | None = None):
    """
    Generate the full engine model figure.
    Pass overrides dict to change any DEFAULT_CONFIG parameter.
    """
    cfg = {**DEFAULT_CONFIG}
    if overrides:
        cfg.update(overrides)

    rpms = np.linspace(cfg["rpm_min"], cfg["rpm_max"], cfg["rpm_points"])
    rho_anchor = air_density(cfg["anchor_alt_m"])
    rho_cruise = air_density(cfg["cruise_alt_m"])

    # Engine curves
    pwr  = power_curve(rpms, cfg)
    trq  = torque_curve(rpms, pwr)
    fuel = fuel_flow(pwr, cfg["bsfc"])

    # Transition simulation
    ts, vs, accs, xs = simulate_transition(cfg)

    sim_d = cfg["sim_prop_d"]
    sim_p = cfg["sim_prop_p"]
    sim_rpm = cfg["sim_rpm"]
    v_target = cfg["v_target"]

    # Force balance vs speed (at cruise altitude)
    v_range = np.linspace(0, 40, 300)
    T_vs_v = np.array([thrust_in_flight(v, sim_rpm, sim_d, sim_p, rho_cruise, cfg)
                        for v in v_range])
    D_vs_v = np.array([aero_drag(v, cfg) for v in v_range])
    F_net  = T_vs_v - D_vs_v

    # Summary metrics
    t_total  = ts[-1]
    x_total  = xs[-1]
    T0_val   = static_thrust(np.array([sim_rpm]), sim_d, sim_p, rho_cruise, cfg)[0]
    T_cruise = thrust_in_flight(v_target, sim_rpm, sim_d, sim_p, rho_cruise, cfg)
    D_cruise = aero_drag(v_target, cfg)
    Fnet_cr  = T_cruise - D_cruise
    reached  = vs[-1] >= v_target * 0.98

    # Find cruise RPM (where thrust = drag at cruise speed with sim prop)
    cruise_rpm = find_cruise_rpm(sim_d, sim_p, rho_cruise, cfg)
    if np.isnan(cruise_rpm):
        print(f"  WARNING: prop {sim_d}x{sim_p} cannot sustain cruise at "
              f"{cfg['v_cruise']} m/s even at max RPM")
        cruise_rpm = cfg["rpm_max"]
    # Cruise fuel: shaft power from thrust requirement, not engine VE curve
    T_at_cruise_rpm = thrust_in_flight(cfg["v_cruise"], cruise_rpm, sim_d, sim_p, rho_cruise, cfg)
    P_shaft_cruise = shaft_power_from_thrust(
        T_at_cruise_rpm, cfg["v_cruise"], cruise_rpm, sim_d, sim_p, cfg)
    eta_cruise = prop_efficiency(cfg["v_cruise"], cruise_rpm, sim_d, sim_p, cfg)
    ff_cruise = float(fuel_flow(np.array([P_shaft_cruise]), cfg["bsfc"])[0])

    # ── Figure layout ───────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 18))
    fig.patch.set_facecolor("#f8f8f6")
    gs = gridspec.GridSpec(4, 2, figure=fig,
                           hspace=0.52, wspace=0.32,
                           left=0.07, right=0.96,
                           top=0.94, bottom=0.04)

    GRID_C  = "#e0dedd"
    TEXT_C  = "#333333"
    SPINE_C = "#cccccc"

    def style_ax(ax, title, xlabel, ylabel):
        ax.set_title(title, fontsize=11, fontweight="bold",
                     color=TEXT_C, pad=7, loc="left")
        ax.set_xlabel(xlabel, fontsize=9, color=TEXT_C)
        ax.set_ylabel(ylabel, fontsize=9, color=TEXT_C)
        ax.set_facecolor("#f8f8f6")
        ax.tick_params(colors=TEXT_C, labelsize=8)
        ax.grid(True, color=GRID_C, linewidth=0.6)
        for sp in ax.spines.values():
            sp.set_color(SPINE_C)

    # ── (0,0) Power & Torque ────────────────────────────────────────────
    anchor_kw  = cfg["rated_hp"] * 0.7457
    anchor_rpm = cfg["rated_rpm"]

    ax = fig.add_subplot(gs[0, 0])
    style_ax(ax, "Power & torque vs RPM", "RPM", "Power (kW)")
    ax.plot(rpms, pwr, color="#378ADD", linewidth=2, label="Power (kW)")
    ax2 = ax.twinx()
    ax2.plot(rpms, trq, color="#D85A30", linewidth=2, linestyle="--",
             label="Torque (Nm)")
    ax2.set_ylabel("Torque (Nm)", fontsize=9, color="#D85A30")
    ax2.tick_params(colors="#D85A30", labelsize=8)
    ax2.spines["right"].set_color(SPINE_C)
    ax.axvline(anchor_rpm, color="#378ADD", linewidth=0.8, linestyle=":")
    ax.annotate(f"{cfg['rated_hp']} HP @ {int(anchor_rpm)} RPM",
                xy=(anchor_rpm, anchor_kw),
                xytext=(anchor_rpm - 2200, anchor_kw * 0.85),
                fontsize=7.5, color="#378ADD",
                arrowprops=dict(arrowstyle="->", color="#378ADD", lw=0.8))
    lines1, lab1 = ax.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, lab1 + lab2, fontsize=7.5, loc="upper left")

    # ── (0,1) Static thrust - multiple props ────────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    style_ax(ax, "Static thrust vs RPM", "RPM", "Thrust (N)")
    for (d, p, col) in cfg["props"]:
        T = static_thrust(rpms, d, p, rho_anchor, cfg)
        ax.plot(rpms, T, color=col, linewidth=1.8, label=f"{d}x{p}")
    ax.scatter([cfg["anchor_thrust_rpm"]], [cfg["anchor_thrust_kg"] * G],
               color="#333333", zorder=5, s=50, marker="^", label="Mfr. anchor")
    ax.legend(fontsize=7.5, loc="upper left")
    ax.set_ylim(bottom=0)

    # ── (1,0) Fuel consumption ──────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    style_ax(ax, "Fuel consumption vs RPM", "RPM", "Fuel flow (ml/hr)")
    ax.plot(rpms, fuel, color="#BA7517", linewidth=2)
    ax.fill_between(rpms, fuel, alpha=0.12, color="#BA7517")
    ax.set_ylim(bottom=0)
    ax.axvline(cruise_rpm, color="#BA7517", linewidth=0.8, linestyle=":")
    ax.annotate(f"~{ff_cruise:.0f} ml/hr\n@ cruise ({cruise_rpm} RPM)",
                xy=(cruise_rpm, ff_cruise),
                xytext=(cruise_rpm + 700, ff_cruise * 0.6),
                fontsize=7.5, color="#BA7517",
                arrowprops=dict(arrowstyle="->", color="#BA7517", lw=0.8))

    # ── (1,1) Endurance vs tank ─────────────────────────────────────────
    tanks = np.linspace(100, 2000, 300)
    endur = tanks / ff_cruise
    target_hr = cfg["endurance_target_hr"]
    tank_target = ff_cruise * target_hr

    ax = fig.add_subplot(gs[1, 1])
    style_ax(ax, f"Endurance vs tank size  (cruise @ {cruise_rpm} RPM)",
             "Tank size (ml)", "Endurance (hr)")
    ax.plot(tanks, endur, color="#7F77DD", linewidth=2)
    ax.fill_between(tanks, endur, alpha=0.12, color="#7F77DD")
    ax.axhline(target_hr, color="#1D9E75", linewidth=0.9, linestyle="--")
    ax.annotate(f"{target_hr:.0f} hr target", xy=(200, target_hr),
                xytext=(250, target_hr + 0.15),
                fontsize=7.5, color="#1D9E75")
    ax.set_ylim(bottom=0)
    ax.axvline(tank_target, color="#7F77DD", linewidth=0.8, linestyle=":")
    ax.annotate(f"{tank_target:.0f} ml\nfor {target_hr:.0f} hr",
                xy=(tank_target, target_hr),
                xytext=(tank_target + 100, target_hr * 0.75),
                fontsize=7.5, color="#7F77DD",
                arrowprops=dict(arrowstyle="->", color="#7F77DD", lw=0.8))

    # ── (2,0) Force balance vs speed ────────────────────────────────────
    prop_label = f"{sim_d}x{sim_p}"
    ax = fig.add_subplot(gs[2, 0])
    style_ax(ax, f"Force balance vs airspeed  ({sim_rpm} RPM, {prop_label})",
             "Airspeed (m/s)", "Force (N)")
    ax.plot(v_range, T_vs_v, color="#1D9E75", linewidth=2, label="Thrust")
    ax.plot(v_range, D_vs_v, color="#D85A30", linewidth=2, label="Drag")
    ax.plot(v_range, F_net,  color="#378ADD", linewidth=1.5,
            linestyle="--", label="Net force")
    ax.axhline(0, color=SPINE_C, linewidth=0.8)
    ax.axvline(v_target, color="#BA7517", linewidth=0.8, linestyle=":",
               label=f"Cruise ({v_target} m/s)")
    ax.annotate(f"T={T_cruise:.1f} N  D={D_cruise:.1f} N\n"
                f"Fnet={Fnet_cr:+.1f} N @ cruise",
                xy=(v_target, T_cruise),
                xytext=(v_target - 14, T_cruise + 15),
                fontsize=7.5, color=TEXT_C,
                arrowprops=dict(arrowstyle="->", color=TEXT_C, lw=0.7))
    ax.legend(fontsize=7.5)
    ax.set_xlim(0, 42)

    # ── (2,1) Speed vs time ─────────────────────────────────────────────
    ax = fig.add_subplot(gs[2, 1])
    lbl = (f"Reached {v_target} m/s in {t_total:.1f} s"
           if reached else "Target not reached")
    style_ax(ax, f"Airspeed vs time - {lbl}", "Time (s)", "Airspeed (m/s)")
    ax.plot(ts, vs, color="#7F77DD", linewidth=2)
    ax.fill_between(ts, vs, alpha=0.1, color="#7F77DD")
    ax.axhline(v_target, color="#BA7517", linewidth=0.9,
               linestyle="--", label=f"Target {v_target} m/s")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7.5)
    if reached:
        ax.annotate(f"t = {t_total:.1f} s", xy=(t_total, v_target),
                    xytext=(t_total * 0.6, v_target * 0.75),
                    fontsize=8, color="#7F77DD",
                    arrowprops=dict(arrowstyle="->", color="#7F77DD", lw=0.8))

    # ── (3,0) Acceleration vs speed ─────────────────────────────────────
    ax = fig.add_subplot(gs[3, 0])
    style_ax(ax, "Acceleration vs airspeed",
             "Airspeed (m/s)", "Acceleration (m/s^2)")
    ax.plot(vs, accs, color="#D85A30", linewidth=2)
    ax.fill_between(vs, accs, alpha=0.1, color="#D85A30")
    ax.axhline(0, color=SPINE_C, linewidth=0.8)
    ax.axvline(v_target, color="#BA7517", linewidth=0.8, linestyle=":")

    # ── (3,1) Distance vs time ──────────────────────────────────────────
    ax = fig.add_subplot(gs[3, 1])
    style_ax(ax, f"Ground distance covered  ({x_total:.0f} m to reach cruise)",
             "Time (s)", "Distance (m)")
    ax.plot(ts, xs, color="#1D9E75", linewidth=2)
    ax.fill_between(ts, xs, alpha=0.1, color="#1D9E75")
    if reached:
        ax.axhline(x_total, color="#1D9E75", linewidth=0.8, linestyle=":")
        ax.annotate(f"{x_total:.0f} m", xy=(t_total, x_total),
                    xytext=(t_total * 0.4, x_total * 0.7),
                    fontsize=8, color="#1D9E75",
                    arrowprops=dict(arrowstyle="->", color="#1D9E75", lw=0.8))

    # ── Super-title ─────────────────────────────────────────────────────
    fig.suptitle(
        f"RCGF Stinger 35cc RE - Engine & Transition Model   "
        f"[MTOW {cfg['mass_kg']:.0f} kg | {sim_rpm} RPM | {prop_label} prop | "
        f"drag {D_cruise:.1f} N @ {v_target:.0f} m/s | eta_pitch {cfg['eta_pitch']}]",
        fontsize=11, fontweight="bold", color=TEXT_C, y=0.97
    )

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "spearhead_engine_model.png")
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved -> {out}")

    # ── Console summary ─────────────────────────────────────────────────
    print(f"\n--- Transition summary ({prop_label} @ {sim_rpm} RPM) ---")
    print(f"  Anchor altitude       : {cfg['anchor_alt_m']:.0f} m  "
          f"(rho = {rho_anchor:.4f} kg/m^3)")
    print(f"  Cruise altitude       : {cfg['cruise_alt_m']:.0f} m  "
          f"(rho = {rho_cruise:.4f} kg/m^3)")
    print(f"  Static thrust ({sim_rpm} RPM): {T0_val:.1f} N  ({T0_val/G:.2f} kg)")
    print(f"  Thrust @ {v_target} m/s ({sim_rpm}): {T_cruise:.1f} N  ({T_cruise/G:.2f} kg)")
    print(f"  Drag   @ {v_target} m/s     : {D_cruise:.1f} N")
    print(f"  Net force ({sim_rpm} RPM)   : {Fnet_cr:+.1f} N")
    if not np.isnan(cruise_rpm):
        T_at_cr = thrust_in_flight(v_target, cruise_rpm, sim_d, sim_p, rho_cruise, cfg)
        print(f"  Cruise RPM            : {cruise_rpm:.0f}  (thrust = drag = {T_at_cr:.1f} N)")
    bat_Wh = cfg["battery_Wh"]
    # Thermal efficiency: eta_th = 3600 / (BSFC [g/kWh] * LHV [kJ/g])
    # Gasoline LHV = 43.4 MJ/kg = 43.4 kJ/g
    LHV_kJ_per_g = 43.4
    eta_thermal = 3600.0 / (cfg["bsfc"] * LHV_kJ_per_g) * 100  # percent
    eta_gen = cfg["eta_generator"]
    print(f"  Thermal efficiency   : {eta_thermal:.1f}%  (BSFC {cfg['bsfc']:.0f} g/kWh)")
    print(f"  Fuel flow (cruise)   : {ff_cruise:.0f} ml/hr  (shaft {P_shaft_cruise:.2f} kW)")
    # Shaft headroom at cruise RPM (available for gen/climb)
    P_avail_cr_kw = float(power_curve(np.array([cruise_rpm]), cfg)[0])
    P_headroom_kw = P_avail_cr_kw - P_shaft_cruise  # max extra shaft power
    P_excess_w = P_headroom_kw * 1000  # watts
    ROC = P_excess_w / (cfg["mass_kg"] * G)
    print(f"  Excess power (cruise): {P_excess_w:.0f} W  ({P_excess_w/1000:.2f} kW)")
    print(f"  Max rate of climb    : {ROC:.1f} m/s  ({ROC*60:.0f} ft/min)")
    if P_headroom_kw > 0:
        P_elec_kw = P_headroom_kw * eta_gen
        t_charge_hr = bat_Wh / (P_elec_kw * 1000)
        P_total_kw = P_avail_cr_kw  # WOT at cruise RPM
        ff_total = float(fuel_flow(np.array([P_total_kw]), cfg["bsfc"])[0])
        extra_fuel_ml = (ff_total - ff_cruise) * t_charge_hr
        print(f"  Fuel flow (charging) : {ff_total:.0f} ml/hr  "
              f"(shaft {P_total_kw:.2f} kW = {P_shaft_cruise:.2f} cruise "
              f"+ {P_headroom_kw:.2f} gen)")
        print(f"  Battery charge       : {t_charge_hr:.1f} hr  |  "
              f"+{extra_fuel_ml:.0f} ml extra fuel  "
              f"({bat_Wh:.0f} Wh @ {P_elec_kw*1000:.0f} W elec, "
              f"gen {eta_gen*100:.0f}%)")
    else:
        print(f"  Battery charge       : N/A  (no headroom at cruise RPM)")
    if reached:
        print(f"  Transition time       : {t_total:.1f} s")
    else:
        print(f"  Target speed NOT reached in {cfg['t_max']:.0f} s")
    print(f"  Distance covered      : {x_total:.0f} m")
    # ── Prop comparison table ─────────────────────────────────────────
    v_cr = cfg["v_cruise"]
    props = cfg["props"]
    print(f"\n--- Prop comparison @ {v_cr:.0f} m/s cruise, {cfg['cruise_alt_m']:.0f} m alt ---")
    bat_Wh = cfg["battery_Wh"]
    print(f"  {'Prop':<10}{'RPM':>7}{'T_stat':>8}{'T_cr':>7}"
          f"{'eta':>6}{'P_shaft':>8}{'Q_req':>7}{'Q_avail':>8}"
          f"{'Margin':>8}{'P_exc':>7}{'ROC':>7}"
          f"{'Fuel':>8}{'Tank':>10}{'t_chg':>7}{'chg_ml':>8}")
    print(f"  {'':10}{'':>7}{'(N)':>8}{'(N)':>7}"
          f"{'(%)':>6}{'(kW)':>8}{'(Nm)':>7}{'(Nm)':>8}"
          f"{'(%)':>8}{'(W)':>7}{'(m/s)':>7}"
          f"{'(ml/h)':>8}{'({:.0f}hr ml)':>10}{'(hr)':>7}{'(ml)':>8}".format(target_hr))
    print(f"  {'-'*118}")

    for (d, p, _col) in props:
        label = f"{d}x{p}"
        cr_rpm = find_cruise_rpm(d, p, rho_cruise, cfg)
        if np.isnan(cr_rpm):
            print(f"  {label:<10}  -- cannot sustain cruise --")
            continue
        T0  = static_thrust(np.array([cr_rpm]), d, p, rho_cruise, cfg)[0]
        T_cr = thrust_in_flight(v_cr, cr_rpm, d, p, rho_cruise, cfg)
        eta = prop_efficiency(v_cr, cr_rpm, d, p, cfg)
        P_sh = shaft_power_from_thrust(T_cr, v_cr, cr_rpm, d, p, cfg)
        # Required torque = shaft power / omega
        omega_cr = cr_rpm * (2 * np.pi / 60)
        Q_req = (P_sh * 1000) / omega_cr  # Nm
        # Available torque from engine at this RPM (WOT)
        P_avail = float(power_curve(np.array([cr_rpm]), cfg)[0])
        Q_avail = (P_avail * 1000) / omega_cr  # Nm
        Q_margin = (Q_avail - Q_req) / Q_avail * 100  # % headroom
        # Excess power = shaft headroom at cruise RPM (available for gen/climb)
        P_exc = (P_avail - P_sh) * 1000  # watts
        roc = P_exc / (cfg["mass_kg"] * G)  # m/s
        ff  = float(fuel_flow(np.array([P_sh]), cfg["bsfc"])[0])
        tnk = ff * target_hr
        # Battery charging: generator uses all available headroom
        P_headroom = P_avail - P_sh  # kW available above cruise load
        P_elec_kw = P_headroom * eta_gen if P_headroom > 0 else 0
        t_chg = bat_Wh / (P_elec_kw * 1000) if P_elec_kw > 0 else float("inf")
        ff_chg = float(fuel_flow(np.array([P_avail]), cfg["bsfc"])[0]) if P_headroom > 0 else ff
        extra_fuel = (ff_chg - ff) * t_chg if t_chg < 99 else float("inf")
        t_chg_s = f"{t_chg:.1f}" if t_chg < 99 else "N/A"
        chg_ml_s = f"+{extra_fuel:.0f}" if extra_fuel < 1e6 else "N/A"
        print(f"  {label:<10}{cr_rpm:>7.0f}{T0:>8.1f}{T_cr:>7.1f}"
              f"{eta*100:>6.1f}{P_sh:>8.2f}{Q_req:>7.2f}{Q_avail:>8.2f}"
              f"{Q_margin:>7.1f}%{P_exc:>7.0f}{roc:>7.1f}"
              f"{ff:>8.0f}{tnk:>10.0f}{t_chg_s:>7}{chg_ml_s:>8}")


if __name__ == "__main__":
    make_plots()
