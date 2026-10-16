"""CG consequences for one aircraft: static margin, elevator trim and hover headroom.

Static margin is the stick-fixed value from the AeroDB at the solved trim:
SM = -dCm/dCL about the CG for a small angle-of-attack change at constant speed,
using wing, tail and fuselage loads only. Trim uses the six-axis solver with the
configured surface travel. Hover headroom is the most loaded rotor's thrust
against its thrust at the firmware's maximum spin (Q_M_SPIN_MAX 0.95).
"""

from dataclasses import replace
from functools import lru_cache
import math

import numpy as np

from .aircraft import load_aircraft_config
from .force_moment import compute_force_moment_breakdown
from .propulsion import hover_trim_control, load_static_performance

TRIM_SPEEDS_M_S = (20., 22., 24.)
MEASURED_ELEVATOR_DEG = 15.
SPIN_MAX = .95


def static_margin(params, state, control, step_rad=1e-3):
    speed = float(np.linalg.norm(state[3:6]))
    alpha0 = math.atan2(state[5], state[3])
    q = .5*params.rho*speed**2
    points = []
    for delta in (-step_rad, step_rad):
        x = np.array(state, dtype=float)
        alpha = alpha0 + delta
        x[3], x[5] = speed*math.cos(alpha), speed*math.sin(alpha)
        x[9:12] = 0
        loads = compute_force_moment_breakdown(0, x, params, control)
        force = sum(loads[k]["force_b"] for k in ("wing", "tail", "fuselage"))
        moment = sum(loads[k]["moment_b"] for k in ("wing", "tail", "fuselage"))
        lift = -force[2]*math.cos(alpha) + force[0]*math.sin(alpha)
        points.append((lift/(q*params.Sref), moment[1]/(q*params.Sref*params.cref)))
    (cl1, cm1), (cl2, cm2) = points
    return -(cm2 - cm1)/(cl2 - cl1)


def hover_headroom(params):
    control = hover_trim_control(params)
    curve = load_static_performance(params.lift_system.performance_path)
    voltage = params.battery.voltage_v
    maximum = curve.lookup(voltage, SPIN_MAX)["thrust_n"]
    loaded = max(curve.lookup(voltage, command)["thrust_n"] for command in control.lift_motor_throttles)
    return loaded/maximum


@lru_cache(maxsize=64)
def cg_card(aircraft_id, cg_x_m, cg_z_m, travel_deg=30., pusher=None):
    """Return a JSON-friendly summary for one CG and pusher option; cached because trims take seconds."""
    from .propulsion import pusher_performance
    from .sitl.physics import full_trim
    params = replace(load_aircraft_config(aircraft_id, pusher=pusher), cg_body_xyz=np.array([cg_x_m, 0., cg_z_m]),
                     control_deflection_policy="linear", control_surface_limits_deg=(travel_deg,)*3)
    trims = {}
    margin = None
    for speed in TRIM_SPEEDS_M_S:
        try:
            state, control = full_trim(params, speed, 100.)
        except ValueError:
            trims[speed] = None
            continue
        trims[speed] = {"elevator_deg": math.degrees(control.de), "pusher_command": control.throttle,
                        "beyond_measured": abs(math.degrees(control.de)) > MEASURED_ELEVATOR_DEG}
        if params.pusher_system is not None:
            performance = pusher_performance(control, params, float(state[3]))
            trims[speed].update(pusher_current_a=performance["current_a"],
                                above_current_limit=bool(performance.get("above_sustainable_current", False)))
        if speed == 22.:
            margin = static_margin(params, state, control)
    try:
        headroom = hover_headroom(params) if params.lift_system is not None else None
    except ValueError:
        headroom = None
    system = params.pusher_system
    return {"cg_body_xyz_m": [cg_x_m, 0., cg_z_m], "mass_kg": params.mass, "static_margin_22": margin,
            "pusher": None if system is None else {"option": system.option, "label": system.label,
                                                   "current_limit_a": system.sustainable_current_a,
                                                   "current_limit_label": system.current_limit_label},
            "trims": trims, "hover_max_rotor_fraction": headroom,
            "rotor_centroid_x_m": None if params.lift_system is None else
            float(np.mean([p[0] for p in params.lift_system.positions_body_m]))}


def cg_card_text(card):
    parts = []
    margin = card["static_margin_22"]
    parts.append("SM at 22 m/s: " + ("no trim" if margin is None else f"{100*margin:+.1f} %"))
    trims = []
    for speed, trim in card["trims"].items():
        if trim is None:
            trims.append(f"{speed:g} m/s no trim")
        else:
            trims.append(f"{speed:g} m/s {trim['elevator_deg']:+.1f}°" + (" *" if trim["beyond_measured"] else ""))
    parts.append("Elevator trim: " + ", ".join(trims))
    cruise = card["trims"].get(22.)
    if card.get("pusher") and cruise and "pusher_current_a" in cruise:
        parts.append(f"Pusher at 22 m/s: {100*cruise['pusher_command']:.0f} % · {cruise['pusher_current_a']:.1f} A"
                     + (f" (above {card['pusher']['current_limit_a']:g} A {card['pusher']['current_limit_label']})"
                        if cruise["above_current_limit"] else ""))
    if card["hover_max_rotor_fraction"] is not None:
        parts.append(f"Hover: busiest rotor at {100*card['hover_max_rotor_fraction']:.0f} % of max thrust")
    text = " · ".join(parts)
    if any(t and t["beyond_measured"] for t in card["trims"].values()):
        text += f" · * beyond the measured ±{MEASURED_ELEVATOR_DEG:g}° elevator data"
    return text
