"""Hover pilot doublets; firmware owns attitude control and motor mixing."""

import json
from pathlib import Path

import numpy as np

from .doublets import doublet_delta
from .transition import confirmed_mode_samples

HOVER_DOUBLETS = {
    "pitch_doublet": {"channel": 2, "axis": 1, "amplitude": 8., "unit": "deg", "label": "Pitch attitude"},
    "roll_doublet": {"channel": 1, "axis": 0, "amplitude": 8., "unit": "deg", "label": "Roll attitude"},
    "yaw_doublet": {"channel": 4, "axis": 2, "amplitude": .2, "unit": "fraction", "label": "Yaw stick"},
    "throttle_doublet": {"channel": 3, "axis": None, "amplitude": .1, "unit": "fraction", "label": "Pilot lift collective"},
}


def hover_doublet_plan(name, params):
    plan = dict(HOVER_DOUBLETS[name], name=name, regime="hover", start_s=10., pulse_s=2., recovery_at_s=15.)
    plan["hover_pwm"] = round(1060+940*params["Q_M_THST_HOVER"])
    plan["trim_value"] = params["Q_M_THST_HOVER"] if name == "throttle_doublet" else 0.
    if plan["unit"] == "deg":
        angle_max = params["Q_A_ANGLE_MAX"]
        if not 8 <= angle_max <= 20:
            raise ValueError("Hover attitude doublet requires the frozen 8–20 degree lean limit")
        plan["pwm_per_unit"] = 500/angle_max
    else:
        plan["pwm_per_unit"] = 940 if name == "throttle_doublet" else 500
    if name == "throttle_doublet" and not .1 < plan["trim_value"] < .9:
        raise ValueError("Hover collective leaves insufficient symmetric doublet headroom")
    plan["neutral_pwm"] = round(1060+940*plan["trim_value"]) if name == "throttle_doublet" else 1500
    plan["excitation_mode"] = 17 if name == "throttle_doublet" else 18
    plan["recovery_mode"] = 18
    return plan


def hover_doublet_command(plan, elapsed, confirmed_mode):
    mode = plan["excitation_mode"] if plan["start_s"]-1 <= elapsed < plan["recovery_at_s"] else 18
    rc = [1500]*8
    if confirmed_mode == 17:
        rc[2] = plan["hover_pwm"]
    # Never deliver excitation in an unconfirmed mode. Reports check actual RC readback.
    if mode == confirmed_mode and plan["start_s"]-1 <= elapsed < plan["recovery_at_s"]:
        rc[plan["channel"]-1] = round(plan["neutral_pwm"]+plan["pwm_per_unit"]*doublet_delta(plan, elapsed))
    return mode, rc


def analyse_hover_doublet(path, release_at, duration, initial, plan):
    rows = [json.loads(line) for line in Path(path).open() if line.strip()]
    rows = [r for r in rows if not r["hold"]]
    base = {"passed": False, "test_kind": "control_doublet", "doublet_plan": plan}
    if not rows or release_at is None:
        return dict(base, reason="No free hover doublet frames")
    times = np.array([r["time_s"] for r in rows]); t = times-release_at
    x = np.array([r["state"] for r in rows]); pwm = np.array([r["pwm"][4:8] for r in rows])
    telemetry = Path(path).with_name("telemetry.jsonl")
    records = [json.loads(line) for line in telemetry.open() if line.strip()] if telemetry.exists() else []
    rc = [(r["fdm_time_s"], r["message"].get(f"chan{plan['channel']}_raw")) for r in records
          if r["message"].get("mavpackettype") == "RC_CHANNELS"]
    rc = [(time, value) for time, value in rc if value is not None and 900 <= value <= 2100]
    rt = np.array([time-release_at for time, _ in rc]); rv = np.array([value for _, value in rc])
    excitation, _ = confirmed_mode_samples(times, telemetry, plan["excitation_mode"])
    recovered, _ = confirmed_mode_samples(times, telemetry, 18)
    pulses = []
    for index, sign in enumerate((1, -1)):
        start = plan["start_s"]+index*plan["pulse_s"]; end = start+plan["pulse_s"]
        selected = (rt >= start+.3) & (rt < end-.1)
        physical = (t >= start+.3) & (t < end-.1)
        target = round(plan["neutral_pwm"]+sign*plan["amplitude"]*plan["pwm_per_unit"])
        delivered = bool(selected.sum() >= 3 and physical.any() and excitation[physical].all() and
                         np.max(abs(rv[selected]-target)) <= 2)
        measured = float((np.median(rv[selected])-plan["neutral_pwm"])/plan["pwm_per_unit"]) if selected.any() else None
        pulses.append({"start_s": start, "end_s": end, "requested": sign*plan["amplitude"],
                       "measured": measured, "rc_pwm_us": target, "delivered": delivered})
    neutral = (rt >= 14.3) & (rt < 14.9)
    restored = bool(neutral.sum() >= 3 and np.max(abs(rv[neutral]-plan["neutral_pwm"])) <= 2)
    attitude = np.rad2deg(x[:, 6:9]-np.asarray(initial)[6:9]); attitude[:, 2] = np.rad2deg(np.unwrap(x[:, 8])-initial[8])
    rates = np.rad2deg(x[:, 9:12]); tail = t >= duration-5
    # QHOVER holds altitude and attitude, not horizontal position or original heading.
    near_trim = bool(tail.any() and np.max(abs(np.rad2deg(x[tail, 6:8]))) < 2 and np.max(abs(rates[tail])) < 1)
    checks = {"complete_duration": bool(t[-1] >= duration-.05), "both_inputs_delivered": all(p["delivered"] for p in pulses),
              "neutral_input_restored": restored, "confirmed_recovery_mode": bool(recovered[-1]),
              "bounded_attitude": bool(np.max(abs(np.rad2deg(x[:, 6:8]))) <= 20),
              "motor_range": bool(np.all((pwm >= 1100) & (pwm <= 1940))),
              "pusher_off": bool(max(r["control"]["throttle"] for r in rows) <= .01)}
    return dict(base, passed=all(checks.values()), experiment_checks=checks, pulses=pulses,
                free_flight_s=float(t[-1]), peak_attitude_change_deg=np.max(abs(attitude), axis=0).tolist(),
                peak_rate_deg_s=np.max(abs(rates), axis=0).tolist(),
                max_altitude_deviation_m=float(np.max(abs(-x[:, 2]+initial[2]))),
                final_altitude_offset_m=float(-x[-1, 2]+initial[2]),
                horizontal_displacement_m=float(np.linalg.norm(x[-1, :2]-x[0, :2])),
                lift_pwm_range_us=[int(pwm.min()), int(pwm.max())], returned_near_trim=near_trim,
                interpretation="Actual RC_CHANNELS proves pilot input delivery; FDM states and lift PWM record the response. QHOVER recovery holds altitude/attitude, not position or original yaw. Collective uses QSTABILIZE and the existing motor linearization; no controller tuning.")
