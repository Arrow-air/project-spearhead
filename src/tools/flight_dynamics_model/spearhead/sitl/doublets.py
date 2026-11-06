"""Finite physical-input doublets, with firmware-controlled recovery."""

import json
from pathlib import Path

import numpy as np

DOUBLETS = {
    "pitch_doublet": {"surface": "de", "channel": "elevator", "axis": 1, "amplitude": 3., "pulse_s": .75, "unit": "deg"},
    "roll_doublet": {"surface": "da", "channel": "aileron", "axis": 0, "amplitude": 5., "pulse_s": .5, "unit": "deg"},
    "yaw_doublet": {"surface": "dr", "channel": "rudder", "axis": 2, "amplitude": 5., "pulse_s": .75, "unit": "deg"},
    "throttle_doublet": {"surface": "throttle", "channel": "pusher", "axis": None, "amplitude": .10, "pulse_s": 2., "unit": "fraction"},
}


def doublet_plan(name, trim, mapping):
    plan = dict(DOUBLETS[name], start_s=10., name=name)
    baseline = getattr(trim, plan["surface"])
    if plan["unit"] == "deg":
        headroom = getattr(mapping, plan["channel"]).travel_deg-abs(np.rad2deg(baseline))
        if headroom < plan["amplitude"]:
            raise ValueError("Doublet does not fit around the selected surface trim")
    else:
        plan["amplitude"] = min(plan["amplitude"], .9*min(baseline, 1-baseline))
        if plan["amplitude"] < .03:
            raise ValueError("Throttle trim leaves less than 3 percentage points of symmetric doublet authority")
    plan["trim_value"] = float(np.rad2deg(baseline) if plan["unit"] == "deg" else baseline)
    plan["recovery_at_s"] = plan["start_s"]+2*plan["pulse_s"]+1
    return plan


def doublet_delta(plan, elapsed):
    if plan["start_s"] <= elapsed < plan["start_s"]+plan["pulse_s"]:
        return plan["amplitude"]
    if plan["start_s"]+plan["pulse_s"] <= elapsed < plan["start_s"]+2*plan["pulse_s"]:
        return -plan["amplitude"]
    return 0.


def doublet_mode(plan, elapsed):
    return 0 if plan["start_s"]-1 <= elapsed < plan["recovery_at_s"] else 5


def analyse_doublet(path, initial, plan, duration, release_at, *, telemetry_path=None):
    from .transition import confirmed_mode_samples
    rows = [json.loads(line) for line in Path(path).open() if line.strip()]
    rows = [row for row in rows if not row["hold"]]
    if not rows or release_at is None:
        return {"passed": False, "test_kind": "control_doublet", "reason": "No free-flight doublet data", "doublet_plan": plan}
    times = np.array([r["time_s"] for r in rows]); t = times-release_at
    states = np.array([r["state"] for r in rows])
    actual = np.array([r["control"][plan["surface"]] for r in rows])
    if plan["unit"] == "deg":
        actual = np.rad2deg(actual)
    delta = actual-plan["trim_value"]
    manual, _ = confirmed_mode_samples(times, telemetry_path or Path(path).with_name("telemetry.jsonl"), 0)
    fbwa, _ = confirmed_mode_samples(times, telemetry_path or Path(path).with_name("telemetry.jsonl"), 5)
    pulses = []
    tolerance = .12 if plan["unit"] == "deg" else .003
    for index, sign in enumerate((1, -1)):
        start = plan["start_s"]+index*plan["pulse_s"]; end = start+plan["pulse_s"]
        selected = (t >= start+.2) & (t < end-.1)
        measured = float(np.median(delta[selected])) if selected.any() else None
        delivered = bool(selected.sum() >= 10 and manual[selected].all() and
                         np.max(abs(delta[selected]-sign*plan["amplitude"])) <= tolerance)
        pulses.append({"start_s": start, "end_s": end, "requested": sign*plan["amplitude"], "measured": measured, "delivered": delivered})
    neutral = (t >= plan["start_s"]+2*plan["pulse_s"]+.2) & (t < plan["recovery_at_s"]-.1)
    restored = bool(neutral.sum() >= 10 and manual[neutral].all() and np.max(abs(delta[neutral])) <= tolerance)
    attitude = np.rad2deg(states[:, 6:9]-initial[6:9]); attitude[:, 2] = np.rad2deg(np.unwrap(states[:, 8])-initial[8])
    rates = np.rad2deg(states[:, 9:12]); speed = np.linalg.norm(states[:, 3:6], axis=1)
    tail = t >= max(plan["recovery_at_s"]+5, duration-5)
    recovery = bool(tail.any() and np.max(abs(attitude[tail, :2])) <= 2 and np.max(abs(rates[tail])) <= 1 and
                    np.max(abs(speed[tail]-np.linalg.norm(initial[3:6]))) <= 1.5)
    checks = {"complete_duration": bool(t[-1] >= duration-.05), "both_inputs_delivered": all(p["delivered"] for p in pulses),
              "confirmed_recovery_mode": bool(fbwa[-1]),
              "neutral_input_restored": restored, "bounded_attitude": bool(np.max(abs(np.rad2deg(states[:, 6]))) <= 20 and
                                                                           np.max(abs(np.rad2deg(states[:, 7]))) <= 15)}
    return {"passed": all(checks.values()), "test_kind": "control_doublet", "experiment_checks": checks,
            "doublet_plan": plan, "pulses": pulses, "free_flight_s": float(t[-1]),
            "peak_attitude_change_deg": np.max(abs(attitude), axis=0).tolist(),
            "peak_rate_deg_s": np.max(abs(rates), axis=0).tolist(), "speed_range_m_s": [float(speed.min()), float(speed.max())],
            "max_altitude_deviation_m": float(np.max(abs(-states[:, 2]+initial[2]))), "returned_near_trim": recovery,
            "outside_control_coverage_fraction": float(np.mean([r.get("control_coverage_exceeded", False) for r in rows])),
            "interpretation": "Input experiment checks; recovery is reported separately. Yaw heading need not return. MANUAL excitation followed by frozen FBWA recovery."}
