"""Acceptance metrics from complete hold-free physics records, in degrees."""

import json
from pathlib import Path

import numpy as np


STEPS = [(10+18*i, axis, amplitude) for i, (axis, amplitude) in enumerate(
    [(1, .2), (1, -.2), (1, .5), (1, -.5), (0, .5), (0, -.5)])]
ELEVATOR_PULSES = [(10., 11.5, 1.), (11.5, 13., -1.)]


def elevator_at(elapsed):
    """Physical elevator increment in degrees; MANUAL, down then up."""
    for start, end, delta in ELEVATOR_PULSES:
        if start <= elapsed < end:
            return delta
    return 0.


def command_at(elapsed):
    """RC increments with exact quantized targets; 12 s steps, 6 s returns."""
    for start, axis, requested in STEPS:
        if start <= elapsed < start+12:
            limit = 30 if axis == 0 else 10
            delta_pwm = round(requested / limit * 500)
            return axis, delta_pwm, delta_pwm/500*limit
    return None, 0, 0.


def step_response(t, values, target, amplitude):
    t, values = np.asarray(t), np.asarray(values)
    if (t.ndim != 1 or values.shape != t.shape or len(t) < 2 or not np.isfinite(t).all()
            or not np.isfinite(values).all() or np.any(np.diff(t) <= 0)
            or not np.isfinite([target, amplitude]).all() or amplitude == 0):
        raise ValueError("Step response requires finite ordered samples and nonzero amplitude")
    band = max(.1, abs(amplitude)*.1)
    outside = np.flatnonzero(abs(values-target) > band)
    settling = 0. if not len(outside) else (None if outside[-1] == len(t)-1 else float(t[outside[-1]+1]))
    overshoot = max(0., float(np.max(np.sign(amplitude)*(values-target)))) / abs(amplitude)*100
    return {"target_deg": float(target), "amplitude_deg": float(amplitude), "band_deg": band,
            "settling_s": settling, "overshoot_percent": overshoot,
            "passed": settling is not None and settling < 5 and overshoot < 30}


def analyse_truth(path, initial, scenario, duration, perturbation=None):
    rows = [json.loads(line) for line in Path(path).open()]
    free = [r for r in rows if not r["hold"]]
    if not free:
        raise ValueError("No hold-free physics frames")
    t = np.array([r["time_s"] for r in free])
    t -= t[0]
    state = np.array([r["state"] for r in free])
    pwm = np.array([r["pwm"][:4] for r in free])
    attitude = np.rad2deg(state[:, 6:8])
    target = np.rad2deg(initial[6:8]) if scenario == "manual" else np.array([0., np.trunc(np.rad2deg(initial[7])*100)/100])
    target_series = np.tile(target, (len(t), 1))
    if scenario == "steps":
        for i, ti in enumerate(t):
            axis, _, amplitude = command_at(ti)
            if axis is not None:
                target_series[i, axis] += amplitude
    error = attitude-target_series
    saturated = np.any((pwm <= 1000) | (pwm >= 2000), axis=1)
    run_lengths = []
    count = 0
    dt = float(np.median(np.diff(t)))
    for flag in saturated:
        if flag:
            count += 1
        else:
            run_lengths.append(count*dt)
            count = 0
    run_lengths.append(count*dt)
    speed = np.linalg.norm(state[:, 3:6], axis=1)
    rms = np.sqrt(np.mean(error**2, axis=0))
    steps = []
    recovery = []
    if perturbation is not None and any(perturbation):
        for axis in range(2):
            outside = np.flatnonzero(abs(error[:, axis]) > .1)
            recovery.append(0. if not len(outside) else
                            (None if outside[-1] == len(t)-1 else float(t[outside[-1]+1])))
    if scenario == "steps":
        for start, axis, _ in STEPS:
            selected = (t >= start) & (t < start+12)
            if not np.any(selected) or t[selected][-1] < start+11.9:
                raise ValueError("Incomplete step test")
            _, _, amplitude = command_at(start)
            response = step_response(t[selected]-start, attitude[selected, axis], target[axis]+amplitude, amplitude)
            steps.append({"start_s": start, "axis": ["roll", "pitch"][axis], **response})
    return {"passed": bool(t[-1] >= duration-.05 and max(rms) < .5 and max(run_lengths) < 1
                            and all(step["passed"] for step in steps)
                            and all(value is not None and value < 5 for value in recovery)),
            "free_flight_s": float(t[-1]), "physics_frames": len(free),
            "target_attitude_deg": target.tolist(), "attitude_rms_error_deg": rms.tolist(),
            "attitude_peak_error_deg": np.max(abs(error), axis=0).tolist(),
            "final_state": state[-1].tolist(), "speed_range_m_s": [float(min(speed)), float(max(speed))],
            "release_altitude_m": float(-state[0, 2]),
            "max_altitude_deviation_m": float(np.max(abs(state[:, 2]-state[0, 2]))),
            "max_altitude_deviation_from_air_start_m": float(np.max(abs(state[:, 2]-initial[2]))),
            "saturation_fraction": float(np.mean(saturated)), "max_saturation_s": max(run_lengths),
            "alpha_range_deg": [min(r["alpha_deg"] for r in free), max(r["alpha_deg"] for r in free)],
            "beta_range_deg": [min(r["beta_deg"] for r in free), max(r["beta_deg"] for r in free)],
            "outside_control_coverage_fraction": float(np.mean([r["control_coverage_exceeded"] for r in free])),
            "extrapolated_deflection_fraction": float(np.mean([bool((r.get("control_coverage") or {}).get("extrapolated_surfaces")) for r in free])),
            "disturbance_recovery_s": recovery, "steps": steps}


def analyse_elevator_truth(path, initial, trim, duration, *, surface_limit_deg=5.):
    """Separate input-demonstration checks, not nominal FBWA tracking criteria."""
    report = analyse_truth(path, initial, "manual", duration)
    rows = [json.loads(line) for line in Path(path).open()]
    free = [r for r in rows if not r["hold"]]
    t = np.array([r["time_s"] for r in free])
    t -= t[0]
    elevator = np.rad2deg([r["control"]["de"] for r in free])
    baseline = float(np.rad2deg(trim.de))
    pitch = np.rad2deg([r["state"][7] for r in free])
    pulses = []
    for start, end, delta in ELEVATOR_PULSES:
        # RC updates at 20 Hz; assess the interior, excluding command latency.
        selected = (t >= start+.2) & (t <= end-.1)
        measured = None if not selected.any() else float(np.median(elevator[selected])-baseline)
        ok = measured is not None and abs(measured-delta) < .03 and selected.sum() >= 100
        pulses.append({"start_s": start, "end_s": end, "requested_delta_deg": delta,
                       "measured_delta_deg": measured, "passed": bool(ok)})
    after = t >= duration-5
    final_error = float(np.max(abs(pitch[after]-np.rad2deg(initial[7])))) if after.any() else None
    excursion = float(np.max(abs(pitch-np.rad2deg(initial[7]))))
    neutral = ((t < 9.8) | (t > 13.2))
    checks = {"complete_duration": report["free_flight_s"] >= duration-.05,
              "both_physical_pulses": all(p["passed"] for p in pulses),
              "neutral_surface_restored": bool(neutral.any() and np.max(abs(elevator[neutral]-baseline)) < .03),
              "visible_pitch_excursion": excursion >= 2.,
              "attitude_envelope": bool(np.max(abs(pitch)) < 15),
              "common_alpha_coverage": report["alpha_range_deg"][0] >= -5 and report["alpha_range_deg"][1] <= 0,
              "surface_travel": bool(np.max(abs(elevator)) <= surface_limit_deg),
              "no_saturation": report["max_saturation_s"] == 0,
              "returned_near_trim": final_error is not None and final_error < .5}
    report.update(passed=bool(all(checks.values())), test_kind="manual_elevator_demo",
                  demonstration_checks=checks, elevator_pulses=pulses,
                  elevator_range_deg=[float(min(elevator)), float(max(elevator))],
                  assumed_elevator_limit_deg=surface_limit_deg,
                  pitch_range_deg=[float(min(pitch)), float(max(pitch))],
                  pitch_excursion_deg=excursion, final_5s_pitch_error_deg=final_error,
                  trim_attitude_deg=report.pop("target_attitude_deg"),
                  attitude_excursion_rms_deg=report.pop("attitude_rms_error_deg"))
    return report
