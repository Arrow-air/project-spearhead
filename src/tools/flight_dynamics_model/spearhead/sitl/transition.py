"""Prototype-only hover-to-FBWA diagnostic and truth-based acceptance checks."""

import json
from pathlib import Path

import numpy as np


def transition_parameters(parameter_set="pt15"):
    # Use firmware's simulated pitot, calibrated during the stationary bootstrap.
    # Keep the current hover tune; forward surface trim is calculated separately.
    if parameter_set == "first_transition_2026_10_09":
        # Frozen earlier proposal, including the verified factory rate gains.
        # These are a comparison baseline, not a new tune or hardware upload.
        return {"ARSPD_TYPE": 100, "ARSPD_USE": 1,
                "AIRSPEED_MIN": 18, "Q_ASSIST_SPEED": 15,
                "AIRSPEED_CRUISE": 25, "AIRSPEED_MAX": 27, "SCALING_SPEED": 25,
                "Q_TRANSITION_MS": 5000, "Q_TRANS_FAIL": 20, "Q_TRANS_FAIL_ACT": 1,
                "Q_OPTIONS": 262177, "ROLL_LIMIT_DEG": 25,
                "PTCH_LIM_MAX_DEG": 10, "PTCH_LIM_MIN_DEG": -10,
                "PTCH2SRV_TCONST": .75, "RLL2SRV_RMAX": 45,
                "PTCH2SRV_RMAX_UP": 22, "PTCH2SRV_RMAX_DN": 22,
                "FLTMODE3": 19, "FLTMODE4": 5, "FLTMODE5": 21, "FLTMODE6": 20,
                "FS_LONG_ACTN": 1, "Q_RTL_MODE": 1,
                "BATT_FS_LOW_ACT": 6, "BATT_FS_CRT_ACT": 4, "BATT_ARM_VOLT": 46.8,
                "RLL_RATE_FF": .345, "RLL_RATE_P": .08, "RLL_RATE_I": .15, "RLL_RATE_D": 0,
                "PTCH_RATE_FF": .345, "PTCH_RATE_P": .04, "PTCH_RATE_I": .15, "PTCH_RATE_D": 0}
    if parameter_set != "pt15":
        raise ValueError("Unknown transition parameter set")
    return {"ARSPD_TYPE": 100, "ARSPD_USE": 1,
            "AIRSPEED_MIN": 20, "AIRSPEED_CRUISE": 24, "AIRSPEED_MAX": 30,
            "SCALING_SPEED": 24, "Q_TRANSITION_MS": 5000,
            "ROLL_LIMIT_DEG": 30, "PTCH_LIM_MIN_DEG": -10,
            "PTCH_LIM_MAX_DEG": 10}


def forward_trim_parameters(params, mapping, altitude, *, include_simulation_gains=True):
    from .physics import full_trim
    from .flight import flight_parameters
    state, control = full_trim(params, 24, altitude)
    proposed = flight_parameters(mapping, state, control)
    # Preserve hover pilot inputs, motor setup and arming configuration.
    keys = {"SERVO1_TRIM", "SERVO2_TRIM", "SERVO4_TRIM", "PTCH_TRIM_DEG",
            "STAB_PITCH_DOWN", "KFF_THR2PTCH", "TRIM_THROTTLE",
            "PTCH_RATE_FF", "PTCH_RATE_P", "PTCH_RATE_I", "PTCH2SRV_TCONST",
            "RLL_RATE_FF", "RLL_RATE_P", "RLL_RATE_I"}
    # AP_Int8 stores a whole percent. Verify exactly what we send/read back.
    proposed["TRIM_THROTTLE"] = round(proposed["TRIM_THROTTLE"])
    if not include_simulation_gains:
        keys -= {key for key in keys if key.startswith(("PTCH_RATE_", "RLL_RATE_")) or key == "PTCH2SRV_TCONST"}
    return {key: proposed[key] for key in keys}, state, control


def transition_command(elapsed, confirmed_mode, forward_confirmed_at=None):
    """Pilot command only; actual surfaces/motors remain firmware-controlled."""
    if elapsed < 10:
        return 18, 1500
    # Never send forward throttle while still in a height-controlled Q mode.
    if confirmed_mode != 5 or forward_confirmed_at is None:
        return 5, 1500
    throttle = round(1500 + 400 * min(1, max(0, (elapsed - forward_confirmed_at) / 5)))
    return 5, throttle


def confirmed_mode_samples(times, telemetry_path, mode, max_age_s=2.):
    """As-of join: only past, fresh FBWA heartbeats confirm a physical sample."""
    heartbeats = []
    if Path(telemetry_path).exists():
        with Path(telemetry_path).open() as stream:
            for line in stream:
                row = json.loads(line)
                message = row.get("message", {})
                if message.get("mavpackettype") == "HEARTBEAT":
                    stamp = row.get("fdm_time_s")
                    if stamp is not None and np.isfinite(stamp):
                        heartbeats.append((stamp, message.get("custom_mode")))
    confirmed = np.zeros(len(times), dtype=bool)
    if heartbeats:
        heartbeats.sort(key=lambda pair: pair[0])
        stamps = np.array([pair[0] for pair in heartbeats])
        modes = np.array([pair[1] == mode for pair in heartbeats])
        indices = np.searchsorted(stamps, times, side="right")-1
        valid = indices >= 0
        previous = np.maximum(indices, 0)
        confirmed = valid & modes[previous] & (times-stamps[previous] <= max_age_s)
    return confirmed, bool(heartbeats)


def confirmed_forward_samples(times, telemetry_path, max_age_s=2.):
    return confirmed_mode_samples(times, telemetry_path, 5, max_age_s)


def analyse_transition(path, release_at, duration, initial_altitude, final_mode, idle_lift_pwm=1185,
                       *, telemetry_path=None):
    rows = []
    with Path(path).open() as stream:
        for line in stream:
            row = json.loads(line)
            if not row["hold"]:
                rows.append(row)
    if not rows or release_at is None:
        return {"passed": False, "transition_completed": False,
                "reason": "No released physics frames"}
    times = np.array([r["time_s"] for r in rows])
    t = times - release_at
    x = np.array([r["state"] for r in rows])
    pwm = np.array([r["pwm"] for r in rows])
    speed = np.linalg.norm(x[:, 3:6], axis=1)
    altitude = -x[:, 2]
    tilt = np.rad2deg(x[:, 6:8])
    confirmed, evidence = confirmed_forward_samples(
        times, telemetry_path or Path(path).with_name("telemetry.jsonl"))
    forward = ((t >= 10) & confirmed & (speed >= 20) & (speed <= 30)
               & (np.max(pwm[:, 4:8], axis=1) <= idle_lift_pwm)
               & (np.max(abs(tilt), axis=1) <= 15)
               & (abs(altitude - initial_altitude) <= 15))
    longest = current = 0.
    for i in range(1, len(t)):
        current = current + t[i] - t[i-1] if forward[i] and forward[i-1] else 0.
        longest = max(longest, current)
    checks = {"free_duration": bool(t[-1] >= duration - .02),
              "confirmed_fbwa": final_mode == 5,
              "sustained_forward_without_lift": bool(longest >= 10)}
    return {"passed": all(checks.values()), "checks": checks,
            "transition_completed": all(checks.values()),
            "free_duration_s": float(t[-1]),
            "longest_forward_without_lift_s": float(longest),
            "max_true_airspeed_m_s": float(speed.max()),
            "final_true_airspeed_m_s": float(speed[-1]),
            "altitude_range_m": [float(altitude.min()), float(altitude.max())],
            "pitch_range_deg": [float(tilt[:, 1].min()), float(tilt[:, 1].max())],
            "peak_abs_pitch_deg": float(abs(tilt[:, 1]).max()),
            "altitude_loss_m": float(max(0, initial_altitude-altitude.min())),
            "hover_peak_abs_pitch_deg": float(abs(tilt[(t >= 0) & (t < 10), 1]).max())
                                       if np.any((t >= 0) & (t < 10)) else None,
            "mode_evidence": {"source": "recorded HEARTBEAT", "available": evidence,
                              "maximum_age_s": 2., "required_during_forward_interval": True},
            "final_lift_pwm_us": pwm[-1, 4:8].tolist(),
            "final_pusher_pwm_us": int(pwm[-1, 2]),
            "idle_lift_pwm_us": idle_lift_pwm,
            "acceptance": f"10 continuous seconds in heartbeat-confirmed FBWA (maximum age 2 s) at 20–30 m/s, lift PWM <= armed-idle {idle_lift_pwm}, roll/pitch <=15 degrees, altitude within 15 m; final FBWA and full requested duration"}
