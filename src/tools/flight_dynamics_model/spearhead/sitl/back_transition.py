"""FW trim to position-holding hover; firmware owns all motor handover."""

import json
from pathlib import Path

import numpy as np

from .transition import confirmed_mode_samples

HOVER_MODE = 19  # QLOITER actively decelerates and holds the stopping position.


def forward_release_ready(pwm, mapping, params, trim_throttle, transition_done):
    idle = np.ceil(params["Q_M_PWM_MIN"]+params["Q_M_SPIN_ARM"]*(params["Q_M_PWM_MAX"]-params["Q_M_PWM_MIN"]))+1
    return bool(transition_done and len(pwm) >= 8 and np.max(pwm[4:8]) <= idle and
                abs(mapping.pusher.motor(pwm)-trim_throttle) <= .02)


def back_transition_command(elapsed, confirmed_mode, trim_throttle, dead_zone=60):
    desired = 5 if elapsed < 10 else HOVER_MODE
    # Retain forward throttle until the requested hover mode is confirmed.
    rc3 = 1500 if elapsed >= 10 and confirmed_mode == HOVER_MODE else round(1000+dead_zone+(1000-dead_zone)*trim_throttle)
    return desired, rc3


def analyse_back_transition(path, release_at, duration, initial_altitude, final_mode, idle_pwm=1185):
    rows = [json.loads(line) for line in Path(path).open() if line.strip()]
    rows = [r for r in rows if not r["hold"]]
    if not rows or release_at is None:
        return {"passed": False, "test_kind": "back_transition", "reason": "No free-flight frames"}
    times = np.array([r["time_s"] for r in rows]); t = times-release_at
    state = np.array([r["state"] for r in rows]); pwm = np.array([r["pwm"] for r in rows])
    speed = np.linalg.norm(state[:, 3:6], axis=1); tilt = np.rad2deg(state[:, 6:8]); altitude = -state[:, 2]
    # The inertial vertical rate comes from the recorded body state and attitude.
    phi, theta = state[:, 6], state[:, 7]
    vertical = -np.sin(theta)*state[:, 3]+np.sin(phi)*np.cos(theta)*state[:, 4]+np.cos(phi)*np.cos(theta)*state[:, 5]
    confirmed, evidence = confirmed_mode_samples(times, Path(path).with_name("telemetry.jsonl"), HOVER_MODE)
    hover = ((t >= 10) & confirmed & (speed <= 2) & (abs(vertical) <= 1) & (np.max(abs(tilt), axis=1) <= 10)
             & (abs(altitude-initial_altitude) <= 15) & (np.min(pwm[:, 4:8], axis=1) > idle_pwm)
             & np.array([r["control"]["throttle"] <= .002 for r in rows]))
    longest = current = 0.
    for i in range(1, len(t)):
        current = current+times[i]-times[i-1] if hover[i] and hover[i-1] else 0.
        longest = max(longest, current)
    checks = {"free_duration": bool(t[-1] >= duration-.02), "confirmed_qloiter": bool(final_mode == HOVER_MODE and confirmed[-1]),
              "sustained_hover": bool(longest >= 10), "final_hover_state": bool(hover[-1])}
    return {"passed": all(checks.values()), "test_kind": "back_transition", "checks": checks,
            "free_duration_s": float(t[-1]), "longest_hover_s": float(longest), "final_true_airspeed_m_s": float(speed[-1]),
            "peak_abs_pitch_deg": float(np.max(abs(tilt[:, 1]))), "altitude_loss_m": float(max(0, initial_altitude-altitude.min())),
            "initial_forward_peak_pitch_deg": float(np.max(abs(tilt[t < 10, 1]))) if np.any(t < 10) else None,
            "final_lift_pwm_us": pwm[-1, 4:8].tolist(), "final_pusher_fraction": float(rows[-1]["control"]["throttle"]),
            "mode_evidence": {"available": evidence, "mode": "QLOITER", "maximum_age_s": 2., "required_during_hover_interval": True}}
