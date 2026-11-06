"""Input-response metrics and figures from SITL doublet runs (hover and forward flight).

Usage (from src/tools/flight_dynamics_model, SITL environment active):

    python scripts/analyse_input_responses.py RUNS_DIR OUT_DIR

RUNS_DIR holds run folders named hover_{pitch,roll,yaw,throttle}_doublet (spearhead.sitl.quad
--doublet ...) and fw_{pitch,roll,yaw,throttle}_doublet (spearhead.sitl.flight <name>_doublet
from a forward trim start); missing runs are skipped. Writes input_response_metrics.json and two
figures. Firmware desired attitude and rates come from each run's DataFlash log.
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymavlink import DFReader

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
RUNS = Path(sys.argv[1])
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
D = np.rad2deg


def truth(run):
    rows = [json.loads(l) for l in (run / "truth.jsonl").open() if l.strip()]
    rows = [r for r in rows if not r["hold"]]
    rep = json.loads((run / "report.json").read_text())
    t0 = rep["release_at_s"]
    t = np.array([r["time_s"] for r in rows]) - t0
    x = np.array([r["state"] for r in rows])
    c = rows
    out = {"t": t, "x": x, "rep": rep, "t0": t0,
           "de": D(np.array([r["control"]["de"] for r in c])),
           "da": D(np.array([r["control"]["da"] for r in c])),
           "dr": D(np.array([r["control"]["dr"] for r in c])),
           "thr": np.array([r["control"]["throttle"] for r in c]),
           "lift": np.array([r["pwm"][4:8] for r in c], dtype=float),
           "tas": np.array([r.get("true_airspeed_m_s") or np.linalg.norm(r["state"][3:6]) for r in c]),
           "alpha": np.array([r.get("alpha_deg") or 0. for r in c]),
           "beta": np.array([r.get("beta_deg") or 0. for r in c]),
           "nz": np.array([-r["sensor"]["imu"]["accel_body"][2]/9.80665 for r in c])}
    return out


def dataflash(run, t0):
    f = sorted(glob.glob(str(run / "logs" / "*.BIN")))[-1]
    log = DFReader.DFReader_binary(f, zero_time_base=False)
    s = {"ATT": [], "RATE": []}
    while True:
        m = log.recv_match(type=list(s))
        if m is None:
            break
        s[m.get_type()].append(m.to_dict())
    o = {}
    for k, rows in s.items():
        tt = np.array([r["TimeUS"] for r in rows])/1e6 - t0
        o[k] = {"t": tt, **{n: np.array([r[n] for r in rows]) for n in rows[0] if n not in ("mavpackettype", "TimeUS")}}
    return o


def settle(t, y, after, band):
    """First time after `after` from which |y| stays within band."""
    m = t >= after
    tt, yy = t[m], y[m]
    out = np.where(np.abs(yy) > band)[0]
    if not out.size:
        return 0.
    if out[-1] == len(yy) - 1:
        return None
    return float(tt[out[-1] + 1] - after)


def rise(t, y, start, target, frac=.9):
    m = (t >= start)
    idx = np.where(m & (np.sign(target)*y >= frac*abs(target)))[0]
    return float(t[idx[0]] - start) if idx.size else None


def hover_attitude(name, axis, dfkey):
    run = RUNS / name
    T = truth(run)
    F = dataflash(run, T["t0"])
    att = F["ATT"]
    base_mask = (att["t"] > 7) & (att["t"] < 9.5)
    a = att[dfkey] - np.median(att[dfkey][base_mask])
    des = att["Des" + dfkey] - np.median(att["Des" + dfkey][base_mask])
    ta = att["t"]
    amp = 8.
    pos = (ta >= 10) & (ta < 12); neg = (ta >= 12) & (ta < 14)
    peak_pos = float(a[pos].max()); peak_neg = float(a[neg].min())
    after = ta >= 14
    # undershoot past zero after neutral is restored
    post = a[after & (ta < 20)]
    other = {"roll": "Pitch", "pitch": "Roll"}[axis]
    o = att[other] - np.median(att[other][base_mask])
    yaw = np.unwrap(np.deg2rad(att["Yaw"])); yaw = D(yaw - np.median(yaw[base_mask]))
    win = (ta >= 9.5) & (ta < 20)
    rate = {"roll": 9, "pitch": 10}[axis]
    m = (T["t"] >= 9.5) & (T["t"] < 20)
    return {
        "run": name, "command": "±8° attitude, 2 s each (QHOVER)",
        "rise_90_s": rise(ta, a, 10, amp),
        "peak_pos_deg": peak_pos, "overshoot_pos_pct": 100*(peak_pos-amp)/amp,
        "peak_neg_deg": peak_neg, "overshoot_neg_pct": 100*(-peak_neg-amp)/amp,
        "rise_90_neg_s": rise(ta, a, 12, -amp),
        "overshoot_past_zero_after_neutral_deg": float(post.max()),
        "settle_1deg_after_neutral_s": settle(ta, a, 14, 1.),
        "max_tracking_error_deg": float(np.max(np.abs(a - des)[win])),
        "peak_rate_deg_s": float(np.max(np.abs(D(T["x"][m, rate])))),
        "cross_axis_peak_deg": float(np.max(np.abs(o[win]))),
        "yaw_coupling_peak_deg": float(np.max(np.abs(yaw[win]))),
        "lift_pwm_range_us": [int(T["lift"][m].min()), int(T["lift"][m].max())],
        "altitude_dev_m": float(np.max(np.abs(-T["x"][m, 2] + T["x"][0, 2]))),
        "drift_45s_m": T["rep"]["horizontal_displacement_m"],
        "passed": T["rep"]["passed"],
    }, (T, F, ta, a, des)


def hover_yaw():
    name = "hover_yaw_doublet"
    run = RUNS / name
    T = truth(run); F = dataflash(run, T["t0"])
    r = F["RATE"]; tr = r["t"]
    y, yd = r["Y"], r["YDes"]
    pos = (tr >= 10.5) & (tr < 12); neg = (tr >= 12.5) & (tr < 14)
    plateau_pos, plateau_neg = float(np.median(y[pos])), float(np.median(y[neg]))
    des_pos, des_neg = float(np.median(yd[pos])), float(np.median(yd[neg]))
    m = (T["t"] >= 9.5) & (T["t"] < 20)
    att = F["ATT"]; base = (att["t"] > 7) & (att["t"] < 9.5)
    win = (att["t"] >= 9.5) & (att["t"] < 20)
    yaw = D(np.unwrap(np.deg2rad(att["Yaw"]))); yaw -= np.median(yaw[base])
    w = (tr >= 9.5) & (tr < 16)
    lags = np.arange(0, 3.0, .01)
    dt = np.median(np.diff(tr))
    corr = [np.corrcoef(yd[w], np.interp(tr[w]+L, tr, y))[0, 1] for L in lags]
    lag = float(lags[int(np.argmax(corr))])
    # Lightly damped heading oscillation after the stick returns to neutral.
    m = (tr >= 14) & (tr < 40)
    tt, yy = tr[m], y[m]
    zc = np.where(np.diff(np.sign(yy)) != 0)[0]
    zc = zc[np.r_[True, np.diff(tt[zc]) > .5]] if zc.size else zc
    peaks = [float(np.max(np.abs(yy[zc[i]:zc[i+1]+1]))) for i in range(len(zc)-1)]
    period = float(2*np.mean(np.diff(tt[zc]))) if len(zc) >= 3 else None
    ratios = [peaks[i+1]/peaks[i] for i in range(len(peaks)-1) if peaks[i] > 0]
    dec = -np.log(np.mean(ratios))*2 if ratios else None
    zeta = float(dec/np.sqrt(4*np.pi**2+dec**2)) if dec else None
    lift_window = (T["t"] >= 9.5) & (T["t"] < 20)
    hs = settle(att["t"], yaw - yaw[-1], 14, 2.)
    return {
        "run": name, "command": "yaw stick ±20 %, 2 s each (QHOVER)",
        "desired_rate_peak_pos_deg_s": float(yd[w].max()), "achieved_rate_peak_pos_deg_s": float(y[w].max()),
        "desired_rate_peak_neg_deg_s": float(yd[w].min()), "achieved_rate_peak_neg_deg_s": float(y[(tr >= 12) & (tr < 20)].min()),
        "lag_s": lag, "yaw_out_peak": float(np.max(np.abs(r["YOut"][w]))),
        "post_input_period_s": period, "post_input_half_cycle_peaks_deg_s": peaks[:5], "post_input_damping_ratio": zeta,
        "heading_settle_2deg_s": hs,
        "heading_change_peak_deg": float(np.max(np.abs(yaw[win]))),
        "roll_pitch_coupling_peak_deg": float(np.max(np.abs(np.r_[att["Roll"][win]-np.median(att["Roll"][base]),
                                                                   att["Pitch"][win]-np.median(att["Pitch"][base])]))),
        "lift_pwm_range_us": [int(T["lift"][lift_window].min()), int(T["lift"][lift_window].max())],
        "passed": T["rep"]["passed"],
        "returned_near_trim": T["rep"]["returned_near_trim"],
        "tail_max_rate_deg_s": float(np.max(np.abs(D(T["x"][T["t"] >= 40, 9:12])))),
    }, (T, F)


def hover_throttle():
    name = "hover_throttle_doublet"
    run = RUNS / name
    T = truth(run)
    t, x = T["t"], T["x"]
    h = -x[:, 2]; h0 = np.median(h[(t > 7) & (t < 9.5)])
    climb = np.gradient(h, t)
    win = (t >= 9.5) & (t < 25)
    pos = (t >= 10) & (t < 12); neg = (t >= 12) & (t < 14)
    tilt = D(np.max(np.abs(x[win, 6:8])))
    return {
        "run": name, "command": "collective ±10 points, 2 s each (QSTABILIZE), QHOVER recovery",
        "vertical_accel_pos_pulse_m_s2": float(np.polyfit(t[(t>=10.3)&(t<11.8)], climb[(t>=10.3)&(t<11.8)], 1)[0]),
        "vertical_accel_neg_pulse_m_s2": float(np.polyfit(t[(t>=12.3)&(t<13.8)], climb[(t>=12.3)&(t<13.8)], 1)[0]),
        "climb_rate_peak_pos_m_s": float(climb[pos].max()), "climb_rate_peak_neg_m_s": float(climb[neg].min()),
        "altitude_peak_m": float((h[win]-h0).max()), "altitude_min_m": float((h[win]-h0).min()),
        "altitude_offset_at_20s_m": float(h[np.argmin(abs(t-20))]-h0),
        "altitude_offset_end_m": float(h[-1]-h0),
        "lift_pwm_range_us": [int(T["lift"][win].min()), int(T["lift"][win].max())],
        "attitude_coupling_peak_deg": float(tilt),
        "passed": T["rep"]["passed"],
    }, T


def fw(name, surface, axis):
    run = RUNS / name
    T = truth(run)
    t, x = T["t"], T["x"]
    rep = T["rep"]; plan = rep["doublet_plan"]
    base = (t > 7) & (t < 9.5)
    end_input = plan["start_s"] + 2*plan["pulse_s"]
    rec = plan["recovery_at_s"]
    win = (t >= 9.5) & (t < 40)
    def rel(y):
        return y - np.median(y[base])
    p, q, r = (D(x[:, 9]), D(x[:, 10]), D(x[:, 11]))
    phi, theta = rel(D(x[:, 6])), rel(D(x[:, 7]))
    psi = rel(D(np.unwrap(x[:, 8])))
    h = rel(-x[:, 2]); V = rel(T["tas"]); al = rel(T["alpha"]); be = rel(T["beta"]); nz = rel(T["nz"])
    o = {"run": name, "passed": rep["passed"], "returned_near_trim": rep["returned_near_trim"],
         "input": f"{plan['channel']} ±{plan['amplitude']:g}{'°' if plan['unit']=='deg' else ''} "
                  f"{plan['pulse_s']:g} s each (MANUAL), FBWA from {rec-plan['start_s']:g} s after input start",
         "trim": plan["trim_value"], "pulse_s": plan["pulse_s"], "fbwa_after_s": rec - plan["start_s"],
         "outside_control_coverage_fraction": rep.get("outside_control_coverage_fraction")}
    pk = lambda y: float(np.max(np.abs(y[win])))
    o.update(peak_p=pk(p), peak_q=pk(q), peak_r=pk(r), peak_dphi=pk(phi), peak_dtheta=pk(theta),
             peak_dpsi=pk(psi), peak_dh=pk(h), peak_dV=pk(V), peak_dalpha=pk(al), peak_dbeta=pk(be), peak_dnz=pk(nz),
             settle_rates_1dps_after_input_s=settle(t, np.max(np.abs(np.c_[p, q, r]), axis=1), end_input, 1.),
             settle_att_2deg_after_input_s=settle(t, np.max(np.abs(np.c_[phi, theta]), axis=1), end_input, 2.),
             altitude_offset_end_m=float(h[-1]), speed_range=rep["speed_range_m_s"])
    if surface in ("dr", "da"):
        # Lateral (Dutch roll) oscillation in r after the input, FBWA active from recovery.
        m = (t >= end_input) & (t < end_input + 12)
        tt, yy = t[m], r[m]
        zc = np.where(np.diff(np.sign(yy)) != 0)[0]
        zc = zc[np.r_[True, np.diff(tt[zc]) > .3]] if zc.size else zc
        if len(zc) >= 3:
            peaks = [float(np.max(np.abs(yy[zc[i]:zc[i+1]+1]))) for i in range(len(zc)-1)]
            o["lateral_osc_period_s"] = float(2*np.mean(np.diff(tt[zc[:5]])))
            ratios = [peaks[i+1]/peaks[i] for i in range(min(4, len(peaks)-1))]
            dec = -np.log(np.mean(ratios))*2
            o["lateral_osc_damping_est"] = float(dec/np.sqrt(4*np.pi**2 + dec**2))
            o["lateral_half_cycle_peaks_deg_s"] = peaks[:5]
    if surface == "de":
        # Per-degree gains during the first pulse.
        m = (t >= plan["start_s"]) & (t < plan["start_s"] + plan["pulse_s"])
        o["q_per_deg_elev"] = float(np.max(np.abs(q[m])))/plan["amplitude"]
        o["nz_per_deg_elev"] = float(np.max(np.abs(nz[m])))/plan["amplitude"]
    if surface == "da":
        m = (t >= plan["start_s"]) & (t < plan["start_s"] + plan["pulse_s"])
        o["p_per_deg_ail"] = float(np.max(np.abs(p[m])))/plan["amplitude"]
        o["yaw_rate_sign_vs_roll_first_pulse"] = "adverse" if np.sign(np.mean(r[m])) != np.sign(np.mean(p[m])) else "proverse"
    series = dict(t=t, p=p, q=q, r=r, phi=phi, theta=theta, psi=psi, h=h, V=V, al=al, be=be, nz=nz,
                  de=T["de"], da=T["da"], dr=T["dr"], thr=T["thr"])
    return o, series


def main():
    results = {"hover": {}, "forward_22": {}}
    hv = {}
    for name, axis, key in (("hover_pitch_doublet", "pitch", "Pitch"), ("hover_roll_doublet", "roll", "Roll")):
        if (RUNS/name/"report.json").exists():
            results["hover"][axis], hv[axis] = hover_attitude(name, axis, key)
    if (RUNS/"hover_yaw_doublet"/"report.json").exists():
        results["hover"]["yaw"], hv["yaw"] = hover_yaw()
    if (RUNS/"hover_throttle_doublet"/"report.json").exists():
        results["hover"]["throttle"], hv["throttle"] = hover_throttle()
    fs = {}
    for name, s, ax in (("fw_pitch_doublet", "de", 1), ("fw_roll_doublet", "da", 0), ("fw_yaw_doublet", "dr", 2),
                        ("fw_throttle_doublet", "throttle", None)):
        if (RUNS/name/"report.json").exists():
            results["forward_22"][name.split("_")[1]], fs[name.split("_")[1]] = fw(name, s, ax)
    (OUT/"input_response_metrics.json").write_text(json.dumps(results, indent=1))

    # Hover figure
    fig, axs = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax_, axis in ((axs[0, 0], "pitch"), (axs[0, 1], "roll")):
        if axis in hv:
            T, F, ta, a, des = hv[axis]
            ax_.plot(ta, des, "k--", lw=1, label="firmware desired")
            ax_.plot(ta, a, "C0", label="attitude")
            ax_.set_title(f"Hover {axis} doublet (±8° attitude)"); ax_.set_ylabel("deg"); ax_.legend(fontsize=8)
    if "yaw" in hv:
        T, F = hv["yaw"]; r = F["RATE"]
        axs[1, 0].plot(r["t"], r["YDes"], "k--", lw=1, label="desired yaw rate")
        axs[1, 0].plot(r["t"], r["Y"], "C0", label="yaw rate")
        axs[1, 0].set_title("Hover yaw doublet (stick ±20 %)"); axs[1, 0].set_ylabel("deg/s"); axs[1, 0].legend(fontsize=8)
    if "throttle" in hv:
        T = hv["throttle"]; t, x = T["t"], T["x"]
        h = -x[:, 2]; h -= np.median(h[(t > 7) & (t < 9.5)])
        axs[1, 1].plot(t, h, "C0", label="altitude change (m)")
        axs[1, 1].plot(t, np.gradient(-x[:, 2], t), "C1", label="climb rate (m/s)")
        axs[1, 1].set_title("Hover collective doublet (±10 points)"); axs[1, 1].legend(fontsize=8)
    for a_ in axs.flat:
        a_.set_xlim(8, 25); a_.grid(alpha=.3); a_.axvspan(10, 14, color="0.9", zorder=0)
    for a_ in axs[1]:
        a_.set_xlabel("time after release (s)")
    fig.suptitle("PT1.5 hover input responses, ArduPlane 4.7.1 SITL, CG x −0.040 m, 21 kg")
    fig.tight_layout(); fig.savefig(OUT/"hover_input_responses.png", dpi=130)

    fig, axs = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    spec = {"pitch": (axs[0, 0], "de", [("q", "q (deg/s)"), ("theta", "Δθ (deg)"), ("nz", "Δnz (g ×10)")]),
            "roll": (axs[0, 1], "da", [("p", "p (deg/s)"), ("phi", "Δφ (deg)"), ("r", "r (deg/s)")]),
            "yaw": (axs[1, 0], "dr", [("r", "r (deg/s)"), ("be", "Δβ (deg)"), ("phi", "Δφ (deg)")]),
            "throttle": (axs[1, 1], "thr", [("V", "ΔTAS (m/s)"), ("h", "Δh (m)"), ("theta", "Δθ (deg)")])}
    for k, (ax_, inp, lines) in spec.items():
        if k not in fs:
            continue
        s = fs[k]; t = s["t"]; base = (t > 7) & (t < 9.5)
        u = s[inp] - np.median(s[inp][base])
        ax_.plot(t, u*(100 if inp == "thr" else 1), "k", lw=1, label="input Δ" + ("(%)" if inp == "thr" else "(deg)"))
        for key, lab in lines:
            ax_.plot(t, s[key]*(10 if key == "nz" else 1), label=lab)
        res = results["forward_22"][k]
        ax_.axvspan(10, 10 + 2*res["pulse_s"], color="0.9", zorder=0)
        ax_.axvline(10 + res["fbwa_after_s"], color="0.5", lw=.8, ls=":")
        ax_.set_title(f"22 m/s {k} doublet"); ax_.legend(fontsize=7, ncol=2); ax_.grid(alpha=.3)
    for a_ in axs.flat:
        a_.set_xlim(8, 35)
    for a_ in axs[1]:
        a_.set_xlabel("time after release (s)")
    fig.suptitle("PT1.5 forward-flight input responses, 22 m/s, MANUAL input then FBWA, CG x −0.040 m")
    fig.tight_layout(); fig.savefig(OUT/"forward_input_responses.png", dpi=130)
    print(json.dumps(results, indent=1))


main()
