"""Web Flight Deck for the external ArduPlane runner, without a second FDM."""

import asyncio
from collections import deque
import math
from pathlib import Path
import time

from .flight_deck import attitude_html
from .plotting import _flight_deck_base_figure
from .flight_setup import (AIRCRAFT, POINTS, MANEUVERS, MIN_DURATION, FlightSetup,
                           TRANSITION_CGS, TRANSITION_CG_X, TRANSITION_CG_Z, TRANSITION_PARAMETERS,
                           FAILURES, GUST_START_S,
                           forward_start, maneuver_options)
from ..aircraft import aircraft_summary, load_aircraft_config
from ..sitl.ground_control import QUADPLANE_MODES, PLANE_MODES, MODE_HELP, MODE_NAMES, session_modes
from ..sitl.doublets import DOUBLETS, doublet_delta
from ..sitl.hover_doublets import HOVER_DOUBLETS
from .flight_setup import HOVER_DOUBLET_LABELS
from .flight_values import VALUE_LABELS, numeric_flight_values, flight_values_status, live_alerts


def info_note(ui, label, text=""):
    """Hover, keyboard focus or click reveals explanatory text."""
    button = ui.button(icon="info_outline").props(f'flat round dense aria-label="{label}"').classes("fd-info")
    with button:
        note = ui.tooltip(text).classes("fd-info-tooltip")
    button.on("focus", lambda: note.run_method("show"))
    button.on("blur", lambda: note.run_method("hide"))
    button.on_click(lambda: note.run_method("show"))
    return note


def arming_badge(view):
    armed = view.get("armed")
    label = "ARMING UNKNOWN" if armed is None else "ARMED" if armed else "DISARMED"
    tone = "unknown" if armed is None else "armed" if armed else "disarmed"
    if view.get("active") and view.get("stale_mavlink") and armed is not None:
        label += " · STALE"
        tone = "unknown"
    elif not view.get("active") and view.get("run_id"):
        label = "RECORDED · " + label
    return label, tone


def startup_activity(view):
    """Describe observed startup state, without driving or bypassing the runner."""
    if not view.get("active") or view.get("free_time_s") is not None or view.get("error"):
        return None
    row = view.get("truth")
    elapsed = view.get("wall_elapsed_s") or 0
    if row is None:
        title, detail = (("Starting flight engine", "Waiting for the first simulation frame.") if view.get("manifest") else
                         ("Preparing aircraft model", "Checking the setup and calculating the initial trim."))
    elif view.get("stale_truth") or view.get("stale_mavlink"):
        title, detail = "Waiting for simulator updates", "The process is running; waiting for fresh physics and firmware telemetry."
    elif view.get("armed"):
        hover = (view.get("manifest") or {}).get("config", {}).get("flight_regime") == "hover"
        title, detail = (("Spooling lift motors", "Armed. Waiting for hover support before releasing the aircraft.") if hover else
                         ("Preparing flight release", "Armed. Completing initialization before flight begins."))
    else:
        health = view.get("health") or {}
        if (health.get("gps_fix") or 0) < 3:
            title, detail = "Acquiring GPS", "Waiting for a simulated 3D GPS fix."
        elif not health.get("ekf_ok"):
            title, detail = "Aligning navigation", "GPS acquired. The attitude and navigation estimate is initializing."
        else:
            title, detail = "Pre-arm checks", "Navigation ready. Firmware checks and normal arming are in progress."
    # Prior local runs released at about 47 simulated seconds (~1–2 min wall time).
    # Refine the estimate from the measured clock, keeping it explicitly approximate.
    sim_time = None if row is None else row["time_s"]
    remaining = None
    if sim_time is not None and sim_time >= 5 and elapsed > 0:
        remaining = max(0, 50-sim_time) * elapsed/sim_time
    estimate = (f"About {math.ceil(remaining/10)*10:g} s to flight · estimate varies" if remaining is not None and remaining >= 10 else
                "Flight should begin shortly · estimate varies" if remaining is not None and remaining > 0 else
                "Usually 1–2 min before flight · may take longer" if remaining is None else
                "Taking longer than estimated · checks are still running")
    return {"title": title, "detail": detail, "elapsed_s": elapsed, "estimate": estimate}


def flight_story(view):
    """Plain-language context from recorded state, never a second test verdict."""
    result = view.get("report") or {}
    if view.get("error"):
        return "Flight needs attention", "The simulator reported an error. Read the message below before starting again.", "error", 2
    if not view.get("active"):
        if view.get("status") == "STOPPED":
            return "Flight stopped", "This recording ends where the run was stopped. It does not confirm a completed test.", "neutral", 2
        if result:
            if result.get("test_kind") == "control_doublet" and result.get("passed"):
                return "Doublet recorded", "Both physical input pulses were delivered. Review the response and the separate recovery result below.", "success", 2
            if result.get("test_kind") == "interactive_session" and result.get("passed"):
                return "Flight complete", "Your session has ended. Review the aircraft motion, flight plots and messages below.", "neutral", 2
            if result.get("passed"):
                if result.get("test_kind") == "quad_transition":
                    detail = (f"{result.get('free_duration_s', 0):.0f} s completed · "
                              f"{result.get('longest_forward_without_lift_s', 0):.1f} s met the forward-flight checks without lift support.")
                elif result.get("test_kind") == "back_transition":
                    detail = f"{result.get('free_duration_s', 0):.0f} s completed · {result.get('longest_hover_s', 0):.1f} s met the hover checks."
                else:
                    detail = "The recorded run met its simulation checks. Review the plots and detailed results below."
                return "Simulation checks passed", detail, "success", 2
            return "Simulation checks failed", "The aircraft did not meet all test criteria. Review the result and messages below.", "error", 2
        if view.get("run_id"):
            return "Recording available", "This run has no completed test report. The displays show its last recorded state.", "neutral", 2
        return "Ready for a flight", "Choose a starting position and flight test, then select Start flight.", "neutral", 0
    if view.get("free_time_s") is None:
        return "Preparing your flight", "The aircraft is held in place while initialization, arming and motor startup finish.", "neutral", 1
    if view.get("stale_truth") or view.get("stale_mavlink"):
        return "Waiting for fresh flight data", "The displays may be out of date. Check the connection warnings below.", "warning", 1
    manifest = view.get("manifest") or {}
    if manifest.get("test_kind") == "back_transition":
        if view.get("mode") == 19:
            return "Braking into hover · QLOITER", "The pusher is stopped by the firmware. Lift motors and aerodynamic controls brake the aircraft and hold position; recorded checks determine whether hover settles.", "neutral", 1
        return "Holding the forward trim", "The aircraft settles in forward flight for 10 s, then requests position-holding hover.", "neutral", 1
    if manifest.get("doublet_plan"):
        plan = manifest["doublet_plan"]
        elapsed = view.get("free_time_s") or 0
        stage = "Settling at trim" if elapsed < 9 else "Applying the doublet" if elapsed < plan["recovery_at_s"] else "Observing recovery · " + ("QHOVER" if plan.get("regime") == "hover" else "FBWA")
        return stage, "A positive pulse, a matching negative pulse, then neutral input. Actual inputs and response are recorded; the recovery controller uses frozen gains.", "neutral", 1
    if manifest.get("test_kind") == "quad_transition":
        if view.get("mode") != 5:
            return "Holding hover", "The lift motors support the aircraft. The test will request forward flight after the initial hover.", "neutral", 1
        pwm = (view.get("truth") or {}).get("pwm") or []
        params = manifest.get("parameter_overlay") or {}
        # Use the runner's armed-idle calculation; this identifies motor state,
        # not stable-forward-flight acceptance or a passed transition.
        idle = None
        if all(k in params for k in ("Q_M_PWM_MIN", "Q_M_PWM_MAX", "Q_M_SPIN_ARM")):
            idle = math.ceil(params["Q_M_PWM_MIN"] + params["Q_M_SPIN_ARM"] *
                             (params["Q_M_PWM_MAX"]-params["Q_M_PWM_MIN"]))
        if idle is not None and len(pwm) >= 8 and all(pwm[i] <= idle for i in (4, 5, 6, 7)):
            return "Forward flight · lift motors at idle", "The test is checking speed, attitude and altitude. A completed handover alone does not mean the transition has passed.", "neutral", 1
        return "Transitioning to forward flight", "The pusher is accelerating the aircraft while the firmware hands over from lift motors to wing-supported flight.", "neutral", 1
    if manifest.get("ground_contact"):
        return "Takeoff, hover and landing", "The automatic mission controls the aircraft. Flight time includes the takeoff and landing sequence.", "neutral", 1
    if view.get("interactive"):
        return "Interactive flight is running", "Watch the aircraft below. You can request a flight mode once the controls are ready.", "neutral", 1
    return "Flight test is running", "The selected test controls the aircraft automatically. Watch its motion and review the results when it finishes.", "neutral", 1


def flight_plan_text(setup):
    text = _flight_plan_base(setup)
    extras = []
    if setup.wind_speed_m_s:
        extras.append(f"Wind {setup.wind_speed_m_s:g} m/s from {setup.wind_from_deg:g}°")
    if setup.gust_m_s:
        extras.append(f"{setup.gust_m_s:g} m/s gust at {GUST_START_S:g} s")
    if setup.failure != "none":
        extras.append(f"{FAILURES[setup.failure]} at {setup.failure_time_s:g} s")
    return text + (" " + " · ".join(extras) + "." if extras else "")


def _flight_plan_base(setup):
    if setup.initial_point == "ground":
        return f"Take off to {setup.altitude_m:g} m, hover for {setup.duration_s:g} s, then land automatically."
    if setup.maneuver == "transition":
        return f"Start hovering at {setup.altitude_m:g} m → hold for 10 s → accelerate into forward flight. {setup.duration_s:g} s of flight."
    if setup.maneuver == "back_transition":
        return f"Start at {setup.speed_m_s:g} m/s forward trim → settle for 10 s → brake into position-holding hover. {setup.duration_s:g} s of flight."
    if setup.maneuver in DOUBLETS and not setup.is_forward:
        return f"Start hovering at {setup.altitude_m:g} m → settle for 10 s → {HOVER_DOUBLET_LABELS[setup.maneuver]}, 2 s per pulse → QHOVER recovery. {setup.duration_s:g} s of flight."
    if setup.maneuver in DOUBLETS:
        plan = DOUBLETS[setup.maneuver]
        amplitude = "up to ±10 throttle points" if plan["unit"] == "fraction" else f"±{plan['amplitude']:g}° {plan['channel']}"
        return f"Start at {setup.speed_m_s:g} m/s trim → {amplitude}, {plan['pulse_s']:g} s per pulse → observe FBWA recovery. {setup.duration_s:g} s of flight."
    if setup.aircraft == "pt15" and not setup.is_forward:
        return f"Start hovering at {setup.altitude_m:g} m. " + (
            f"Change flight modes during a {setup.duration_s:g} s session." if setup.maneuver == "interactive" else
            f"Hold hover automatically for {setup.duration_s:g} s.")
    return f"Start flying at {setup.speed_m_s:g} m/s and {setup.altitude_m:g} m. Run the selected test for {setup.duration_s:g} s."


MODE_LABELS = {5: "Forward flight", 17: "Level hover", 18: "Altitude-hold hover",
               19: "Position-hold hover", 20: "Vertical landing", 21: "Return and land",
               10: "Automatic mission", 0: "Manual control", 16: "Initializing"}


class SimulationPace:
    """Observe the simulation clock against wall time; never change its rate."""

    def __init__(self):
        self.samples = deque(maxlen=40)
        self.run_id = None

    def update(self, view, now):
        if view.get("run_id") != self.run_id or not view.get("active"):
            self.samples.clear()
            self.run_id = view.get("run_id")
        if not view.get("active") or view.get("time_s") is None:
            return None
        self.samples.append((now, view["time_s"]))
        while len(self.samples) > 2 and self.samples[1][0] < now-3:
            self.samples.popleft()
        elapsed = now-self.samples[0][0]
        if elapsed < 1:
            return None
        return max(0., (view["time_s"]-self.samples[0][1])/elapsed)


def control_coverage_text(aircraft):
    data = aircraft.get("control_data") or {}
    coverage = data.get("coverage") or {}
    parts = []
    for name in ("elevator", "aileron", "rudder"):
        surface = coverage.get(name)
        if surface:
            lo, hi = surface["deflection_deg"]
            speeds = surface.get("measured_velocities_m_s")
            parts.append(f"{name} {lo:g}…{hi:g}°" +
                         (" at " + "/".join(f"{v:g}" for v in speeds) + " m/s" if speeds else ""))
    return ("Measured bounds: " + " · ".join(parts) +
            f" · β={data.get('beta_deg', 0):g}°; sparse coverage depends on α and speed."
            if parts else "Measured coverage unavailable; consult the recorded manifest.")


def transition_result_text(result):
    if result.get("test_kind") == "back_transition":
        outcome = "PASSED" if result.get("passed") else "FAILED"
        return f"{outcome} · Accepted hover: {result.get('longest_hover_s', 0):.1f} s (requires 10 s) · Final speed {result.get('final_true_airspeed_m_s', 0):.2f} m/s · Pusher {100*result.get('final_pusher_fraction', 0):.1f}%"
    done = "yes" if result.get("firmware_transition_done") else "no"
    outcome = "PASSED" if result.get("passed") else "FAILED"
    return (f"{outcome} · Firmware handover: {done} · "
            f"Forward acceptance without lift: {result.get('longest_forward_without_lift_s', 0):.1f} s (requires 10 s) · "
            f"Free flight: {result.get('free_duration_s', 0):.1f} s · "
            f"Final TAS: {result.get('final_true_airspeed_m_s', 0):.2f} m/s")


def transition_metrics(view):
    """Recorded/report metrics, separate from the sustained-interval verdict."""
    report = view.get("report") or {}
    observed = view.get("transient_metrics") or {}
    values = {key: report.get(key, observed.get(key)) for key in
            ("longest_forward_without_lift_s", "peak_abs_pitch_deg", "altitude_loss_m", "hover_peak_abs_pitch_deg")}
    if (view.get("manifest") or {}).get("test_kind") == "back_transition":
        values["longest_forward_without_lift_s"] = report.get("longest_hover_s")
        values["hover_peak_abs_pitch_deg"] = report.get("initial_forward_peak_pitch_deg", observed.get("initial_forward_peak_pitch_deg"))
    return values


def doublet_result_text(result):
    plan = result.get("doublet_plan", {})
    pulses = result.get("pulses", [])
    scale = 100 if plan.get("unit") == "fraction" else 1
    unit = "points" if scale == 100 else "°"
    inputs = " / ".join("missing" if p.get("measured") is None else f"{scale*p['measured']:+.2f}{unit}" for p in pulses)
    rates = result.get("peak_rate_deg_s", [])
    peaks = " · ".join(f"{axis} {value:.1f}°/s" for axis, value in zip("pqr", rates))
    if plan.get("regime") == "hover":
        return (f"Measured pilot inputs: {inputs} · Peak rates {peaks} · Attitude settled: {'yes' if result.get('returned_near_trim') else 'no'} · "
                f"Final altitude offset {result.get('final_altitude_offset_m', 0):+.2f} m · Horizontal drift {result.get('horizontal_displacement_m', 0):.1f} m. QHOVER does not hold position. Input delivery is assessed separately.")
    return f"Measured input increments: {inputs} · Peak rates {peaks} · Returned near trim: {'yes' if result.get('returned_near_trim') else 'no'} · Input delivery and recovery are assessed separately."


def lift_pwm_limits(manifest):
    """Use the recording's motor endpoints, not SERVOx limits or observed peaks."""
    manifest = manifest or {}
    params = manifest.get("parameter_overlay", {})
    low, high = params.get("Q_M_PWM_MIN"), params.get("Q_M_PWM_MAX")
    if low is None or high is None:
        channels = manifest.get("config", {}).get("channels", {}).get("lift", [])
        if not channels or any("minimum" not in c or "maximum" not in c for c in channels):
            return None
        low = min(c["minimum"] for c in channels)
        high = max(c["maximum"] for c in channels)
    if not all(isinstance(v, (int, float)) and math.isfinite(v) for v in (low, high)) or low >= high:
        return None
    return low, high


def comparison_figures(result):
    """Aircraft log (solid) against the SITL run (dotted), both from takeoff / release."""
    traces = result.get("traces")
    if not traces:
        return {}
    sim, real = traces["simulation"], traces["aircraft"]
    figures = {}
    go, figure = _flight_deck_base_figure("Lift outputs · aircraft vs SITL")
    for name, color in (("FR", "#38bdf8"), ("RL", "#a78bfa"), ("FL", "#fbbf24"), ("RR", "#34d399")):
        key = f"lift_{name}"
        if key in real:
            figure.add_trace(go.Scatter(x=real["time_s"], y=real[key], name=f"{name} aircraft", line={"color": color}))
        figure.add_trace(go.Scatter(x=sim["time_s"], y=sim[key], name=f"{name} SITL", line={"color": color, "dash": "dot"}))
    figure.update_yaxes(title_text="PWM · µs")
    figures["lift"] = figure
    go, figure = _flight_deck_base_figure("Attitude · aircraft vs SITL")
    for key, label, color in (("roll_deg", "Roll", "#38bdf8"), ("pitch_deg", "Pitch", "#fbbf24")):
        if key in real:
            figure.add_trace(go.Scatter(x=real["time_s"], y=real[key], name=f"{label} aircraft", line={"color": color}))
        figure.add_trace(go.Scatter(x=sim["time_s"], y=sim[key], name=f"{label} SITL", line={"color": color, "dash": "dot"}))
    figure.update_yaxes(title_text="deg")
    figures["attitude"] = figure
    go, figure = _flight_deck_base_figure("Altitude and pack voltage · aircraft vs SITL")
    for key, label, color in (("altitude_m", "Altitude · m", "#34d399"), ("battery_v", "Pack · V", "#f87171")):
        if key in real:
            figure.add_trace(go.Scatter(x=real["time_s"], y=real[key], name=f"{label} aircraft", line={"color": color}))
        figure.add_trace(go.Scatter(x=sim["time_s"], y=sim[key], name=f"{label} SITL", line={"color": color, "dash": "dot"}))
    figures["altitude"] = figure
    for figure in figures.values():
        figure.update_xaxes(title_text="s from takeoff / release")
    return figures


def comparison_text(result):
    lines = [f"{result['log']}: {result['aircraft_parameter_count']} aircraft parameters; SITL values from the {result['simulation_parameter_source']}."]
    flights = result.get("flights") or []
    if flights:
        chosen = result.get("takeoff_s")
        listing = ", ".join(f"#{i} at {f['takeoff_s']:.0f} s ({f['airborne_s']:.0f} s)" + (" ← shown" if f["takeoff_s"] == chosen else "")
                            for i, f in enumerate(flights, 1))
        lines.append("Flights in log: " + listing + ".")
    hover = result.get("hover") or {}
    if hover.get("aircraft") and hover.get("simulation"):
        a, s = hover["aircraft"], hover["simulation"]
        lines.append("Median lift PWM after 5 s · aircraft FR/RL/FL/RR " + "/".join(f"{a['median_pwm'][k]:.0f}" for k in ("FR", "RL", "FL", "RR"))
                     + " vs SITL " + "/".join(f"{s['median_pwm'][k]:.0f}" for k in ("FR", "RL", "FL", "RR"))
                     + f" · front share of PWM above idle {100*a['front_share_pwm']:.1f} % vs {100*s['front_share_pwm']:.1f} %.")
    elif result.get("traces") is None and flights == []:
        lines.append("Parameter comparison only (no flight data in this file).")
    lines.append(f"{len(result['parameter_diff'])} flight-relevant parameters differ or exist only in SITL "
                 "(calibration, telemetry, hardware and simulator settings are left out).")
    return " ".join(lines)


def flight_time_origin(snapshot):
    """Simulator time of release (ground starts: free physics start); None before release.

    The run duration counts from this instant, so plots and messages use it as zero and the
    firmware startup before it appears at negative times.
    """
    release = (snapshot.get("report") or {}).get("release_at_s")
    if release is None and snapshot.get("free_time_s") is not None and snapshot.get("time_s") is not None:
        release = snapshot["time_s"]-snapshot["free_time_s"]
    if release is None:
        release = next((r["time_s"] for r in snapshot.get("history") or [] if not r.get("hold", True)), None)
    return release


def sitl_figures(snapshot):
    truth = snapshot["history"]
    figures = {}
    origin = flight_time_origin(snapshot)
    zero = origin or 0.
    times = [r["time_s"]-zero for r in truth]
    for axis, index in [("Roll", 6), ("Pitch", 7)]:
        go, figure = _flight_deck_base_figure(f"{axis} · FDM")
        figure.add_trace(go.Scatter(x=times,
                                   y=[math.degrees(r["state"][index]) for r in truth],
                                   name="FDM", line={"color": "#38bdf8"}))
        figure.update_yaxes(title_text="deg")
        figures[axis] = figure
    go, figure = _flight_deck_base_figure("Firmware outputs")
    for index, label in enumerate(["Aileron", "Elevator", "Pusher", "Rudder"]):
        figure.add_trace(go.Scatter(x=times,
                                   y=[r["pwm"][index] for r in truth], name=label))
    figure.update_yaxes(title_text="PWM · µs")
    figures["PWM"] = figure
    go, figure = _flight_deck_base_figure("Physical elevator")
    figure.add_trace(go.Scatter(x=times,
                               y=[math.degrees(r["control"]["de"]) for r in truth],
                               name="Elevator · deg", line={"color": "#fbbf24"}))
    manifest = snapshot.get("manifest")
    if manifest is not None:
        trim_control = (manifest.get("forward_trim") or {}).get("control", manifest["trim_control"])
        figure.add_hline(y=math.degrees(trim_control["de"]),
                         line_dash="dot", line_color="#94a3b8", annotation_text="Trim")
    figure.update_yaxes(title_text="Surface angle · deg")
    figures["Elevator"] = figure
    go, figure = _flight_deck_base_figure("FDM speed / altitude")
    for label, values in [("TAS · m/s", [math.sqrt(sum(v*v for v in r["state"][3:6])) for r in truth]),
                          ("Altitude · m", [-r["state"][2] for r in truth])]:
        figure.add_trace(go.Scatter(x=times, y=values, name=label))
    figures["Energy"] = figure
    go, figure = _flight_deck_base_figure("Angle of attack · FDM")
    figure.add_trace(go.Scatter(x=times,
        y=[math.degrees(math.atan2(r["state"][5], r["state"][3]))
           if math.sqrt(sum(v*v for v in r["state"][3:6])) >= 1 else None for r in truth], name="Alpha"))
    figure.update_yaxes(title_text="deg · hidden below 1 m/s")
    figures["Alpha"] = figure
    go, figure = _flight_deck_base_figure("Pitch rate · FDM")
    figure.add_trace(go.Scatter(x=times,
        y=[math.degrees(r["state"][10]) for r in truth], name="q"))
    figure.update_yaxes(title_text="deg/s")
    figures["Pitch rate"] = figure
    go, figure = _flight_deck_base_figure("Pusher command · firmware output mapped to FDM")
    figure.add_trace(go.Scatter(x=times,
        y=[100*r["control"]["throttle"] if r["control"].get("throttle") is not None else None
           for r in truth], name="Pusher"))
    figure.update_yaxes(title_text="Actuator command · %", range=[0, 105])
    figures["Pusher"] = figure
    go, figure = _flight_deck_base_figure("Surface travel · physical FDM inputs")
    for name, field in [("Elevator", "de"), ("Aileron", "da"), ("Rudder", "dr")]:
        figure.add_trace(go.Scatter(x=times,
            y=[math.degrees(r["control"][field]) if r["control"].get(field) is not None else None
               for r in truth], name=name))
    limits = (snapshot.get("manifest") or {}).get("control_model", {}).get("assumed_travel_deg", {})
    for angle in sorted({v for value in limits.values() for v in (-value, value)}):
        figure.add_hline(y=angle, line_dash="dot", line_color="#94a3b8", annotation_text="Assumed limit")
    figure.update_yaxes(title_text="deg")
    figures["Surface travel"] = figure
    go, figure = _flight_deck_base_figure("Body rates · FDM")
    for index, name in [(9, "p · roll"), (10, "q · pitch"), (11, "r · yaw")]:
        figure.add_trace(go.Scatter(x=times, y=[math.degrees(r["state"][index]) for r in truth], name=name))
    figure.update_yaxes(title_text="deg/s")
    figures["Body rates"] = figure
    go, figure = _flight_deck_base_figure("Heading change / sideslip · FDM")
    initial_yaw = (manifest or {}).get("initial_state", [0]*12)[8]
    heading = []; previous = initial_yaw; accumulated = 0
    for row in truth:
        accumulated += math.remainder(row["state"][8]-previous, 2*math.pi)
        heading.append(math.degrees(accumulated)); previous = row["state"][8]
    figure.add_trace(go.Scatter(x=times, y=heading, name="Heading change"))
    figure.add_trace(go.Scatter(x=times,
        y=[math.degrees(math.asin(max(-1, min(1, r["state"][4]/math.sqrt(sum(v*v for v in r["state"][3:6]))))))
           if math.sqrt(sum(v*v for v in r["state"][3:6])) >= 1 else None for r in truth], name="Sideslip β"))
    figure.update_yaxes(title_text="deg · measured β=0 only")
    figures["Yaw / sideslip"] = figure
    go, figure = _flight_deck_base_figure("Doublet input · increment from trim")
    plan = (manifest or {}).get("doublet_plan")
    release = origin
    if plan and release is not None:
        scale = 100 if plan["unit"] == "fraction" else 1
        figure.add_trace(go.Scatter(x=times,
            y=[None if r["hold"] else scale*doublet_delta(plan, r["time_s"]-release) for r in truth],
            name="Requested", line={"dash": "dot", "shape": "hv"}))
        if plan.get("regime") == "hover":
            pilot = snapshot.get("pilot_inputs", [])
            figure.add_trace(go.Scatter(x=[r["time_s"]-zero for r in pilot],
                y=[None if r["time_s"] < release or r["channels"][plan["channel"]-1] is None or
                   (plan["name"] == "throttle_doublet" and not plan["start_s"]-1 <= r["time_s"]-release < plan["recovery_at_s"]) else
                   scale*(r["channels"][plan["channel"]-1]-plan["neutral_pwm"])/plan["pwm_per_unit"] for r in pilot],
                name="Measured pilot input · RC", line={"shape": "hv"}))
            if plan["unit"] == "deg":
                figure.add_trace(go.Scatter(x=times,
                    y=[None if r["hold"] else math.degrees(r["state"][6+plan["axis"]]) for r in truth], name="Actual attitude"))
            figure.update_layout(title_text=plan["label"] + " · command and RC readback")
            figure.update_yaxes(title_text="Attitude · deg" if plan["unit"] == "deg" else "Pilot input percentage points")
        else:
            figure.add_trace(go.Scatter(x=times,
                y=[None if r["hold"] else scale*((math.degrees(r["control"][plan["surface"]]) if plan["unit"] == "deg" else r["control"][plan["surface"]])-plan["trim_value"]) for r in truth], name="Actual"))
            figure.update_yaxes(title_text="Throttle percentage points" if scale == 100 else "Surface degrees")
        figure.update_xaxes(range=[plan["start_s"]-2, plan["recovery_at_s"]+3])
    figures["Doublet input"] = figure
    if (manifest or {}).get("parameter_overlay", {}).get("Q_ENABLE") or (manifest or {}).get("config", {}).get("flight_regime") == "hover":
        go, figure = _flight_deck_base_figure("Quad lift motors · firmware outputs")
        for index, label in [(4, "FR · M1"), (7, "RR · M4"), (5, "RL · M2"), (6, "FL · M3")]:
            # Zero means no PWM signal during startup; retain raw records but
            # leave a gap on the motor plot instead of stretching its scale.
            figure.add_trace(go.Scatter(x=times,
                y=[r["pwm"][index] or None for r in truth], name=label, connectgaps=False))
        figure.update_yaxes(title_text="Lift PWM · µs")
        limits = lift_pwm_limits(manifest)
        if limits is not None:
            low, high = limits
            figure.update_yaxes(range=[low, high], autorange=False,
                                tickmode="array", tickvals=[low, (low+high)/2, high])
            for value, label, position in ((low, "Min", "top right"), (high, "Max", "bottom right")):
                figure.add_hline(y=value, line_dash="dot", line_color="#94a3b8",
                                 annotation_text=f"{label} · {value:g} µs", annotation_position=position)
            outside = sum(1 for row in truth for pwm in row["pwm"][4:8]
                          if pwm != 0 and not low <= pwm <= high)
            if outside:
                figure.add_annotation(x=.5, y=.5, xref="paper", yref="paper", showarrow=False,
                    text=f"{outside} motor samples outside limits · inspect raw recording",
                    font={"color": "#f87171"}, bgcolor="#0f172a")
        figures["PWM"] = figure
    for figure in figures.values():
        figure.update_layout(margin={"l": 42, "r": 14, "t": 72, "b": 34},
                             title={"y": .98, "x": .03})
        if origin is not None:
            figure.update_xaxes(title_text="Flight time [s] · 0 = release, startup negative")
            figure.add_vline(x=0, line_dash="dot", line_color="#475569")
    return figures


def register_sitl_page(ui, manager, render_shell):
    current_models = {name: aircraft_summary(load_aircraft_config(name, aerodynamics="tail_update_2026_10_07"))
                      for name in AIRCRAFT}
    @ui.page("/embed")
    @ui.page("/")
    @ui.page("/sitl")
    def sitl_page():
        render_shell(ui)
        preset = FlightSetup(**(manager.snapshot(include_history=False).get("flight_setup") or {}))
        with ui.dialog() as history_dialog, ui.card().classes("fd-panel text-slate-200 w-full max-w-xl p-5 gap-4"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label("Recorded flights").classes("text-lg font-semibold")
                ui.button(icon="close", on_click=history_dialog.close).props('flat round dense aria-label="Close recordings"')
            ui.label("Review the last physical state, messages and plots. Recording controls are read-only.").classes("fd-muted text-sm")
            recording = ui.select(manager.completed_runs(), label="Recorded flights").props("dark filled dense").classes("w-full")
            open_recording = ui.button("Open recording", icon="history")
        embedded = ui.context.client.request.url.path == "/embed"
        if embedded:
            ui.add_body_html('''<script>
            (() => {
              let previous = 0;
              const origins = ['https://www.arrowair.com', 'https://arrowair.com', location.origin];
              if (window.parent === window) return;
              let parentOrigin;
              try { parentOrigin = new URL(document.referrer || location.origin).origin; } catch { return; }
              if (!origins.includes(parentOrigin)) return;
              const connect = () => {
                const root = document.querySelector('.fd-embed');
                if (!root) { requestAnimationFrame(connect); return; }
                new ResizeObserver(() => {
                  const height = Math.ceil(root.getBoundingClientRect().height);
                  if (height === previous) return;
                  previous = height;
                  window.parent.postMessage({type:'arrow-flight-deck:size',height}, parentOrigin);
                }).observe(root);
              };
              connect();
            })();
            </script>''')
        with ui.column().classes("fd-console fd-workspace w-full gap-4" + (" fd-embed" if embedded else "")):
            with ui.row().classes("w-full items-center justify-between fd-intro"):
                with ui.column().classes("gap-1"):
                    ui.label("PROJECT SPEARHEAD / SIMULATION").classes("fd-eyebrow")
                    ui.label("Flight Deck").classes("fd-title")
                    ui.label("Choose a flight. Follow the aircraft. Understand the result.").classes("fd-muted text-sm")
                    with ui.row().classes("items-center gap-1 fd-firmware-note"):
                        ui.icon("memory", size="16px")
                        ui.label("Ardupilot 4.7.1")
                        info_note(ui, "ArduPilot firmware information", "Runs the official ArduPlane 4.7.1 executable on this computer using Software-in-the-Loop (SITL). ArduPilot runs its state estimator, flight modes and controllers, and supplies the motor and servo outputs. The Spearhead model supplies simulated sensors and aircraft physics. Binary identity and firmware version are checked before flight; the exact commit and binary hash are recorded. The aircraft model remains provisional.")
                    cg_summary = ui.label("").classes("fd-muted text-xs")
                with ui.column().classes("items-end gap-1"):
                    status = ui.label("READY").classes("fd-session-status")
                    progress = ui.label("").classes("fd-muted text-xs")
                    with ui.row().classes("items-center gap-1") as pace_panel:
                        pace_readout = ui.label("").classes("fd-muted text-xs")
                        info_note(ui, "Simulation speed information", "1.00× real time means one simulated second takes one real second. Below 1.00×, physics is advancing more slowly. This measures the recent simulation clock against wall time; it does not change physics or controller timing.")
                    history_button = ui.button("Recorded flights", icon="history", on_click=history_dialog.open).props("flat dense color=blue-grey-3").classes("text-xs")
            with ui.element("div").classes("fd-workflow w-full").props('aria-label="Flight workflow"'):
                workflow_steps = []
                for number, title, detail in [(1, "Plan", "Choose your flight"), (2, "Fly", "Follow live progress"), (3, "Review", "Understand the result")]:
                    with ui.row().classes("fd-workflow-step items-center gap-3") as step:
                        ui.label(str(number)).classes("fd-step-number")
                        with ui.column().classes("gap-0"):
                            ui.label(title).classes("font-medium text-sm")
                            ui.label(detail).classes("fd-muted text-xs fd-step-detail")
                    workflow_steps.append(step)
            with ui.expansion("Plan your flight", value=True, icon="flight_takeoff").classes("w-full fd-panel fd-plan").props('header-class="fd-plan-header"') as plan_panel:
                with ui.column().classes("w-full gap-3"):
                    with ui.element("div").classes("fd-setup-grid"):
                        aircraft = ui.select(AIRCRAFT, value=preset.aircraft, label="Aircraft").props("dark outlined dense").classes("w-full")
                        initial_point = ui.select(POINTS[preset.aircraft], value=preset.initial_point, label="Start position").props("dark outlined dense").classes("w-full")
                        maneuver = ui.select(maneuver_options(preset.aircraft, preset.initial_point),
                                             value=preset.maneuver, label="Flight test").props("dark outlined dense").classes("w-full")
                        duration = ui.number("Flight duration (s)", value=preset.duration_s, min=15, max=600, step=5).props("dark outlined dense").classes("w-full")
                    with ui.row().classes("items-center gap-1"):
                        plan_summary = ui.label("").classes("text-sm text-slate-300")
                        setup_hint = info_note(ui, "Flight setup information")
                    with ui.expansion("Starting conditions & model settings", icon="tune").classes("w-full fd-advanced") as advanced_settings:
                        ui.label("Adjust the starting point or compare saved model assumptions. Controller parameters stay frozen.").classes("fd-muted text-xs mb-3")
                        with ui.element("div").classes("fd-condition-grid w-full"):
                            altitude = ui.number("Initial altitude (m)", value=preset.altitude_m, min=6, max=110, step=1).props("dark outlined dense").classes("w-full")
                            speed = ui.number("Initial airspeed (m/s)", value=preset.speed_m_s, min=0, max=30, step=1).props("dark outlined dense").classes("w-full")
                        with ui.row().classes("w-full gap-3 mt-3") as transition_choices:
                            cg_case = ui.select(TRANSITION_CGS, value=preset.cg_case, label="Transition CG").props("dark outlined dense").classes("grow min-w-0")
                            transition_parameters = ui.select(TRANSITION_PARAMETERS, value=preset.transition_parameters, label="Transition parameters").props("dark outlined dense").classes("grow min-w-0")
                        with ui.row().classes("w-full items-center gap-1 mt-2") as cg_card_row:
                            cg_card_label = ui.label("").classes("text-xs text-sky-100")
                            info_note(ui, "CG summary information", "For the selected CG with the October AeroDB: stick-fixed static margin at 22 m/s trim, elevator trim at 20/22/24 m/s with ±30° travel, and the busiest lift rotor's thrust as a fraction of its thrust at 95 % spin in static hover. Elevator beyond ±15° relies on extrapolated control data. These are model estimates, not a measured CG envelope.")
                        ui.label("Wind and failures · times count from release (ground starts: from motor spin-up)").classes("fd-muted text-xs mt-3")
                        with ui.element("div").classes("fd-condition-grid w-full") as condition_choices:
                            wind_speed = ui.number("Wind (m/s)", value=preset.wind_speed_m_s, min=0, max=15, step=1).props("dark outlined dense").classes("w-full")
                            wind_from = ui.number("Wind from (°)", value=preset.wind_from_deg, min=0, max=359, step=15).props("dark outlined dense").classes("w-full")
                            gust = ui.number(f"Gust at {GUST_START_S:g} s (m/s)", value=preset.gust_m_s, min=0, max=10, step=1).props("dark outlined dense").classes("w-full")
                        with ui.row().classes("w-full gap-3 mt-3"):
                            failure = ui.select(FAILURES, value=preset.failure, label="Failure").props("dark outlined dense").classes("grow min-w-0")
                            failure_time = ui.number("Failure at (s)", value=preset.failure_time_s, min=0, max=600, step=1).props("dark outlined dense").classes("w-40")
                    with ui.row().classes("w-full items-center gap-3 fd-flight-actions"):
                        start = ui.button("Start flight", icon="play_arrow")
                        ui.space()
                        ui.label("Startup runs before the flight clock begins.").classes("fd-muted text-xs")
                        info_note(ui, "Flight timing information", "Startup usually takes 1–2 minutes before the flight clock starts. Flight duration is simulated time; real time depends on computer speed. Closing this tab keeps the simulator running. Stop flight ends the run.")
                    setup_error = ui.label("").classes("text-amber-200 text-sm").props('role="status"')
            with ui.row().classes("w-full fd-startup items-center gap-3 p-3") as startup_panel:
                ui.spinner(size="24px", color="light-blue-3")
                with ui.column().classes("grow gap-1"):
                    startup_title = ui.label("").classes("text-sm font-semibold text-sky-100").props('role="status" aria-live="polite"')
                    startup_detail = ui.label("").classes("text-xs text-slate-300")
                with ui.column().classes("gap-1 items-end fd-startup-time"):
                    startup_elapsed = ui.label("").classes("text-xs font-mono text-slate-300")
                    startup_estimate = ui.label("").classes("text-xs text-sky-200")
            with ui.column().classes("w-full fd-panel fd-flight-status p-4 gap-3") as flight_panel:
                with ui.row().classes("w-full items-start justify-between gap-3"):
                    with ui.column().classes("grow gap-1"):
                        flight_title = ui.label("").classes("fd-section-title").props('role="status" aria-live="polite"')
                        flight_detail = ui.label("").classes("text-sm text-slate-300")
                    display_scope = ui.label("").classes("fd-eyebrow")
                flight_duration = ui.label("").classes("fd-muted text-xs")
                flight_progress = ui.linear_progress(value=0, color="primary", show_value=False).classes("w-full")
                with ui.row().classes("w-full items-center gap-3"):
                    with ui.row().classes("items-center gap-2 mr-3"):
                        current_mode = ui.label("Flight mode: —").classes("text-sm font-medium text-sky-100")
                        mode_help = info_note(ui, "Flight mode information", MODE_HELP[18])
                    arm_state = ui.label("ARMING UNKNOWN").classes("fd-arm-badge fd-arm-unknown").props('role="status" aria-live="polite"')
                    ui.space()
                    stop = ui.button("Stop flight", icon="stop").props("outline color=red-4")
                    plan_again = ui.button("Plan another flight", icon="add", on_click=lambda: plan_panel.set_value(True)).props("outline")
                    with ui.row().classes("items-center gap-3") as mode_controls:
                        mode_select = ui.select(QUADPLANE_MODES, value=18, label="Requested mode").props("dark outlined dense").classes("w-64")
                        apply_mode = ui.button("Change mode", icon="swap_horiz")
                    command_status = ui.label("").classes("text-sm text-sky-200")
                control_hint = ui.label("").classes("text-xs text-slate-300")
                live_alert = ui.label("").classes("w-full fd-alert text-amber-200 text-sm").props('role="status" aria-live="polite"')
            with ui.row().classes("w-full fd-health gap-x-5 gap-y-1 items-center px-1") as readiness_panel:
                link_status = ui.label("Link: offline").classes("text-slate-300")
                gps_status = ui.label("GPS: unknown").classes("text-slate-300")
                ekf_status = ui.label("EKF: unknown").classes("text-slate-300")
                system_status = ui.label("System: unknown").classes("text-slate-300")
                battery_status = ui.label("Simulated battery: unknown").classes("text-slate-300")
                info_note(ui, "Aircraft readiness information", "Connection shows whether fresh simulator and firmware data are arriving. GPS provides position; navigation readiness comes from the firmware's EKF state estimator. Check messages for the reasons behind failed pre-arm checks. Battery values appear only when the simulated firmware supplies them.")
            error = ui.label("").classes("w-full fd-alert fd-error text-red-200 text-sm whitespace-pre-wrap").props('role="alert"')
            warning = ui.label("").classes("w-full fd-alert text-amber-200 text-xs").props('role="status"')
            with ui.element("div").classes("fd-instruments fd-flight-instruments w-full gap-3"):
                with ui.column().classes("fd-panel p-3 w-full min-w-0 gap-2"):
                    with ui.row().classes("items-center gap-1"):
                        ui.label("Aircraft attitude").classes("fd-section-title")
                        info_note(ui, "Attitude display information", "The horizon shows pitch (nose up or down) and roll (bank left or right) from the physical simulation. Airspeed is speed through the air; altitude is height above the simulated reference ground.")
                    fdm_display = ui.html(attitude_html(pitch_deg=0, roll_deg=0)).classes("w-full")
                    with ui.column().classes("w-full items-center justify-center gap-2 rounded-lg bg-slate-950").style("height:330px") as attitude_empty:
                        ui.icon("flight", size="40px").classes("text-slate-600")
                        ui.label("Ready for your next flight").classes("text-slate-300")
                        ui.label("Start a session or open a recording").classes("fd-muted text-xs")
                    truth_readout = ui.label("Waiting for flight data").classes("font-mono text-sm text-slate-300")
                    ui.label("Pitch & roll from the simulated aircraft").classes("fd-muted text-xs")
                with ui.column().classes("fd-panel fd-values-panel p-3 w-full min-w-0 gap-2") as values_panel:
                    with ui.row().classes("w-full items-center justify-between gap-1"):
                        ui.label("Flight values").classes("fd-section-title")
                        info_note(ui, "Flight values information", "Airspeed, altitude, heading and body rotation rates come from the physical simulation, matching the horizon. Altitude is above the simulated ground datum. Positive climb rate means moving upward in earth axes; p/q/r are body roll/pitch/yaw rates. Pusher command is the firmware output mapped to actuator percentage, not thrust percentage. Lift PWM is the actual firmware output, arranged as front left/front right above rear left/rear right. A dash means unavailable or no PWM signal. Saved and stale readings are labeled; these are not hardware measurements or separate EKF estimates.")
                    values_status = ui.label("WAITING FOR FLIGHT DATA").classes("fd-eyebrow")
                    with ui.element("div").classes("fd-values-grid w-full"):
                        value_readouts = {}
                        for key, label in VALUE_LABELS.items():
                            with ui.column().classes("fd-value gap-0"):
                                ui.label(label).classes("fd-value-label")
                                value_readouts[key] = ui.label("—").classes("fd-value-number")
            with ui.column().classes("w-full fd-panel gap-0"):
                with ui.tabs().props("no-caps inline-label align=left active-color=light-blue-3").classes("w-full") as detail_tabs:
                    overview_tab = ui.tab("Overview", icon="fact_check")
                    messages_tab = ui.tab("Messages", icon="chat_bubble_outline")
                    plots_tab = ui.tab("Flight plots", icon="show_chart")
                    model_tab = ui.tab("Technical details", icon="tune")
                    compare_tab = ui.tab("Compare with aircraft", icon="compare_arrows")
                with ui.tab_panels(detail_tabs, value=overview_tab).props("animated=false").classes("w-full bg-transparent text-slate-200"):
                    with ui.tab_panel(overview_tab).classes("p-4 gap-3"):
                        ui.label("Understanding this flight").classes("fd-section-title")
                        overview_context = ui.label("").classes("text-sm text-slate-300")
                        report = ui.label("").classes("text-sky-100 text-sm")
                        with ui.element("div").classes("fd-result-metrics w-full") as metrics_panel:
                            metric_readouts = {}
                            metric_labels = {}
                            for key, label in [("longest_forward_without_lift_s", "Accepted forward interval"),
                                               ("peak_abs_pitch_deg", "Peak pitch magnitude"),
                                               ("altitude_loss_m", "Maximum altitude loss"),
                                               ("hover_peak_abs_pitch_deg", "Initial hover peak pitch")]:
                                with ui.column().classes("gap-1"):
                                    metric_labels[key] = ui.label(label).classes("fd-muted text-xs")
                                    metric_readouts[key] = ui.label("—").classes("text-lg")
                        metrics_note = ui.label("Forward acceptance covers a sustained interval. The pitch and altitude extremes cover the whole released flight; initial hover covers its first 10 s.").classes("fd-muted text-xs")
                        with ui.row().classes("items-center gap-1"):
                            ui.label("How the test is assessed").classes("fd-muted text-sm")
                            acceptance_note = info_note(ui, "Test assessment information")
                    with ui.tab_panel(messages_tab).classes("p-3 gap-3"):
                        with ui.row().classes("w-full items-center justify-between"):
                            message_count = ui.label("No firmware messages yet.").classes("fd-muted text-xs")
                            message_filter = ui.select(["All messages", "Warnings and errors"], value="Warnings and errors", label="Show").props("dark filled dense").classes("w-56")
                        message_table = ui.table(columns=[
                            {"name": "time", "label": "Time", "field": "time", "align": "left"},
                            {"name": "level", "label": "Severity", "field": "level", "align": "left"},
                            {"name": "source", "label": "Source", "field": "source", "align": "left"},
                            {"name": "text", "label": "Message", "field": "text", "align": "left", "style": "white-space: normal"},
                        ], rows=[], row_key="id", pagination=5).props("dark flat dense wrap-cells").classes("w-full")
                        message_table.add_slot('body-cell-level', '<q-td :props="props"><q-badge :color="props.row.severity <= 3 ? \'red\' : props.row.severity === 4 ? \'amber-9\' : \'blue-grey\'">{{ props.value }}</q-badge></q-td>')
                    with ui.tab_panel(plots_tab).classes("p-3"):
                        with ui.element("div").classes("grid grid-cols-1 lg:grid-cols-2 w-full gap-3"):
                            from .visible_plotly import VisiblePlotly
                            plots = {name: VisiblePlotly(figure).classes("w-full")
                                     for name, figure in sitl_figures(manager.snapshot()).items()}
                            plots["Doublet input"].set_visibility(False)
                    with ui.tab_panel(compare_tab).classes("p-4 gap-3"):
                        ui.label("Compare this run with the real aircraft").classes("fd-section-title")
                        ui.label("Give a DataFlash log (.BIN) to overlay lift outputs, attitude and altitude from takeoff and to compare parameters, "
                                 "or a .param file for parameters only. Open the SITL run to compare first (current or Recorded flights).").classes("fd-muted text-xs")
                        with ui.row().classes("w-full items-end gap-3"):
                            compare_path = ui.input("File path on this computer", placeholder="~/Downloads/00000017.BIN").props("dark outlined dense").classes("grow min-w-0")
                            compare_flight = ui.number("Flight in log (0 = last)", value=0, min=0, max=50, step=1).props("dark outlined dense").classes("w-48")
                            compare_button = ui.button("Compare", icon="compare_arrows")
                        compare_upload = ui.upload(label="…or upload a .BIN / .param file", auto_upload=True,
                                                   max_file_size=600_000_000).props("dark flat bordered accept=.bin,.BIN,.param,.parm,.params,.txt").classes("w-full")
                        compare_status = ui.label("").classes("text-sm text-sky-100")
                        with ui.element("div").classes("grid grid-cols-1 lg:grid-cols-2 w-full gap-3"):
                            from .visible_plotly import VisiblePlotly
                            compare_plots = {name: VisiblePlotly({"data": [], "layout": {}}).classes("w-full") for name in ("lift", "attitude", "altitude")}
                        for plot in compare_plots.values():
                            plot.set_visibility(False)
                        compare_table = ui.table(columns=[
                            {"name": "name", "label": "Parameter", "field": "name", "align": "left", "sortable": True},
                            {"name": "aircraft", "label": "Aircraft", "field": "aircraft", "align": "right"},
                            {"name": "simulation", "label": "SITL run", "field": "simulation", "align": "right"},
                            {"name": "status", "label": "", "field": "status", "align": "left"},
                        ], rows=[], row_key="name", pagination=15).props("dark flat dense").classes("w-full")
                        compare_table.set_visibility(False)
                    with ui.tab_panel(model_tab).classes("p-4 gap-3"):
                        ui.label("Simulator & data versions").classes("font-semibold")
                        aircraft_readout = ui.label("").classes("fd-muted text-sm break-all")
                        ui.separator()
                        with ui.row().classes("items-center gap-1"):
                            ui.label("Aircraft & model assumptions").classes("font-semibold")
                            mode_limits = info_note(ui, "Aircraft model information")
                        mode_scope = ui.label("").classes("text-sm text-slate-300")
                        cg_readout = ui.label("").classes("text-sm font-mono text-slate-300")
                        trim_readout = ui.label("").classes("text-sm text-sky-100")
                        with ui.row().classes("items-center gap-1"):
                            ui.label("Control travel & measured coverage").classes("text-sm text-slate-300")
                            travel_readout = info_note(ui, "Control coverage information")
                        aero_readout = ui.label("").classes("text-sm font-mono text-slate-300")
                        ui.separator()
                        ui.label("Run information").classes("font-semibold")
                        connection = ui.label("Firmware not running").classes("font-mono text-xs text-slate-300")
                        run_label = ui.label("").classes("fd-muted text-xs break-all")

        launch_pending = [False]
        launch_started = [None]

        async def start_clicked():
            launch_pending[0] = True
            launch_started[0] = time.monotonic()
            refresh()
            try:
                await asyncio.to_thread(manager.start_configured, selected_setup())
                plan_panel.value = False
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            finally:
                launch_pending[0] = False
            refresh()

        def selected_setup():
            return FlightSetup(aircraft.value, initial_point.value, maneuver.value,
                               duration.value, altitude.value, speed.value,
                               cg_case.value, transition_parameters.value,
                               wind_speed.value or 0, wind_from.value or 0, gust.value or 0,
                               failure.value, failure_time.value or 0).validate()

        syncing_setup = [False]
        setup_run = [None]
        displayed_directory = [None]

        async def compare_file(path):
            from nicegui import run
            from ..flight_log import compare_log_with_run
            from pathlib import Path as _Path
            directory = displayed_directory[0]
            if directory is None or not (_Path(directory) / "truth.jsonl").is_file():
                compare_status.text = "Open a SITL run first (start one or pick it from Recorded flights)."
                return
            path = _Path(path).expanduser()
            if not path.is_file():
                compare_status.text = f"File not found: {path}"
                return
            flight = int(compare_flight.value or 0)
            compare_button.disable()
            compare_status.text = f"Reading {path.name}… large logs take a minute."
            try:
                result = await run.cpu_bound(compare_log_with_run, str(path), str(directory), flight - 1 if flight else -1)
            except Exception as exc:
                compare_status.text = f"Could not compare {path.name}: {exc}"
                return
            finally:
                compare_button.enable()
            compare_status.text = comparison_text(result)
            figures = comparison_figures(result)
            for name, plot in compare_plots.items():
                if name in figures:
                    plot.figure = figures[name]
                    plot.update()
                plot.set_visibility(name in figures)
            compare_table.rows = result["parameter_diff"]
            compare_table.update()
            compare_table.set_visibility(True)

        async def compare_uploaded(event):
            target_dir = manager.root / "aircraft_logs"
            target_dir.mkdir(parents=True, exist_ok=True)
            target = target_dir / Path(event.file.name).name
            await event.file.save(target)
            compare_path.value = str(target)
            await compare_file(target)
        cg_card_key = [None]

        async def update_cg_card():
            from nicegui import run
            from ..cg_card import cg_card, cg_card_text
            if aircraft.value != "pt15":
                return
            key = (aircraft.value, TRANSITION_CG_X[cg_case.value], TRANSITION_CG_Z)
            if key == cg_card_key[0]:
                return
            cg_card_key[0] = key
            cg_card_label.text = "Computing static margin, elevator trim and hover headroom for this CG…"
            try:
                card = await run.cpu_bound(cg_card, *key)
                text = cg_card_text(card)
            except Exception as exc:
                text = f"CG summary unavailable: {exc}"
            if cg_card_key[0] == key:
                cg_card_label.text = text

        def schedule_cg_card():
            from nicegui import background_tasks
            background_tasks.create(update_cg_card(), name="cg-card")

        def aircraft_changed():
            if syncing_setup[0]:
                return
            syncing_setup[0] = True
            try:
                hover = aircraft.value == "pt15"
                initial_point.set_options(POINTS[aircraft.value], value="hover20" if hover else "trim28")
                maneuver.set_options(maneuver_options(aircraft.value, initial_point.value), value="interactive")
                altitude.value, speed.value = (20, 0) if hover else (100, 28)
            finally:
                syncing_setup[0] = False
            refresh()

        def point_changed():
            if syncing_setup[0]:
                return
            syncing_setup[0] = True
            try:
                point = initial_point.value
                options = maneuver_options(aircraft.value, point)
                maneuver.set_options(options, value=maneuver.value if maneuver.value in options else next(iter(options)))
                if point == "ground":
                    altitude.value = min(30, max(3, altitude.value or 20))
                    speed.value = 0
                elif forward_start(point):
                    if not speed.value:
                        altitude.value, speed.value = 100, 24
                else:
                    if speed.value:
                        altitude.value = 20
                    speed.value = 0
                if point.startswith("hover"):
                    altitude.value = int(point[5:])
                elif point.startswith("trim"):
                    speed.value = int(point[4:])
            finally:
                syncing_setup[0] = False
            refresh()

        def conditions_changed():
            if syncing_setup[0]:
                return
            point = initial_point.value
            if point.startswith("hover") and altitude.value != int(point[5:]):
                initial_point.value = "custom_hover"
            elif point.startswith("trim") and speed.value != int(point[4:]):
                initial_point.value = "custom_trim"
            refresh()

        def maneuver_changed():
            if syncing_setup[0]:
                return
            minimum = MIN_DURATION.get(maneuver.value, 15)
            if duration.value is not None and duration.value < minimum:
                duration.value = minimum
            refresh()

        async def stop_clicked():
            stop.disable()
            await asyncio.to_thread(manager.stop)
            refresh()

        start.on_click(start_clicked)
        stop.on_click(stop_clicked)
        async def mode_clicked():
            try:
                manager.set_mode(mode_select.value)
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            refresh()
        apply_mode.on_click(mode_clicked)
        mode_select.on_value_change(lambda e: mode_help.set_text(MODE_HELP.get(e.value, "Choose a flight mode.")))
        async def recording_clicked():
            if recording.value is None:
                return
            try:
                await asyncio.to_thread(manager.load_completed, recording.value)
                history_dialog.close()
                plan_panel.value = False
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            refresh()
        open_recording.on_click(recording_clicked)
        last_plot = [None, None]
        last_messages = [None]
        last_attitude = [None]
        pace = SimulationPace()

        def refresh():
            view = manager.snapshot(include_history=False)
            active = view["active"]
            clock_speed = pace.update(view, time.monotonic())
            pace_panel.set_visibility(clock_speed is not None)
            if clock_speed is not None:
                pace_readout.text = f"Simulation speed · {clock_speed:.2f}× real time"
            locked = active or launch_pending[0]
            activity = ({"title": "Preparing your flight", "detail": "Validating the setup and starting the local simulator.",
                         "elapsed_s": time.monotonic()-launch_started[0], "estimate": "Usually 1–2 min before flight · may take longer"}
                        if launch_pending[0] else startup_activity(view))
            startup_panel.set_visibility(activity is not None)
            if activity:
                startup_title.text = activity["title"]
                startup_detail.text = activity["detail"]
                startup_elapsed.text = f"Elapsed {activity['elapsed_s']:.0f} s"
                if not launch_pending[0] and view.get("time_s") is not None:
                    startup_elapsed.text += f" · sim clock {view['time_s']:.1f} s"
                startup_estimate.text = activity["estimate"]
            start.text = "Starting…" if activity else "Start flight"
            plan_panel.props('label="Plan your next flight"' if not active and view.get("run_id") else 'label="Plan your flight"')
            # Show the actual recipe when a shared session or recording changes.
            # After synchronization, completed-run fields can be edited for a new flight.
            if view.get("flight_setup") and setup_run[0] != view.get("run_id"):
                recipe = FlightSetup(**view["flight_setup"]).validate()
                syncing_setup[0] = True
                try:
                    aircraft.value = recipe.aircraft
                    initial_point.set_options(POINTS[recipe.aircraft], value=recipe.initial_point)
                    maneuver.set_options(maneuver_options(recipe.aircraft, recipe.initial_point), value=recipe.maneuver)
                    altitude.max = 1000 if recipe.is_forward else 110
                    speed.min = 22 if recipe.is_forward else 0
                    duration.min = MIN_DURATION.get(recipe.maneuver, 15)
                    duration.value = recipe.duration_s
                    altitude.value = recipe.altitude_m
                    speed.value = recipe.speed_m_s
                    cg_case.value = recipe.cg_case
                    transition_parameters.value = recipe.transition_parameters
                    wind_speed.value, wind_from.value, gust.value = recipe.wind_speed_m_s, recipe.wind_from_deg, recipe.gust_m_s
                    failure.value, failure_time.value = recipe.failure, recipe.failure_time_s
                    setup_run[0] = view.get("run_id")
                    plan_panel.value = False
                finally:
                    syncing_setup[0] = False
            prototype_setup = aircraft.value == "pt15"
            hover_setup = not forward_start(initial_point.value)
            prototype_forward = prototype_setup and not hover_setup
            ground_setup = initial_point.value == "ground"
            try:
                plan_summary.text = flight_plan_text(selected_setup())
                setup_error.text = ""
            except ValueError as exc:
                setup_error.text = str(exc)
                plan_summary.text = "Check the flight settings below."
            setup_error.set_visibility(bool(setup_error.text))
            start.set_enabled(not locked and not setup_error.text)
            stop.set_enabled(active)
            stop.set_visibility(active)
            for field in (aircraft, initial_point, maneuver, duration, altitude, cg_case, transition_parameters,
                          wind_speed, wind_from, gust, failure, failure_time):
                field.set_enabled(not locked)
            forward_choices = aircraft.value == "pt15"
            transition_choices.set_visibility(forward_choices)
            cg_card_row.set_visibility(forward_choices)
            cg_case.label = "Simulation CG (target, not measured)"
            transition_parameters.label = "Forward parameters"
            caption = (TRANSITION_CGS[cg_case.value] + " · " + TRANSITION_PARAMETERS[transition_parameters.value]
                       if forward_choices else "Altitude and airspeed")
            if wind_speed.value or gust.value:
                caption += f" · wind {wind_speed.value or 0:g} m/s" + (f" + {gust.value:g} m/s gust" if gust.value else "")
            if failure.value != "none":
                caption += " · " + FAILURES[failure.value]
            advanced_settings.props(f'caption="{caption}"')
            speed.set_enabled(not locked and not hover_setup)
            # Keep numeric limits numeric; props strings also reach server-side sanitization.
            altitude.min = 3 if ground_setup else 6
            altitude.max = 30 if ground_setup else 110 if hover_setup else 1000
            altitude.label = "Takeoff altitude (m)" if ground_setup else "Initial altitude (m)"
            duration.label = "Hover time (s)" if ground_setup else "Flight duration (s)"
            speed.min = 0 if hover_setup else 22
            speed.max = 26 if prototype_forward else 30
            duration.min = MIN_DURATION.get(maneuver.value, 15)
            setup_hint.text = ("Air-start hover · KDE lift · MAD 150KV / APC 20×12 pusher estimate." if hover_setup else
                               "Air-start at six-axis trim · provisional pusher · Spearhead hover awaits MAD performance data.")
            if prototype_forward:
                setup_hint.text = "PT1.5 forward air start at computed six-axis trim for the selected speed and simulation CG. ±30° proposed travel; frozen forward parameters. Lift motors remain off; fixed-wing modes are available. Physical CG is unconfirmed."
            if maneuver.value == "interactive":
                setup_hint.text += " Interactive means you can request mode changes. Pilot sticks remain neutral; no scripted input is applied."
            if maneuver.value == "back_transition":
                setup_hint.text = "Start in the selected FW trim, settle for 10 s, then request QLOITER. Firmware stops the pusher, brakes and hands over to lift motors. Requires at least 10 s of confirmed low-speed hover. Frozen gains; physical CG/inertia remain provisional."
            if maneuver.value in DOUBLETS and hover_setup:
                setup_hint.text = "Hover pilot doublet: 10 s settling, positive and negative 2 s pulses, then neutral QHOVER recovery. Pitch/roll command attitude; yaw commands stick input. Lift collective uses QSTABILIZE around calculated hover thrust; motor PWM is nonlinear. QHOVER does not hold horizontal position. Frozen logged gains; minimum 45 s."
            if maneuver.value in DOUBLETS and not hover_setup:
                setup_hint.text += f" Doublets use MANUAL only for the physical input, then FBWA for recovery. Each pulse lasts {DOUBLETS[maneuver.value]['pulse_s']:g} s. Throttle amplitude is capped symmetrically by available trim headroom and recorded. Input-delivery checks do not establish stability."
            setup_hint.text += " Simulation time is free-flight time; calibration and arming run first."
            if hover_setup and not ground_setup:
                setup_hint.text += " Hover attitude uses lift motor mixing. ArduPilot also commands ailerons/elevator; rudder is centered in QHOVER/QSTABILIZE. Surface authority depends on airflow. The lift PWM plot uses Q_M_PWM_MIN/MAX; zero startup values mean no signal."
            if ground_setup:
                setup_hint.text = "Start on the ground → VTOL takeoff → timed hover → landing. Climb and landing are additional to hover time. Flat ground and provisional landing gear."
            if maneuver.value == "transition":
                setup_hint.text = ("PT1.5 diagnostic · 10 s hover → FBWA → 5 s pusher ramp. "
                                   "CG uses body axes: x forward, y right, z down; source moments at (0,0,0). "
                                   "The base export is numerically restored to its source origin before transfer to CG. Earlier transition passes used the old reference handling. "
                                   "The default CG is 10 mm forward; z remains a provisional 150 mm below the origin. Corrected-model transition acceptance is recorded separately. Physical CG remains unconfirmed. "
                                   "Frozen parameters; trim is recomputed for the selected CG.")
            if maneuver.value in MIN_DURATION:
                setup_hint.text += f" This maneuver sequence needs at least {MIN_DURATION[maneuver.value]} s."
            setup_hint.text += " Closing this tab keeps the run active; Stop flight ends it."
            recording.set_enabled(not locked)
            open_recording.set_enabled(not locked and recording.value is not None)
            history_button.set_enabled(not locked)
            options = manager.completed_runs()
            if options != recording.options:
                recording.set_options(options)
            status.text = view["status"]
            manifest = view.get("manifest") or {}
            displayed_aircraft = manifest.get("aircraft") or current_models[aircraft.value]
            aircraft_readout.text = "ArduPlane 4.7.1 · October AeroDB · Provisional simulation"
            if manifest:
                aircraft_readout.text += " · Displaying " + manifest.get("aircraft", {}).get("label", "recorded aircraft")
            aircraft_readout.text += " · " + displayed_aircraft.get("aero_database_version", "version unknown")
            cg = displayed_aircraft["cg_body_xyz_m"] if "cg_body_xyz_m" in displayed_aircraft else None
            if not manifest and forward_choices:
                cg = [TRANSITION_CG_X[cg_case.value], 0, TRANSITION_CG_Z]
            origin = displayed_aircraft.get("moment_center_body_xyz_m")
            cg_readout.text = f"Simulation CG body xyz (m, target, not measured): {cg} · Source moment center (m): {origin} · x forward / z down"
            cg_summary.set_visibility(cg is not None and (bool(manifest) or forward_choices))
            if cg is not None:
                cg_summary.text = ("Recorded simulation CG" if manifest and not active else "Simulation CG (target, not measured)") + ": " + \
                    " · ".join(f"{axis} {1000*value:+g} mm" for axis, value in zip("xyz", cg)) + " · body axes"
            is_back = manifest.get("test_kind") == "back_transition"
            is_transition = manifest.get("test_kind") in ("quad_transition", "back_transition")
            metric_labels["longest_forward_without_lift_s"].text = "Accepted hover interval" if is_back else "Accepted forward interval"
            metric_labels["hover_peak_abs_pitch_deg"].text = "Initial forward peak pitch" if is_back else "Initial hover peak pitch"
            metrics_note.text = "Acceptance covers a sustained interval; pitch and altitude extremes cover the whole released flight. The initial peak covers its first 10 s."
            metrics_panel.set_visibility(is_transition)
            metrics_note.set_visibility(is_transition)
            for key, value in transition_metrics(view).items():
                unit = "s" if key.endswith("_s") else "m" if key.endswith("_m") else "°"
                metric_readouts[key].text = "—" if value is None else f"{value:.1f}{'' if unit == '°' else ' '}{unit}"
            forward_trim = manifest.get("forward_trim") or {}
            trim_readout.text = (f"Forward trim at {forward_trim['speed_m_s']:g} m/s: elevator "
                                 f"{math.degrees(forward_trim['control']['de']):+.2f}° · "
                                 f"pusher {100*forward_trim['control']['throttle']:.1f}% · "
                                 f"parameters: {manifest.get('transition_plan', {}).get('parameter_set', manifest.get('forward_parameter_set', 'unknown'))}"
                                 if forward_trim else "Forward trim is computed for the selected speed and CG when the flight initializes.")
            trim_readout.set_visibility(bool(forward_trim) or forward_choices)
            context = manifest if manifest else {"parameter_overlay": {"Q_ENABLE": int(hover_setup or maneuver.value == "back_transition")}}
            mode_options = session_modes(context)
            if mode_select.options != mode_options:
                mode_select.set_options(mode_options, value=18 if 18 in mode_options else 5)
            if not view.get("interactive") and view.get("mode") in mode_options:
                mode_select.value = view["mode"]
            mode_help.text = MODE_HELP.get(mode_select.value, "Choose a flight mode.")
            forward_mode = 18 in mode_options and (mode_select.value in PLANE_MODES or mode_select.value == 25)
            if forward_mode:
                mode_help.text += " PT1.5 electric pusher is estimated; forward-flight results depend on the simulation CG."
            mass = manifest.get("aircraft", {}).get("mass_kg", 21)
            mode_scope.text = (f"PT1.5 Quad-X · {mass:g} kg · four KDE lift motors · neutral pilot sticks" if 18 in mode_options else
                               f"PT1.5 fixed-wing · {mass:g} kg · lift motors off · neutral pilot sticks" if displayed_aircraft.get("id") == "pt15" else
                               "Fixed-wing reference aircraft · neutral pilot sticks")
            mode_limits.text = ("PT1.5 pusher: simulation estimate; physical motor and thrust-line placement are not finalized. CG and inertia remain unconfirmed. Servo travel is proposed at ±30°, with uncalibrated linkage. Earlier forward-CG passes used ±20° travel. "
                                "No ground contact or mission upload. AUTOLAND is unavailable for QuadPlane." if 18 in mode_options else
                                "Air-start simulation: no runway/contact model, mission upload or guided-target entry.")
            if 18 not in mode_options and displayed_aircraft.get("id") == "pt15":
                mode_limits.text = "PT1.5 forward trim air start · lift off · fixed-wing modes only. Pusher hardware, CG and inertia remain provisional; no runway/contact model."
            if manifest.get("ground_contact"):
                mode_limits.text = "AUTO flies the built-in takeoff, hover and landing mission. Flat-ground gear is provisional; forward transition remains unvalidated."
                if mode_select.value in (10, 13, 20, 21):
                    mode_help.text = "This ground session includes physical gear contact and an uploaded VTOL mission. Mode changes are disabled while the automatic sequence runs."
            mode_limits.text += " Simulation CG (target, not measured) and provisional aircraft model. Control data cover β=0 only; flight outside that plane uses endpoint assumptions."
            error.text = view.get("error") or ""
            error.set_visibility(bool(error.text))
            historical = not active and view.get("run_id") is not None
            plan_again.set_visibility(historical)
            mode_name = view.get("mode_name")
            current_mode.text = (f"{MODE_LABELS.get(view.get('mode'), 'Flight mode')} · {mode_name}" if mode_name else "Flight mode: awaiting firmware")
            arm_state.text, arm_tone = arming_badge(view)
            arm_state.classes(replace=f"fd-arm-badge fd-arm-{arm_tone}")
            command = view.get("command") or {}
            pending = command.get("status") in ("queued", "sent", "accepted")
            enabled = (active and view.get("interactive") and view["status"] == "FREE FLIGHT"
                       and not view.get("stale_truth") and not view.get("stale_mavlink") and not view.get("error")
                       and view.get("mode_controls_ready", True))
            mode_select.set_enabled(enabled and not pending)
            apply_mode.set_enabled(enabled and not pending)
            mode_controls.set_visibility(bool(view.get("interactive")))
            command_text = command.get("text", "")
            if command.get("status") in ("rejected", "timeout"):
                target = MODE_NAMES.get(command.get("mode"), "Mode")
                command_text = f"{target} request {command['status']}: {command_text}"
            command_status.text = command_text
            command_status.classes(replace="text-sm text-amber-200" if command.get("status") in ("rejected", "timeout") else "text-sm text-sky-200")
            control_hint.text = ("Choose a mode, then select Change mode. The firmware confirms each request." if enabled else
                                 "Saved recording · displays show the last sample. Edit the plan above to start a new flight." if historical else
                                 "This test controls flight modes automatically." if active and manifest and not view.get("interactive") else
                                 "Calibrating and arming — mode controls unlock when ready." if active else
                                 "Motors cannot arm until a flight is started.")
            stale = view.get("stale_mavlink") or view.get("stale_truth")
            link_status.text = "Connection: recorded" if historical else "Connection: waiting" if active and view.get("truth") is None else "Connection: stale" if active and stale else "Connection: live" if active else "Connection: offline"
            health = view.get("health") or {}
            fix, satellites = health.get("gps_fix"), health.get("satellites")
            gps_status.text = "GPS: " + ("unknown" if fix is None else f"{'3D fix' if fix >= 3 else 'no 3D fix'} · {satellites if satellites is not None else '—'} satellites")
            ekf = health.get("ekf_ok")
            ekf_status.text = "Navigation: " + ("unknown" if ekf is None else "ready" if ekf else "initializing")
            ekf_status.classes(replace="text-amber-200" if ekf is False else "text-slate-300")
            faults = health.get("sensor_faults", [])
            system_status.text = "Checks: " + (", ".join(faults) if faults else "healthy" if health.get("system_known") else "pending")
            system_status.classes(replace="text-amber-200" if faults else "text-slate-300")
            link_status.classes(replace="text-red-300" if active and stale else "text-slate-300")
            voltage, remaining = health.get("battery_v"), health.get("battery_pct")
            battery_status.text = "Battery (sim): " + ("unknown" if voltage is None else f"{voltage:.1f} V") + ("" if remaining is None else f" · {remaining}%")
            battery_status.set_visibility(voltage is not None)
            title, detail, tone, step_index = flight_story(view)
            if launch_pending[0]:
                title, detail, tone, step_index = "Preparing your flight", activity["detail"], "neutral", 1
            flight_title.text, flight_detail.text = title, detail
            readiness_panel.set_visibility(bool(view.get("run_id")) or active)
            flight_panel.classes(remove="fd-result-neutral fd-result-success fd-result-error fd-result-warning", add=f"fd-result-{tone}")
            flight_panel.set_visibility(bool(view.get("run_id")) or active or launch_pending[0])
            display_scope.text = "SAVED RECORDING" if historical else "LIVE SIMULATION" if active else "FLIGHT STATUS"
            for index, step in enumerate(workflow_steps):
                step.classes(replace="nicegui-row fd-workflow-step items-center gap-3 " + ("fd-step-current" if index == step_index else "fd-step-done" if index < step_index else ""))
                if index == step_index:
                    step.props('aria-current="step"')
                else:
                    step.props(remove="aria-current")
            status.text = ("SAVED RECORDING" if historical else "NEEDS ATTENTION" if view.get("error") else "STARTING" if activity else "FLYING" if active else "READY")
            recorded_setup = view.get("flight_setup") or {}
            target_time = recorded_setup.get("duration_s") or manifest.get("scenario_options", {}).get("duration_s")
            flight_time = view.get("free_time_s")
            timed = flight_time is not None and target_time is not None and target_time > 0
            flight_progress.set_visibility(timed)
            flight_duration.set_visibility(timed)
            if timed:
                flight_progress.value = min(1, max(0, flight_time/target_time))
                flight_duration.text = f"Flight time {flight_time:.1f} / {target_time:g} s" + (" · recorded" if historical else " · simulated time")
            overview_context.text = detail
            transition_run = manifest.get("test_kind") == "quad_transition"
            acceptance_note.text = ("Transition passes only after the full flight duration and at least 10 continuous seconds in firmware-confirmed FBWA at 20–30 m/s, roll and pitch within ±15°, altitude within ±15 m, and all lift motors at armed idle. Recorded heartbeats must remain fresh throughout that interval. Firmware handover alone is not a pass." if transition_run else
                                    "Each automatic test checks its recorded duration and scenario-specific limits. Hover checks include attitude, altitude and lift motor range. Interactive sessions record completion; they do not establish flight stability.")
            if is_back:
                acceptance_note.text = "Forward-to-hover requires full duration and 10 continuous seconds with fresh QLOITER confirmation, TAS ≤2 m/s, vertical speed ≤1 m/s, roll/pitch within ±10°, altitude within ±15 m, active lift motors and pusher off. Mode entry alone is not a pass."
            if manifest.get("doublet_plan"):
                acceptance_note.text = "Doublet checks verify both actual physical input lobes, neutral input restoration, recorded duration and attitude envelope. Returned-near-trim response is reported separately. Nonzero sideslip uses unmeasured β assumptions."
            acceptance_note.text += " These are simulation results from a provisional model; the physical aircraft and its CG still need validation."
            if transition_run and view.get("report") and "mode_evidence" not in view["report"]:
                acceptance_note.text += " This archived verdict uses its original checks, which confirmed FBWA at the end; stronger heartbeat reanalysis is saved separately."
            messages = view.get("notifications", [])
            visible = [m for m in messages if message_filter.value == "All messages" or m["severity"] <= 4]
            stamp_messages = (view.get("run_id"), message_filter.value, repr(visible))
            if last_messages[0] != stamp_messages:
                zero = flight_time_origin(view) or 0.
                message_table.rows = [{**m, "id": i, "time": f"T{m['time_s']-zero:+.1f} s"} for i, m in reversed(list(enumerate(visible)))]
                message_table.update()
                last_messages[0] = stamp_messages
            warning_count = sum(m['severity'] <= 4 for m in messages)
            warning_label = "warning" if warning_count == 1 else "warnings"
            messages_tab.props(f'label="Messages · {warning_count} recorded {warning_label}"' if warning_count and historical else f'label="Messages · {warning_count} {warning_label}"' if warning_count else 'label="Messages"')
            message_count.text = f"{len(visible)} shown · {len(messages)} recent messages" + (" · past events from this recording" if historical else " · latest 200 retained")
            row = view["truth"]
            values_status.text, values_stale = flight_values_status(view)
            values_panel.classes(add="fd-values-stale" if values_stale else "", remove="" if values_stale else "fd-values-stale")
            for key, value in numeric_flight_values(view).items():
                value_readouts[key].text = value
            model = (view.get("manifest") or {}).get("control_model")
            if model:
                travel = model["assumed_travel_deg"]
                travel_readout.text = (f"Assumed travel: elevator ±{travel['elevator']:g}° · aileron ±{travel['aileron']:g}° · rudder ±{travel['rudder']:g}° · "
                                       f"{'linear extrapolation' if model['deflection_policy'] == 'linear' else 'endpoint hold'} · uncalibrated. " + control_coverage_text(displayed_aircraft))
            elif row is not None:
                travel_readout.text = "Recorded travel setup: see manifest. " + control_coverage_text(displayed_aircraft)
            else:
                travel_readout.text = (f"Assumed surface travel ±{30 if prototype_setup else 20}° · linear extrapolation · uncalibrated. "
                                       + control_coverage_text(displayed_aircraft))
            attitude_empty.set_visibility(row is None)
            alerts = live_alerts(row)
            live_alert.text = ("" if not alerts else ("Recorded · " if historical else "") + " · ".join(alerts))
            live_alert.set_visibility(bool(alerts))
            if row is None:
                progress.text = activity["title"] if activity else "Select a flight test to begin."
                connection.text = "Waiting for firmware" if active else "Firmware not running"
                warning.text = ""
                warning.set_visibility(False)
                truth_readout.text = "Waiting for flight data"
                aero_readout.text = ""
                fdm_display.set_visibility(False)
                report.text = ""
                run_label.text = "" if view["run_id"] is None else f"Run: {view['run_id']}"
                return
            t = row["time_s"]
            free = view["free_time_s"]
            progress.text = f"Simulator clock {t:.1f} s · " + ("calibrating / arming on ground" if free is None and manifest.get("ground_contact") else
                                              "calibrating / arming" if free is None else
                                              f"{free:.1f} s since takeoff command" if manifest.get("ground_contact") else f"{free:.1f} s free flight")
            x = row["state"]
            fdm_display.set_visibility(True)
            roll, pitch = math.degrees(x[6]), math.degrees(x[7])
            attitude_stamp = (view["run_id"], t)
            if last_attitude[0] != attitude_stamp:
                fdm_display.content = attitude_html(pitch_deg=pitch, roll_deg=roll)
                fdm_display.update()
                last_attitude[0] = attitude_stamp
            tas = row.get("true_airspeed_m_s") or math.sqrt(sum(v*v for v in x[3:6]))
            truth_readout.text = f"Airspeed {tas:.2f} m/s · Altitude {-x[2]:.2f} m"
            aero_readout.text = f"α {row['alpha_deg']:.3f}° · β {row['beta_deg']:.3f}° · Elevator {math.degrees(row['control']['de']):+.2f}°"
            connection.text = (("LIVE · " if active else "RUN ENDED · Last sample: ") +
                               f"{'ARMED' if view['armed'] else 'DISARMED'} · "
                               f"Mode {view.get('mode_name') or view['mode']} · "
                               f"EKF flags {view['ekf_flags']} · Physics {row['physics_rate_hz']} Hz · "
                               f"{'BOOTSTRAP HOLD' if row['hold'] else 'FIRMWARE CONTROL'}")
            notices = []
            if active and view["stale_truth"]:
                notices.append("FDM connection stale")
            if active and view["stale_mavlink"]:
                notices.append("MAVLink connection stale")
            if row["control_coverage_exceeded"]:
                notices.append("Model uses estimates outside measured aerodynamic data")
            extrapolated = (row.get("control_coverage") or {}).get("extrapolated_surfaces", [])
            if extrapolated:
                notices.append("Control deflection extrapolated: " + ", ".join(extrapolated))
            warning.text = " · ".join(notices)
            if historical and notices:
                warning.text = "Recorded model warning · " + warning.text
            warning.set_visibility(bool(warning.text))
            run_label.text = f"Run: {view['run_id']} · Logs: {view['directory']}"
            displayed_directory[0] = view.get("directory")
            result = view["report"]
            if result is not None:
                rms = result.get("attitude_rms_error_deg")
                if result.get("test_kind") in ("quad_transition", "back_transition"):
                    report.text = transition_result_text(result)
                elif result.get("test_kind") == "control_doublet":
                    report.text = doublet_result_text(result)
                elif result.get("test_kind") == "quad_hover" and rms:
                    report.text = (f"{view['status']} · Nominal hover · Roll / pitch RMS {rms[0]:.3f}° / {rms[1]:.3f}° · "
                                   f"Altitude error ≤{result['settled_max_altitude_error_m']:.2f} m · "
                                   f"Horizontal drift {result['horizontal_displacement_m']:.1f} m (QHOVER)")
                elif result.get("test_kind") == "manual_elevator_demo":
                    report.text = (f"{view['status']} · MANUAL elevator input · "
                                   f"Pitch excursion {result['pitch_excursion_deg']:.2f}° · "
                                   f"Final pitch error {result['final_5s_pitch_error_deg']:.2f}°")
                else:
                    report.text = (f"{view['status']} · " + (f"FDM target tracking · Roll / pitch RMS {rms[0]:.4f}° / {rms[1]:.4f}°" if rms else result.get("error", "Run ended")))
            else:
                report.text = ""
            # Display updates only; no physics or MAVLink calls from a browser timer.
            stamp = (view["run_id"], int(t/2))
            if detail_tabs.value in (plots_tab, "Flight plots") and tuple(last_plot) != stamp:
                for name, figure in sitl_figures(manager.snapshot()).items():
                    plots[name].figure = figure
                    plots[name].update()
                plots["Doublet input"].set_visibility(bool(manifest.get("doublet_plan")))
                last_plot[:] = stamp

        aircraft.on_value_change(lambda: aircraft_changed())
        initial_point.on_value_change(lambda: point_changed())
        maneuver.on_value_change(lambda: maneuver_changed())
        altitude.on_value_change(lambda: conditions_changed())
        speed.on_value_change(lambda: conditions_changed())
        duration.on_value_change(lambda: refresh())
        cg_case.on_value_change(lambda: (refresh(), schedule_cg_card()))
        compare_button.on_click(lambda: compare_file(compare_path.value or ""))
        compare_upload.on_upload(compare_uploaded)
        ui.timer(.5, schedule_cg_card, once=True)
        transition_parameters.on_value_change(lambda: refresh())
        for field in (wind_speed, wind_from, gust, failure, failure_time):
            field.on_value_change(lambda: refresh())
        refresh()
        # The passive observer publishes at 10 Hz. Hidden plots are not rebuilt,
        # and an unchanged attitude sample is not sent to the browser again.
        detail_tabs.on_value_change(lambda: refresh())
        ui.timer(.1, refresh)
