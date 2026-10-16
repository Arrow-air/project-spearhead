"""Compare a real DataFlash log, or a parameter file, with a recorded SITL run.

Alignment: the aircraft log starts at its first takeoff (all four lift outputs above
1250 µs and altitude above 1 m); the SITL run starts at release (air starts) or at
the first lift spin-up (ground starts). Traces are resampled to 10 Hz.
"""

import json
import math
from pathlib import Path
import re

import numpy as np

LOG_TYPES = ("PARM", "RCOU", "ATT", "BAT", "QTUN", "MODE")
# Lift outputs: ArduPilot motors 1-4 on SERVO5-8 = FR, RL, FL, RR.
LIFT_CHANNELS = {"FR": 5, "RL": 6, "FL": 7, "RR": 8}
TAKEOFF_PWM_US, TAKEOFF_ALT_M = 1250, 1.
SAMPLE_HZ = 10.
# Hardware, calibration, telemetry and simulator settings that legitimately differ.
IGNORED_PREFIXES = ("SIM_", "STAT_", "SYSID_", "INS_ACC", "INS_GYR", "INS_TCAL", "COMPASS_", "BARO", "GND_",
                    "AHRS_TRIM", "SR0_", "SR1_", "SR2_", "SR3_", "SR4_", "SR5_", "SR6_", "SERIAL", "CAN_",
                    "BRD_", "GPS", "LOG_", "FORMAT_VERSION", "STAT", "RC_", "NTF_", "SCR_", "MIS_", "FENCE_",
                    "RALLY_", "TERRAIN_", "OSD", "MSP_", "DID_", "RPM", "EAHRS_", "VISO", "AFS_", "ADSB_",
                    "MAV", "INS_LOG", "INS_FAST", "INS_STILL", "EK2_", "GEN_", "KDE_", "FLTMODE", "INITIAL_MODE",
                    "ARMING_", "TELEM", "BATT2_", "BATT3_", "LGR_", "RNGFND", "PRX", "WINCH", "MNT", "CAM",
                    "SID_", "NET_", "PUP_", "RELAY", "TEMP", "SERVO_ROB", "SERVO_VOLZ", "RSSI", "Q_AUTOTUNE")
# Pilot RC channels 5-16 carry switches, not stick ranges that change the flight.
IGNORED_PATTERN = re.compile(r"RC(?:[5-9]|1[0-6])_")


def read_dataflash(path, types=LOG_TYPES):
    """Return (parameters, {type: [row dicts]}) from a DataFlash .BIN file."""
    from pymavlink import DFReader
    log = DFReader.DFReader_binary(str(path), zero_time_base=False)
    parameters, series = {}, {kind: [] for kind in types if kind != "PARM"}
    while True:
        message = log.recv_match(type=list(types))
        if message is None:
            break
        kind = message.get_type()
        if kind == "PARM":
            parameters[message.Name] = float(message.Value)
        else:
            row = message.to_dict()
            row.pop("mavpackettype", None)
            series[kind].append(row)
    return parameters, series


def read_parameter_file(path):
    """Read a Mission Planner / QGC .param or .parm file (last value wins)."""
    parameters = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        fields = re.split(r"[\s,]+", line)
        # QGC format: vehicle component name value type.
        if len(fields) >= 5 and fields[0].isdigit() and fields[1].isdigit():
            fields = fields[2:4]
        try:
            parameters[fields[0]] = float(fields[1])
        except (IndexError, ValueError):
            continue
    return parameters


def run_parameters(directory):
    """Best available SITL parameter set for a run: its firmware log, else parm files plus overlay."""
    directory = Path(directory)
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    # SITL writes its own DataFlash logs under logs/; the last PARM values are what flew.
    for log_path in sorted((directory / "logs").glob("*.BIN"), key=lambda p: p.stat().st_mtime, reverse=True):
        parameters, _ = read_dataflash(log_path, types=("PARM",))
        if parameters:
            return parameters, "SITL firmware log"
    parameters = {}
    for parm in sorted(directory.glob("*.parm")):
        parameters.update(read_parameter_file(parm))
    for key, value in (manifest.get("parameter_overlay") or {}).items():
        if isinstance(value, (int, float)):
            parameters[key] = float(value)
    return parameters, "parameter files and overlay"


def parameter_diff(aircraft, simulation, *, include_ignored=False, tolerance=1e-4):
    """Rows for parameters set differently on the aircraft and in SITL, plus SITL-only keys."""
    rows = []
    for name in sorted(set(aircraft) | set(simulation)):
        if not include_ignored and (name.startswith(IGNORED_PREFIXES) or IGNORED_PATTERN.match(name)):
            continue
        a, s = aircraft.get(name), simulation.get(name)
        if a is None and s is None:
            continue
        if a is None:
            rows.append({"name": name, "aircraft": "—", "simulation": f"{s:g}", "status": "SITL only"})
        elif s is None:
            continue    # thousands of aircraft-only defaults; not useful in a diff
        elif abs(a - s) > tolerance*max(1., abs(a), abs(s)):
            rows.append({"name": name, "aircraft": f"{a:g}", "simulation": f"{s:g}", "status": "different"})
    return rows


def _times(rows):
    return np.array([row["TimeUS"] for row in rows], dtype=float)/1e6


def log_takeoffs(series, gap_s=5.):
    """(takeoff time, airborne duration) for each flight in the log, in log seconds."""
    rcou, qtun = series.get("RCOU") or [], series.get("QTUN") or []
    if not rcou or not qtun:
        return []
    t = _times(rcou)
    lift = np.array([[row.get(f"C{c}", 0) for c in LIFT_CHANNELS.values()] for row in rcou], dtype=float)
    altitude = np.interp(t, _times(qtun), [row.get("Alt", 0.) for row in qtun])
    airborne = t[(lift.min(axis=1) > TAKEOFF_PWM_US) & (altitude > TAKEOFF_ALT_M)]
    if not airborne.size:
        return []
    breaks = np.where(np.diff(airborne) > gap_s)[0]
    starts = np.r_[airborne[0], airborne[breaks + 1]]
    ends = np.r_[airborne[breaks], airborne[-1]]
    return [(float(s), float(e - s)) for s, e in zip(starts, ends)]


def log_takeoff_time(series, index=-1):
    flights = log_takeoffs(series)
    return flights[index][0] if flights and -len(flights) <= index < len(flights) else None


def log_traces(series, start_s, duration_s):
    """10 Hz traces from takeoff: lift PWM by position, roll/pitch (deg), altitude (m), pack voltage."""
    grid = np.arange(0., duration_s, 1/SAMPLE_HZ)
    traces = {"time_s": grid}
    rcou = series.get("RCOU") or []
    if rcou:
        t = _times(rcou) - start_s
        for name, channel in LIFT_CHANNELS.items():
            traces[f"lift_{name}"] = np.interp(grid, t, [row.get(f"C{channel}", np.nan) for row in rcou])
    for kind, key, out in (("ATT", "Roll", "roll_deg"), ("ATT", "Pitch", "pitch_deg"),
                           ("QTUN", "Alt", "altitude_m"), ("BAT", "Volt", "battery_v")):
        rows = series.get(kind) or []
        if rows:
            traces[out] = np.interp(grid, _times(rows) - start_s, [row.get(key, np.nan) for row in rows])
    return traces


def run_traces(directory, duration_s=None):
    """10 Hz traces from a SITL run's truth.jsonl, starting at release or ground spin-up."""
    rows = []
    with (Path(directory) / "truth.jsonl").open() as stream:
        for line in stream:
            row = json.loads(line)
            if not row.get("hold"):
                rows.append(row)
    if not rows:
        return None
    pwm = np.array([row["pwm"][:8] for row in rows], dtype=float)
    t = np.array([row["time_s"] for row in rows])
    spinning = pwm[:, 4:8].min(axis=1) > TAKEOFF_PWM_US
    start = t[np.argmax(spinning)] if spinning.any() else t[0]
    end = t[-1] - start if duration_s is None else duration_s
    grid = np.arange(0., max(end, 0.), 1/SAMPLE_HZ)
    state = np.array([row["state"] for row in rows])
    traces = {"time_s": grid}
    for name, channel in LIFT_CHANNELS.items():
        traces[f"lift_{name}"] = np.interp(grid, t - start, pwm[:, channel-1])
    traces["roll_deg"] = np.interp(grid, t - start, np.degrees(state[:, 6]))
    traces["pitch_deg"] = np.interp(grid, t - start, np.degrees(state[:, 7]))
    traces["altitude_m"] = np.interp(grid, t - start, -state[:, 2])
    voltage = [((row.get("sensor") or {}).get("battery") or {}).get("voltage", np.nan) for row in rows]
    traces["battery_v"] = np.interp(grid, t - start, np.asarray(voltage, dtype=float))
    return traces


def hover_summary(traces, settle_s=5.):
    """Median lift PWM per rotor and the front/rear PWM split after the climb settles."""
    if traces is None or "lift_FR" not in traces:
        return None
    mask = traces["time_s"] >= settle_s
    if not mask.any():
        return None
    medians = {name: float(np.nanmedian(traces[f"lift_{name}"][mask])) for name in LIFT_CHANNELS}
    above_idle = {name: value - 1100 for name, value in medians.items()}
    total = sum(above_idle.values())
    front = (above_idle["FR"] + above_idle["FL"])/total if total > 0 else math.nan
    return {"median_pwm": medians, "front_share_pwm": front,
            "roll_deg_median": float(np.nanmedian(traces["roll_deg"][mask])),
            "pitch_deg_median": float(np.nanmedian(traces["pitch_deg"][mask]))}


def compare_log_with_run(log_path, run_directory, flight=-1):
    """Everything the Review tab shows; parameter files give the diff only.

    ``flight`` picks the takeoff within the log (default: the last one).
    """
    log_path = Path(log_path)
    simulation, source = run_parameters(run_directory)
    if log_path.suffix.lower() in (".param", ".parm", ".params", ".txt"):
        aircraft, series = read_parameter_file(log_path), {}
    else:
        aircraft, series = read_dataflash(log_path)
    result = {"log": log_path.name, "simulation_parameter_source": source,
              "parameter_diff": parameter_diff(aircraft, simulation), "aircraft_parameter_count": len(aircraft)}
    sim = run_traces(run_directory)
    flights = log_takeoffs(series) if series else []
    result["flights"] = [{"takeoff_s": s, "airborne_s": d} for s, d in flights]
    takeoff = log_takeoff_time(series, flight) if series else None
    result["takeoff_s"] = takeoff
    if sim is not None and takeoff is not None:
        real = log_traces(series, takeoff, float(sim["time_s"][-1]) + 1/SAMPLE_HZ)
        result["traces"] = {"simulation": {k: v.tolist() for k, v in sim.items()},
                            "aircraft": {k: v.tolist() for k, v in real.items()}}
        result["hover"] = {"simulation": hover_summary(sim), "aircraft": hover_summary(real)}
    return result
