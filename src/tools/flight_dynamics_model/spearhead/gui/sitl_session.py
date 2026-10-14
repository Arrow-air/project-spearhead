"""Passive log observer and owner of one local SITL runner; no physics clock."""

from collections import deque
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time
import uuid

import yaml

from ..sitl.backend import CONFIG_ROOT, verify_firmware
from ..sitl.ground_control import MODE_NAMES, session_modes
from .ground_status import StatusTexts, firmware_health
from .flight_setup import FlightSetup


SCENARIOS = {
    **{"quad_"+name: "Hover · "+name.replace("_", " ") for name in ("pitch_doublet", "roll_doublet", "yaw_doublet", "throttle_doublet")},
    "back_transition": "PT1.5 · forward trim to hover",
    "pitch_doublet": "Pitch doublet · physical elevator",
    "roll_doublet": "Roll doublet · physical aileron",
    "yaw_doublet": "Yaw doublet · physical rudder",
    "throttle_doublet": "Throttle doublet · physical pusher command",
    "quad_transition": "PT1.5 · hover to forward flight diagnostic",
    "quad_ground": "PT1.5 · ground takeoff, hover and landing",
    "quad_interactive": "PT1.5 QuadPlane · all modes · 10 min",
    "quad_hover": "PT1.5 QuadPlane hover · KDE · 60 s",
    "interactive": "Interactive flight · mode controls · 10 min",
    "fbwa": "Nominal FBWA · 60 s",
    "steps": "Pitch / roll steps · 118 s",
    "plus": "+0.5° roll / pitch disturbance · 60 s",
    "minus": "−0.5° roll / pitch disturbance · 60 s",
    "manual": "MANUAL direction checks · 60 s",
    "elevator": "Elevator ±1° doublet · MANUAL · 32 m/s · 45 s",
}


class JsonlTail:
    """Read complete appended lines, retaining a partial UTF-8 line for next poll."""

    def __init__(self, path):
        self.path, self.offset, self.pending = Path(path), 0, b""

    def read(self):
        if not self.path.exists():
            return []
        with self.path.open("rb") as stream:
            if os.fstat(stream.fileno()).st_size < self.offset:
                raise ValueError("Observed SITL log was truncated")
            stream.seek(self.offset)
            chunk = stream.read(2*1024*1024)
            self.offset += len(chunk)
        lines = (self.pending+chunk).split(b"\n")
        self.pending = lines.pop()
        if len(self.pending) > 2*1024*1024:
            raise ValueError("Oversized SITL log record")
        return [json.loads(line) for line in lines if line.strip()]


def body_attitude(message, pitch_trim):
    """Undo Plane's ATTITUDE pitch display trim (FLIGHT_OPTIONS=0)."""
    return {"roll_deg": math.degrees(message["roll"]),
            "pitch_deg": math.degrees(message["pitch"])+pitch_trim,
            "heading_deg": math.degrees(message["yaw"])}


class RunObserver:
    """Incremental, bounded display history from the runner's authoritative logs."""

    def __init__(self, directory, history_limit=3000):
        # About 150 seconds at 20 Hz keeps preset inputs visible after completion.
        self.directory = Path(directory)
        self.truth_tail = JsonlTail(self.directory / "truth.jsonl")
        self.telemetry_tail = JsonlTail(self.directory / "telemetry.jsonl")
        self.control_tail = JsonlTail(self.directory / "control-events.jsonl")
        self.notifications = deque(maxlen=200)
        self.status_texts = StatusTexts(self.notifications)
        self.last_command = None
        self.last_heartbeat_wall = None
        self.history = deque(maxlen=history_limit)
        self.estimates = deque(maxlen=history_limit)
        self.pilot_inputs = deque(maxlen=history_limit)
        self.latest = None
        self.messages = {}
        self.manifest = None
        self.report = None
        self.last_truth_wall = None
        self.free_start = None
        self.rows_seen = 0
        self.flight_setup = None
        self.recipe_loaded = False
        self.phase = None
        self.started_wall = time.monotonic()
        self.transient_metrics = {}

    def poll(self):
        phase_path = self.directory / "phase.json"
        if phase_path.exists():
            self.phase = json.loads(phase_path.read_text())
        manifest = self.directory / "manifest.json"
        if self.manifest is None and manifest.exists():
            try:
                self.manifest = json.loads(manifest.read_text())
            except json.JSONDecodeError:
                pass  # Writer may be replacing the startup manifest.
        recipe = self.directory.parent / (self.directory.name+".session.json")
        if not self.recipe_loaded and recipe.exists():
            self.recipe_loaded = True
            try:
                saved = json.loads(recipe.read_text()).get("flight_setup")
                if saved is not None:
                    if saved.get("maneuver") == "transition" and "transition_parameters" not in saved:
                        # Older GUI transitions used the PT1.5 overlay; do not
                        # relabel those recordings with the new proposal default.
                        saved["transition_parameters"] = (self.manifest or {}).get("transition_plan", {}).get("parameter_set", "pt15")
                    if saved.get("initial_point") == "hover30":
                        saved["initial_point"] = "custom_hover"
                    self.flight_setup = FlightSetup(**saved).validate().as_dict()
            except (OSError, TypeError, ValueError):
                pass  # Optional labels must not prevent observing/stopping a run.
        for row in self.truth_tail.read():
            self.latest = row
            self.rows_seen += 1
            self.last_truth_wall = time.monotonic()
            if not row["hold"] and self.free_start is None and not (self.manifest or {}).get("ground_contact"):
                self.free_start = row["time_s"]
            if not row["hold"] and (self.manifest or {}).get("test_kind") in ("quad_transition", "back_transition"):
                pitch = abs(math.degrees(row["state"][7]))
                initial_altitude = -self.manifest["initial_state"][2]
                self.transient_metrics["peak_abs_pitch_deg"] = max(
                    self.transient_metrics.get("peak_abs_pitch_deg", 0), pitch)
                self.transient_metrics["altitude_loss_m"] = max(
                    self.transient_metrics.get("altitude_loss_m", 0), initial_altitude+row["state"][2])
                if row["time_s"]-self.free_start < 10:
                    key = "initial_forward_peak_pitch_deg" if self.manifest.get("test_kind") == "back_transition" else "hover_peak_abs_pitch_deg"
                    self.transient_metrics[key] = max(self.transient_metrics.get(key, 0), pitch)
            if not self.history or row["time_s"]-self.history[-1]["time_s"] >= .049:
                # Display history needs only these fields. Full sensor/coverage
                # records remain in truth.jsonl and the current-state readout.
                self.history.append({key: row[key] for key in ("time_s", "hold", "state", "pwm", "control")})
        if self.phase is not None:
            self.free_start = self.phase.get("flight_start_s")
        for row in self.telemetry_tail.read():
            message = row["message"]
            kind = message["mavpackettype"]
            self.messages[kind] = row
            if kind == "RC_CHANNELS" and (self.manifest or {}).get("doublet_plan", {}).get("regime") == "hover":
                self.pilot_inputs.append({"time_s": row["fdm_time_s"], "channels": [message.get(f"chan{i}_raw") for i in range(1, 5)]})
            if kind == "HEARTBEAT":
                self.last_heartbeat_wall = time.monotonic()
            if kind == "STATUSTEXT":
                self.status_texts.add(message, row["fdm_time_s"])
            if kind == "ATTITUDE" and self.manifest is not None:
                trim = self.manifest["parameter_overlay"]["PTCH_TRIM_DEG"]
                self.estimates.append({"time_s": row["fdm_time_s"], **body_attitude(message, trim)})
        for row in self.control_tail.read():
            self.last_command = row
            self.notifications.append({"time_s": row["time_s"], "source": "Flight Deck",
                                       "severity": 4 if row["status"] in ("rejected", "timeout") else 6,
                                       "level": "Warning" if row["status"] in ("rejected", "timeout") else "Info",
                                       "text": row["text"]})
        report = self.directory / "report.json"
        if report.exists():
            try:
                self.report = json.loads(report.read_text())
            except json.JSONDecodeError:
                pass

    def snapshot(self, *, active=False, stopped=False, error=None):
        now = time.monotonic()
        truth = self.latest
        t = None if truth is None else truth["time_s"]
        heartbeat = self.messages.get("HEARTBEAT")
        attitude = self.messages.get("ATTITUDE")
        stale_truth = active and (self.last_truth_wall is None or now-self.last_truth_wall > 2)
        stale_mavlink = active and (heartbeat is None or t is None or t-heartbeat["fdm_time_s"] > 2.5
                                   or self.last_heartbeat_wall is None or now-self.last_heartbeat_wall > 2.5)
        estimate = self.estimates[-1] if self.estimates else None
        if active and (attitude is None or t is None or t-attitude["fdm_time_s"] > .75):
            estimate = None
        status = "STARTING" if active and truth is None else "READY"
        if active and truth is not None:
            status = "BOOTSTRAP HOLD" if truth["hold"] else "FREE FLIGHT"
            if self.phase:
                status = self.phase["phase"]
        if not active and truth is not None:
            status = "PASSED" if self.report and self.report.get("passed") else "FAILED"
            if self.report and self.report.get("test_kind") == "interactive_session" and self.report.get("passed"):
                status = "COMPLETED"
        if stopped:
            status = "STOPPED"
        if error:
            status = "ERROR"
        ekf = self.messages.get("EKF_STATUS_REPORT", {}).get("message", {})
        return deepcopy({"status": status, "active": active, "run_id": self.directory.name,
                         "directory": str(self.directory), "truth": truth,
                         "history": list(self.history), "pilot_inputs": list(self.pilot_inputs), "estimates": list(self.estimates),
                         "estimate": estimate, "manifest": self.manifest, "report": self.report,
                         "transient_metrics": self.transient_metrics,
                         "time_s": t, "free_time_s": None if self.free_start is None else t-self.free_start,
                         "wall_elapsed_s": now-self.started_wall if active else None,
                         "truth_age_s": None if self.last_truth_wall is None else now-self.last_truth_wall,
                         "stale_truth": stale_truth, "stale_mavlink": stale_mavlink,
                         "armed": None if heartbeat is None else bool(heartbeat["message"]["base_mode"] & 128),
                         "mode": None if heartbeat is None else heartbeat["message"]["custom_mode"],
                         "mode_name": None if heartbeat is None else MODE_NAMES.get(heartbeat["message"]["custom_mode"], f"Mode {heartbeat['message']['custom_mode']}"),
                         "interactive": bool((self.manifest or {}).get("interactive")),
                         "mode_controls_ready": (not (self.manifest or {}).get("requires_controls_ready")
                                                 or (self.directory / "controls-ready.json").exists()),
                         "notifications": sorted(self.notifications, key=lambda r: r["time_s"]),
                         "flight_setup": self.flight_setup,
                         "command": self.last_command,
                         "health": firmware_health({} if stale_mavlink else self.messages, t, live=active),
                         "ekf_flags": ekf.get("flags"), "error": error, "rows_seen": self.rows_seen})


class SITLManager:
    """Shared server-owned runner; browser connections only read cached snapshots."""

    def __init__(self, *, binary=None, root=None, python=None, port=19602, mavlink_port=16360):
        self.binary = Path(binary or CONFIG_ROOT.parent / ".sitl/arduplane-4.7.1").resolve()
        self.root = Path(root or CONFIG_ROOT.parent / ".sitl/flight-deck").resolve()
        self.python = str(python or sys.executable)
        self.port, self.mavlink_port = port, mavlink_port
        self.process = self.worker = self.observer = None
        self.lock = threading.RLock()
        self.stopped = False
        self.error = None
        self.pending_command = None
        self.cached = {"status": "READY", "active": False, "run_id": None,
                       "truth": None, "history": [], "estimates": [], "estimate": None,
                       "manifest": None, "report": None, "error": None,
                       "interactive": False, "notifications": [], "command": None, "health": {}}

    def start_configured(self, setup):
        if not isinstance(setup, FlightSetup):
            raise ValueError("Choose a valid flight setup")
        setup.validate()
        return self.start(setup.scenario, setup=setup)

    def start(self, scenario, *, setup=None):
        if scenario not in SCENARIOS:
            raise ValueError("Choose a supported SITL scenario")
        if setup is not None:
            setup.validate()
            if scenario != setup.scenario:
                raise ValueError("Flight setup and runner scenario must agree")
        if os.name != "posix":
            raise RuntimeError("Managed SITL processes currently require macOS or Linux")
        with self.lock:
            if ((self.process is not None and self.process.poll() is None)
                    or (self.worker is not None and self.worker.is_alive())):
                raise RuntimeError("A SITL run is already active")
            verify_firmware(self.binary)
            self.root.mkdir(parents=True, exist_ok=True)
            directory = self.root / (time.strftime("%Y%m%d-%H%M%S")+"-"+uuid.uuid4().hex[:8])
            command = [self.python, "-m", "spearhead.sitl.flight",
                       scenario if scenario in ("manual", "elevator", "steps", "back_transition", "pitch_doublet", "roll_doublet", "yaw_doublet", "throttle_doublet") else "fbwa",
                       "--binary", str(self.binary), "--output", str(directory),
                       "--port", str(self.port), "--mavlink-port", str(self.mavlink_port)]
            if scenario in ("quad_hover", "quad_interactive", "quad_transition") or scenario.startswith("quad_") and scenario.endswith("_doublet"):
                command = [self.python, "-m", "spearhead.sitl.quad", "--binary", str(self.binary),
                           "--output", str(directory), "--port", str(self.port),
                           "--mavlink-port", str(self.mavlink_port)]
                if scenario.endswith("_doublet"):
                    command += ["--doublet", scenario.removeprefix("quad_")]
                if scenario == "quad_transition":
                    command += ["--transition", "--transition-parameters",
                                setup.transition_parameters if setup is not None else "first_transition_2026_10_09"]
            if scenario == "quad_ground":
                command = [self.python, "-m", "spearhead.sitl.ground_flight", "--binary", str(self.binary),
                           "--output", str(directory), "--port", str(self.port),
                           "--mavlink-port", str(self.mavlink_port)]
            if scenario in ("plus", "minus"):
                delta = ".5" if scenario == "plus" else "-.5"
                command += ["--perturb-roll", delta, "--perturb-pitch", delta]
            if scenario in ("interactive", "quad_interactive"):
                command += ["--interactive", "--duration", "600"]
            if setup is not None:
                config_path = self.root / (directory.name+".setup.yaml")
                config_path.write_text(yaml.safe_dump(setup.config(), sort_keys=False))
                command += ["--config", str(config_path)]
                if setup.aircraft == "pt15" and setup.is_forward:
                    command += ["--forward-parameters", setup.transition_parameters]
                if "--duration" in command:
                    command[command.index("--duration")+1] = str(setup.duration_s)
                else:
                    command += ["--duration", str(setup.duration_s)]
            env = os.environ.copy()
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env["OPENBLAS_NUM_THREADS"] = "1"
            env["PYTHONPATH"] = str(CONFIG_ROOT.parent)+os.pathsep+env.get("PYTHONPATH", "")
            self.observer = RunObserver(directory)
            self.observer.flight_setup = None if setup is None else setup.as_dict()
            self.stopped, self.error = False, None
            self.pending_command = None
            (self.root / (directory.name+".session.json")).write_text(json.dumps(
                {"scenario": scenario, "cancelled": False,
                 "flight_setup": None if setup is None else setup.as_dict()}))
            with (self.root / (directory.name+".console.log")).open("xb") as log:
                self.process = subprocess.Popen(command, cwd=CONFIG_ROOT.parent, env=env,
                                                stdout=log, stderr=subprocess.STDOUT,
                                                start_new_session=True)
            self.cached = self.observer.snapshot(active=True)
            self.worker = threading.Thread(target=self._watch, daemon=True)
            self.worker.start()
            return directory

    def _watch(self):
        try:
            while True:
                with self.lock:
                    self.observer.poll()
                    active = self.process.poll() is None
                    if not active:
                        self.observer.poll()  # Include buffers flushed at runner shutdown.
                    if not active and self.process.returncode != 0 and not self.stopped:
                        try:
                            os.killpg(self.process.pid, signal.SIGTERM)
                        except ProcessLookupError:
                            pass
                        self.error = (self._console_error() if self.observer.report is None
                                      else self.observer.report.get("error"))
                    self.cached = self.observer.snapshot(active=active, stopped=self.stopped, error=self.error)
                    if self.pending_command is not None:
                        event = self.cached.get("command")
                        if event and event.get("command_id") == self.pending_command["command_id"]:
                            self.pending_command = event
                        self.cached["command"] = deepcopy(self.pending_command)
                if not active:
                    return
                time.sleep(.1)
        except Exception as exc:
            with self.lock:
                self.error = str(exc)
                self.cached = self.observer.snapshot(active=True, error=self.error)
            # Observer errors must not abandon a child firmware process.
            self.stop()
            with self.lock:
                self.cached = self.observer.snapshot(active=False, stopped=True, error=self.error)

    def _console_error(self):
        path = self.root / (self.observer.directory.name+".console.log")
        with path.open("rb") as stream:
            stream.seek(max(0, path.stat().st_size-1600))
            return stream.read().decode(errors="replace") or "SITL runner exited without a report"

    def snapshot(self, *, include_history=True):
        with self.lock:
            if include_history:
                return deepcopy(self.cached)
            return deepcopy({**self.cached, "history": self.cached["history"][-600:], "estimates": []})

    def set_mode(self, mode):
        """Queue only to this manager's current interactive simulation runner."""
        with self.lock:
            view = self.cached
            if mode not in session_modes(view.get("manifest")):
                raise ValueError("Choose a supported mode for this aircraft configuration")
            if (self.process is None or self.process.poll() is not None or not view.get("interactive")
                    or view.get("status") != "FREE FLIGHT" or view.get("stale_truth") or view.get("stale_mavlink")
                    or view.get("error") or not view.get("mode_controls_ready", True)):
                raise RuntimeError("Mode controls require a live interactive flight after calibration and arming")
            if self.pending_command and self.pending_command["status"] in ("queued", "sent", "accepted"):
                raise RuntimeError("Wait for the current mode change to complete")
            ident = uuid.uuid4().hex
            directory = self.observer.directory / "commands"
            directory.mkdir(exist_ok=True)
            request = {"id": ident, "kind": "set_mode", "mode": mode, "created_at": time.monotonic()}
            path = directory / (ident+".tmp")
            path.write_text(json.dumps(request))
            path.replace(directory / (ident+".json"))
            self.pending_command = {"command_id": ident, "mode": mode, "status": "queued",
                                    "text": f"Queued {MODE_NAMES[mode]}"}
            self.cached["command"] = deepcopy(self.pending_command)
            return ident

    def _recording_paths(self):
        roots = {"": self.root}
        # CLI diagnostics are siblings of the default managed-session root.
        # Explicit test/custom roots stay isolated from repository recordings.
        if self.root == (CONFIG_ROOT.parent / ".sitl" / "flight-deck").resolve():
            roots["diagnostic/"] = self.root.parent
        return {prefix + path.parent.name: path.parent for prefix, root in roots.items()
                for path in sorted(root.glob("*/report.json"), reverse=True)
                if path.parent.resolve().parent == root.resolve()}

    def completed_runs(self):
        return {key: key for key in self._recording_paths()}

    def load_completed(self, run_id):
        """Passively view a finished managed session or local CLI diagnostic."""
        with self.lock:
            if self.process is not None and self.process.poll() is None:
                raise RuntimeError("Stop the active run before opening a recording")
            directory = self._recording_paths().get(run_id)
            if directory is None:
                raise ValueError("Unknown recorded SITL run")
            observer = RunObserver(directory)
            while True:
                observer.poll()
                if all(not tail.path.exists() or tail.offset == tail.path.stat().st_size
                       for tail in [observer.truth_tail, observer.telemetry_tail]):
                    break
            manifest = observer.manifest or {}
            if (manifest.get("config", {}).get("firmware") != "Plane-4.7.1"
                    or manifest.get("config", {}).get("aerodynamics") != "tail_update_2026_10_07"):
                raise ValueError("Recording requires 4.7.1 and the October AeroDB")
            self.observer = observer
            self.pending_command = None
            self.stopped, self.error = False, None
            session_path = directory.parent / (directory.name+".session.json")
            cancelled = session_path.exists() and json.loads(session_path.read_text())["cancelled"]
            self.cached = observer.snapshot(stopped=cancelled)

    def stop(self):
        with self.lock:
            process = self.process
            if process is None or process.poll() is not None:
                return
            self.stopped = True
            path = self.root / (self.observer.directory.name+".session.json")
            try:
                session = json.loads(path.read_text())
                session["cancelled"] = True
                path.write_text(json.dumps(session))
            except (OSError, ValueError):
                pass  # A damaged/full log filesystem must not prevent stop.
        # Signal only our new process group, including the firmware child.
        for sig, timeout in [(signal.SIGINT, 5), (signal.SIGTERM, 2), (signal.SIGKILL, 2)]:
            try:
                os.killpg(process.pid, sig)
            except ProcessLookupError:
                break
            try:
                process.wait(timeout=timeout)
                break
            except subprocess.TimeoutExpired:
                continue
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        if self.worker is not None and self.worker is not threading.current_thread():
            self.worker.join(timeout=3)

    close = stop
