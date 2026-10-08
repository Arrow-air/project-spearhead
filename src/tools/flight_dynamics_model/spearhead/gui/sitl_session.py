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

from ..sitl.backend import CONFIG_ROOT, verify_firmware


SCENARIOS = {
    "fbwa": "Nominal FBWA · 60 s",
    "steps": "Pitch / roll steps · 118 s",
    "plus": "+0.5° roll / pitch disturbance · 60 s",
    "minus": "−0.5° roll / pitch disturbance · 60 s",
    "manual": "MANUAL direction checks · 60 s",
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

    def __init__(self, directory, history_limit=600):
        self.directory = Path(directory)
        self.truth_tail = JsonlTail(self.directory / "truth.jsonl")
        self.telemetry_tail = JsonlTail(self.directory / "telemetry.jsonl")
        self.history = deque(maxlen=history_limit)
        self.estimates = deque(maxlen=history_limit)
        self.latest = None
        self.messages = {}
        self.manifest = None
        self.report = None
        self.last_truth_wall = None
        self.free_start = None
        self.rows_seen = 0

    def poll(self):
        manifest = self.directory / "manifest.json"
        if self.manifest is None and manifest.exists():
            try:
                self.manifest = json.loads(manifest.read_text())
            except json.JSONDecodeError:
                pass  # Writer may be replacing the startup manifest.
        for row in self.truth_tail.read():
            self.latest = row
            self.rows_seen += 1
            self.last_truth_wall = time.monotonic()
            if not row["hold"] and self.free_start is None:
                self.free_start = row["time_s"]
            if not self.history or row["time_s"]-self.history[-1]["time_s"] >= .049:
                self.history.append(row)
        for row in self.telemetry_tail.read():
            message = row["message"]
            kind = message["mavpackettype"]
            self.messages[kind] = row
            if kind == "ATTITUDE" and self.manifest is not None:
                trim = self.manifest["parameter_overlay"]["PTCH_TRIM_DEG"]
                self.estimates.append({"time_s": row["fdm_time_s"], **body_attitude(message, trim)})
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
        stale_mavlink = active and (heartbeat is None or t is None or t-heartbeat["fdm_time_s"] > 2.5)
        estimate = self.estimates[-1] if self.estimates else None
        if active and (attitude is None or t is None or t-attitude["fdm_time_s"] > .75):
            estimate = None
        status = "STARTING" if active and truth is None else "READY"
        if active and truth is not None:
            status = "BOOTSTRAP HOLD" if truth["hold"] else "FREE FLIGHT"
        if not active and truth is not None:
            status = "PASSED" if self.report and self.report.get("passed") else "FAILED"
        if stopped:
            status = "STOPPED"
        if error:
            status = "ERROR"
        ekf = self.messages.get("EKF_STATUS_REPORT", {}).get("message", {})
        return deepcopy({"status": status, "active": active, "run_id": self.directory.name,
                         "directory": str(self.directory), "truth": truth,
                         "history": list(self.history), "estimates": list(self.estimates),
                         "estimate": estimate, "manifest": self.manifest, "report": self.report,
                         "time_s": t, "free_time_s": None if self.free_start is None else t-self.free_start,
                         "stale_truth": stale_truth, "stale_mavlink": stale_mavlink,
                         "armed": None if heartbeat is None else bool(heartbeat["message"]["base_mode"] & 128),
                         "mode": None if heartbeat is None else heartbeat["message"]["custom_mode"],
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
        self.cached = {"status": "READY", "active": False, "run_id": None,
                       "truth": None, "history": [], "estimates": [], "estimate": None,
                       "manifest": None, "report": None, "error": None}

    def start(self, scenario):
        if scenario not in SCENARIOS:
            raise ValueError("Choose a supported SITL scenario")
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
                       scenario if scenario in ("manual", "steps") else "fbwa",
                       "--binary", str(self.binary), "--output", str(directory),
                       "--port", str(self.port), "--mavlink-port", str(self.mavlink_port)]
            if scenario in ("plus", "minus"):
                delta = ".5" if scenario == "plus" else "-.5"
                command += ["--perturb-roll", delta, "--perturb-pitch", delta]
            env = os.environ.copy()
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env["OPENBLAS_NUM_THREADS"] = "1"
            env["PYTHONPATH"] = str(CONFIG_ROOT.parent)+os.pathsep+env.get("PYTHONPATH", "")
            self.observer = RunObserver(directory)
            self.stopped, self.error = False, None
            (self.root / (directory.name+".session.json")).write_text(json.dumps(
                {"scenario": scenario, "cancelled": False}))
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

    def snapshot(self):
        with self.lock:
            return deepcopy(self.cached)

    def completed_runs(self):
        return {path.parent.name: path.parent.name for path in
                sorted(self.root.glob("*/report.json"), reverse=True)}

    def load_completed(self, run_id):
        """View only a finished run in this manager's own log directory."""
        with self.lock:
            if self.process is not None and self.process.poll() is None:
                raise RuntimeError("Stop the active run before opening a recording")
            if (self.root / run_id).resolve().parent != self.root or run_id not in self.completed_runs():
                raise ValueError("Unknown recorded SITL run")
            observer = RunObserver(self.root / run_id)
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
            self.stopped, self.error = False, None
            session_path = self.root / (run_id+".session.json")
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
