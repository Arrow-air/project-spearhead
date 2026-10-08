"""SITL display isolation, incremental logs and owned process lifecycle."""

import json
import math
from pathlib import Path

import pytest

from spearhead.gui import sitl_session
from spearhead.gui.sitl_page import sitl_figures
from spearhead.gui.sitl_session import JsonlTail, RunObserver, SITLManager, body_attitude


def truth(t=1, hold=False):
    return {"time_s": t, "state": [0, 0, -100, 28, 0, 0, .01, -.02, 0, 0, 0, 0],
            "hold": hold, "pwm": [1500, 1801, 1254, 1495], "physics_rate_hz": 200,
            "alpha_deg": -.1, "beta_deg": 0, "control_coverage_exceeded": False}


def append(path, row):
    with path.open("a") as stream:
        stream.write(json.dumps(row)+"\n")


def test_tail_preserves_partial_utf8_records_and_rejects_truncation(tmp_path):
    path = tmp_path / "stream.jsonl"
    tail = JsonlTail(path)
    assert tail.read() == []
    encoded = '{"text":"β"}\n'.encode()
    split = encoded.index("β".encode())+1
    path.write_bytes(encoded[:split])
    assert tail.read() == []
    with path.open("ab") as stream:
        stream.write(encoded[split:])
    assert tail.read() == [{"text": "β"}]
    assert tail.read() == []
    path.write_bytes(b"")
    with pytest.raises(ValueError, match="truncated"):
        tail.read()


def test_observer_is_passive_bounded_and_snapshot_cannot_mutate_backend(tmp_path):
    for i in range(20):
        append(tmp_path / "truth.jsonl", truth(i*.05, hold=i < 5))
    observer = RunObserver(tmp_path, history_limit=3)
    observer.poll()
    before = (tmp_path / "truth.jsonl").read_bytes()
    one = observer.snapshot(active=True)
    assert len(one["history"]) == 3 and one["rows_seen"] == 20
    assert one["free_time_s"] == pytest.approx(.7)
    one["truth"]["state"][0] = 999
    for _ in range(5):
        observer.snapshot()
        observer.poll()
    assert observer.latest["state"][0] == 0
    assert observer.rows_seen == 20
    assert (tmp_path / "truth.jsonl").read_bytes() == before


def test_observer_distinguishes_telemetry_trim_stale_estimate_and_hold(tmp_path, monkeypatch):
    (tmp_path / "manifest.json").write_text(json.dumps({"parameter_overlay": {"PTCH_TRIM_DEG": -3}}))
    append(tmp_path / "truth.jsonl", truth(1, hold=True))
    append(tmp_path / "telemetry.jsonl", {"fdm_time_s": 1, "message": {
        "mavpackettype": "ATTITUDE", "roll": .1, "pitch": .2, "yaw": .3}})
    append(tmp_path / "telemetry.jsonl", {"fdm_time_s": 1, "message": {
        "mavpackettype": "HEARTBEAT", "base_mode": 128, "custom_mode": 5}})
    observer = RunObserver(tmp_path)
    observer.poll()
    snapshot = observer.snapshot(active=True)
    assert snapshot["status"] == "BOOTSTRAP HOLD" and snapshot["armed"]
    assert snapshot["estimate"]["pitch_deg"] == pytest.approx(math.degrees(.2)-3)
    append(tmp_path / "truth.jsonl", truth(4))
    observer.poll()
    snapshot = observer.snapshot(active=True)
    assert snapshot["status"] == "FREE FLIGHT" and snapshot["estimate"] is None
    assert snapshot["stale_mavlink"]
    monkeypatch.setattr(sitl_session.time, "monotonic", lambda: observer.last_truth_wall+3)
    assert observer.snapshot(active=True)["stale_truth"]


def test_finished_report_and_cancelled_run_are_distinct(tmp_path):
    append(tmp_path / "truth.jsonl", truth())
    (tmp_path / "report.json").write_text('{"passed":true}')
    observer = RunObserver(tmp_path)
    observer.poll()
    assert observer.snapshot()["status"] == "PASSED"
    assert observer.snapshot(stopped=True)["status"] == "STOPPED"
    observer.report["passed"] = False
    assert observer.snapshot()["status"] == "FAILED"


def test_plot_body_pitch_uses_converted_estimate_and_actual_pwm():
    pytest.importorskip("plotly")
    estimate = body_attitude({"roll": .01, "pitch": .02, "yaw": 0}, -3)
    figures = sitl_figures({"history": [truth()], "estimates": [{"time_s": 1, **estimate}]})
    assert figures["Pitch"].data[1].y[0] == pytest.approx(math.degrees(.02)-3)
    assert figures["PWM"].data[1].y[0] == 1801


def test_manager_rejects_unknown_scenario_and_wrong_firmware(tmp_path):
    manager = SITLManager(binary=tmp_path / "dev", root=tmp_path / "runs")
    with pytest.raises(ValueError, match="supported"):
        manager.start("shell command")
    manager.binary.write_bytes(b"ArduPlane V4.8.0-dev (a4d4ec06)\x00")
    with pytest.raises(ValueError, match="4.7.1"):
        manager.start("fbwa")
    assert manager.process is None
    assert not manager.root.exists()


def test_manager_owns_process_group_and_cannot_start_a_second_run(tmp_path, monkeypatch):
    calls = []
    class Process:
        pid = 987654
        returncode = None
        def poll(self):
            return self.returncode
        def wait(self, timeout):
            self.returncode = 130
            return self.returncode
    process = Process()
    def popen(command, **kwargs):
        calls.append((command, kwargs))
        return process
    signals = []
    monkeypatch.setattr(sitl_session, "verify_firmware", lambda path: {})
    monkeypatch.setattr(sitl_session.subprocess, "Popen", popen)
    monkeypatch.setattr(sitl_session.os, "killpg", lambda pid, sig: signals.append((pid, sig)))
    manager = SITLManager(root=tmp_path / "runs")
    directory = manager.start("plus")
    try:
        with pytest.raises(RuntimeError, match="already active"):
            manager.start("fbwa")
        command, kwargs = calls[0]
        assert command[1:4] == ["-m", "spearhead.sitl.flight", "fbwa"]
        assert command[-4:] == ["--perturb-roll", ".5", "--perturb-pitch", ".5"]
        assert kwargs["start_new_session"] is True and "shell" not in kwargs
        assert directory.parent == manager.root
    finally:
        # A corrupted optional recording label must not strand firmware.
        (manager.root / (directory.name+".session.json")).write_text("{")
        manager.stop()
    assert signals[0] == (process.pid, sitl_session.signal.SIGINT)
    assert not manager.worker.is_alive()
    assert manager.snapshot()["status"] == "STOPPED"


def test_recording_load_is_read_only_and_cannot_escape_log_root(tmp_path):
    directory = tmp_path / "saved"
    directory.mkdir()
    append(directory / "truth.jsonl", truth())
    (directory / "manifest.json").write_text(json.dumps({
        "config": {"firmware": "Plane-4.7.1", "aerodynamics": "tail_update_2026_10_07"},
        "parameter_overlay": {"PTCH_TRIM_DEG": 0}}))
    (directory / "report.json").write_text('{"passed":true}')
    manager = SITLManager(root=tmp_path)
    before = (directory / "truth.jsonl").read_bytes()
    manager.load_completed("saved")
    assert manager.snapshot()["status"] == "PASSED"
    assert manager.process is None
    assert (directory / "truth.jsonl").read_bytes() == before
    with pytest.raises(ValueError, match="Unknown"):
        manager.load_completed("../saved")
    (tmp_path / "saved.session.json").write_text('{"cancelled":true}')
    manager.load_completed("saved")
    assert manager.snapshot()["status"] == "STOPPED"
