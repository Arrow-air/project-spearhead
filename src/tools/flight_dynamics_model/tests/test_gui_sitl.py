"""SITL display isolation, incremental logs and owned process lifecycle."""

import json
import math
from pathlib import Path

import pytest
import numpy as np

from spearhead.gui import sitl_session
from spearhead.gui.sitl_page import (sitl_figures, control_coverage_text, transition_result_text,
                                     arming_badge, startup_activity, flight_story, SimulationPace, transition_metrics)
from spearhead.gui.sitl_session import JsonlTail, RunObserver, SITLManager, body_attitude
from spearhead.gui.external_view import camera_pose, scene_pose, trail_points


def truth(t=1, hold=False):
    return {"time_s": t, "state": [0, 0, -100, 28, 0, 0, .01, -.02, 0, 0, 0, 0],
            "control": {"de": math.radians(-3), "da": 0., "dr": 0., "throttle": .3},
            "hold": hold, "pwm": [1500, 1801, 1254, 1495], "physics_rate_hz": 200,
            "alpha_deg": -.1, "beta_deg": 0, "control_coverage_exceeded": False}


@pytest.mark.parametrize("angles, body_axis, expected", [
    ([0, 0, 0], [1, 0, 0], [1, 0, 0]),
    ([0, 30, 0], [1, 0, 0], [np.sqrt(3)/2, 0, .5]),  # nose up
    ([30, 0, 0], [0, 1, 0], [0, -np.sqrt(3)/2, -.5]),  # right wing down
    ([0, 0, 90], [1, 0, 0], [0, -1, 0]),  # east is scene -Y
])
def test_external_view_preserves_physical_attitude_signs(angles, body_axis, expected):
    state = truth()["state"]
    state[6:9] = np.deg2rad(angles).tolist()
    matrix = np.array(scene_pose(state)["rotation"])
    np.testing.assert_allclose(matrix @ body_axis, expected, atol=1e-12)
    np.testing.assert_allclose(matrix.T @ matrix, np.eye(3), atol=1e-12)
    assert np.linalg.det(matrix) == pytest.approx(1)


def test_external_view_rebases_large_positions_without_changing_altitude_or_path():
    rows = []
    for i in range(1000):
        row = truth(i*.05, hold=i < 5)
        row["state"][:3] = [1e7+i, 2e7+2*i, -100-i*.01]
        rows.append(row)
    current = rows[-1]
    pose = scene_pose(current["state"])
    points = trail_points(rows, current)
    assert pose["position"] == pytest.approx([0, 0, 109.99])
    assert len(points) <= 121
    assert points[-1] == pytest.approx(pose["position"])
    assert all(p[0] <= 0 and p[1] >= 0 for p in points)
    assert np.linalg.norm(np.array(points[-1])-points[0]) < 1500
    assert all(0 <= v < 25 for v in pose["ground_offset"][:2])
    assert trail_points(rows[:5], rows[4]) == []
    # Neither coordinate conversion nor display history changes FDM samples.
    assert current["state"][:3] == [10000999., 20001998., -109.99]


def test_external_camera_follows_heading_with_world_up_and_rejects_bad_state():
    state = truth()["state"]
    state[8] = 0
    north = camera_pose(state, "chase")
    state[8] = math.pi/2
    east = camera_pose(state, "chase")
    assert east == pytest.approx([north[1], -north[0], north[2]])
    assert east[2] > -state[2]
    state[0] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        scene_pose(state)


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


def test_wall_time_is_live_only_and_independent_of_simulation_clock(tmp_path, monkeypatch):
    monkeypatch.setattr(sitl_session.time, "monotonic", lambda: 100)
    observer = RunObserver(tmp_path)
    monkeypatch.setattr(sitl_session.time, "monotonic", lambda: 120)
    append(tmp_path / "truth.jsonl", truth(1, hold=True))
    observer.poll()
    live = observer.snapshot(active=True)
    assert live["wall_elapsed_s"] == 20 and live["time_s"] == 1
    assert live["truth_age_s"] == 0
    assert observer.snapshot()["wall_elapsed_s"] is None


@pytest.mark.parametrize("changes,title", [
    ({"truth": None}, "Preparing aircraft model"),
    ({"truth": None, "manifest": {"config": {}}}, "Starting flight engine"),
    ({"health": {"gps_fix": 1}}, "Acquiring GPS"),
    ({"health": {"gps_fix": 3, "ekf_ok": False}}, "Aligning navigation"),
    ({"health": {"gps_fix": 3, "ekf_ok": True}}, "Pre-arm checks"),
    ({"armed": True, "manifest": {"config": {"flight_regime": "hover"}}}, "Spooling lift motors"),
    ({"stale_truth": True}, "Waiting for simulator updates"),
])
def test_startup_feedback_follows_observed_state(changes, title):
    view = {"active": True, "truth": truth(20, hold=True), "wall_elapsed_s": 50,
            "free_time_s": None, "armed": False, **changes}
    activity = startup_activity(view)
    assert activity["title"] == title
    assert activity["elapsed_s"] == 50 and "estimate" in activity


@pytest.mark.parametrize("changes", [{"active": False}, {"free_time_s": 0}, {"error": "Failed"}])
def test_startup_feedback_stops_at_release_completion_or_failure(changes):
    assert startup_activity({"active": True, "free_time_s": None, **changes}) is None


def test_overdue_startup_does_not_show_a_zero_countdown():
    activity = startup_activity({"active": True, "truth": truth(60, hold=True), "wall_elapsed_s": 140})
    assert "Taking longer" in activity["estimate"]


def test_arming_badge_preserves_unknown_stale_and_recorded_states():
    assert arming_badge({}) == ("ARMING UNKNOWN", "unknown")
    assert arming_badge({"active": True, "armed": False}) == ("DISARMED", "disarmed")
    assert arming_badge({"active": True, "armed": True}) == ("ARMED", "armed")
    assert arming_badge({"active": True, "armed": True, "stale_mavlink": True}) == ("ARMED · STALE", "unknown")
    assert arming_badge({"armed": True, "run_id": "saved"}) == ("RECORDED · ARMED", "armed")


def test_flight_story_never_treats_handover_or_stopped_run_as_a_pass():
    view = {"active": True, "free_time_s": 20, "mode": 5, "truth": {"pwm": [1500]*8},
            "manifest": {"test_kind": "quad_transition", "parameter_overlay": {
                "Q_M_PWM_MIN": 1100, "Q_M_PWM_MAX": 1940, "Q_M_SPIN_ARM": .1}}}
    assert flight_story(view)[0] == "Transitioning to forward flight"
    view["truth"]["pwm"][4:] = [1184]*4
    title, detail, tone, step = flight_story(view)
    assert "lift motors at idle" in title
    assert "does not mean" in detail and tone == "neutral" and step == 1
    view.update(active=False, status="STOPPED", report={"passed": True})
    assert flight_story(view)[0] == "Flight stopped"


def test_flight_story_needs_measured_idle_parameters_and_fresh_telemetry():
    view = {"active": True, "free_time_s": 20, "mode": 5, "truth": {"pwm": [1100]*8},
            "manifest": {"test_kind": "quad_transition"}}
    assert flight_story(view)[0] == "Transitioning to forward flight"
    view["stale_mavlink"] = True
    assert flight_story(view)[0] == "Waiting for fresh flight data"
    view["error"] = "Firmware exited"
    assert flight_story(view)[0] == "Flight needs attention"


def test_recorded_interactive_completion_is_distinct_from_test_acceptance():
    view = {"active": False, "report": {"test_kind": "interactive_session", "passed": True}}
    assert flight_story(view)[0] == "Flight complete"
    view["report"] = {"test_kind": "quad_transition", "passed": False, "firmware_transition_done": True}
    assert flight_story(view)[0] == "Simulation checks failed"
    view["report"].update(passed=True, free_duration_s=60, longest_forward_without_lift_s=33.5)
    title, detail, tone, step = flight_story(view)
    assert title == "Simulation checks passed" and "33.5 s" in detail
    assert tone == "success" and step == 2


def test_simulation_pace_measures_wall_time_and_resets_between_runs():
    pace = SimulationPace()
    view = {"active": True, "run_id": "a", "time_s": 10}
    assert pace.update(view, 100) is None
    assert pace.update({**view, "time_s": 10.5}, 102) == .25
    assert pace.update({**view, "run_id": "b", "time_s": 0}, 103) is None
    assert pace.update({**view, "active": False}, 104) is None


def test_simulation_pace_reports_a_stalled_clock_without_guessing_progress():
    pace = SimulationPace()
    view = {"active": True, "run_id": "a", "time_s": 10}
    pace.update(view, 100)
    assert pace.update(view, 102) == 0
    assert pace.update({**view, "time_s": 12}, 103) == pytest.approx(2/3)


def test_plot_shows_only_fdm_attitude_and_actual_pwm():
    pytest.importorskip("plotly")
    estimate = body_attitude({"roll": .01, "pitch": .02, "yaw": 0}, -3)
    figures = sitl_figures({"history": [truth()], "estimates": [{"time_s": 1, **estimate}]})
    assert len(figures["Pitch"].data) == len(figures["Roll"].data) == 1
    assert figures["Pitch"].data[0].y[0] == pytest.approx(math.degrees(-.02))
    assert figures["PWM"].data[1].y[0] == 1801


def test_transition_plot_uses_forward_trim_instead_of_hover_neutral():
    pytest.importorskip("plotly")
    figures = sitl_figures({"history": [truth()], "manifest": {
        "trim_control": {"de": 0}, "forward_trim": {"control": {"de": math.radians(-8.2)}}}})
    assert figures["Elevator"].layout.shapes[0].y0 == pytest.approx(-8.2)


def test_transition_summary_does_not_equate_firmware_handover_with_pass():
    text = transition_result_text({"passed": False, "firmware_transition_done": True,
                                  "longest_forward_without_lift_s": .135, "free_duration_s": 29.045,
                                  "final_true_airspeed_m_s": 19})
    assert "FAILED" in text and "handover: yes" in text and "requires 10 s" in text
    assert "29.0 s" in text


def test_coverage_readout_uses_recorded_surface_bounds_and_velocities():
    text = control_coverage_text({"control_data": {"beta_deg": 0, "coverage": {
        "elevator": {"deflection_deg": [-15, 15], "measured_velocities_m_s": [10, 15, 25]},
        "aileron": {"deflection_deg": [-30, 30], "measured_velocities_m_s": [25]}}}})
    assert "elevator -15…15° at 10/15/25 m/s" in text
    assert "aileron -30…30° at 25 m/s" in text and "depends on α and speed" in text


def test_recording_symlink_cannot_escape_root(tmp_path):
    root = tmp_path / "runs"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "report.json").write_text('{"passed":true}')
    (root / "linked").symlink_to(outside, target_is_directory=True)
    manager = SITLManager(root=root)
    assert manager.completed_runs() == {}
    with pytest.raises(ValueError, match="Unknown"):
        manager.load_completed("linked")


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


def test_quad_interactive_preset_uses_owned_quad_runner(tmp_path, monkeypatch):
    from types import SimpleNamespace
    calls = []
    process = SimpleNamespace(pid=987654, returncode=None, poll=lambda: None)
    monkeypatch.setattr(sitl_session, "verify_firmware", lambda path: {})
    monkeypatch.setattr(sitl_session.subprocess, "Popen", lambda command, **kwargs: calls.append((command, kwargs)) or process)
    monkeypatch.setattr(sitl_session.threading.Thread, "start", lambda self: None)
    manager = SITLManager(root=tmp_path / "runs")
    manager.start("quad_interactive")
    command, kwargs = calls[0]
    assert command[1:4] == ["-m", "spearhead.sitl.quad", "--binary"]
    assert command[-3:] == ["--interactive", "--duration", "600"]
    assert kwargs["start_new_session"] is True


def test_quad_controls_ready_and_interactive_completion_are_distinct_from_validation(tmp_path):
    append(tmp_path / "truth.jsonl", truth())
    (tmp_path / "manifest.json").write_text(json.dumps({"interactive": True, "requires_controls_ready": True}))
    observer = RunObserver(tmp_path); observer.poll()
    assert not observer.snapshot(active=True)["mode_controls_ready"]
    (tmp_path / "controls-ready.json").write_text('{"time_s": 1}')
    assert observer.snapshot(active=True)["mode_controls_ready"]
    (tmp_path / "report.json").write_text('{"passed":true,"test_kind":"interactive_session","validation_run":false}')
    observer.poll()
    assert observer.snapshot()["status"] == "COMPLETED"


def test_transition_extremes_survive_bounded_display_history(tmp_path):
    manifest = {'test_kind':'quad_transition','initial_state':[0,0,-80], 'parameter_overlay':{}}
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    rows = []
    for t, pitch, altitude in [(1,18,80),(2,5,79),(12,-24,65),(30,1,80)]:
        row = truth(t)
        row['state'][7] = math.radians(pitch)
        row['state'][2] = -altitude
        rows.append(row)
    (tmp_path/'truth.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    observer = RunObserver(tmp_path,history_limit=1)
    observer.poll()
    view = observer.snapshot()
    assert len(view['history']) == 1
    metrics = transition_metrics(view)
    assert metrics['peak_abs_pitch_deg'] == pytest.approx(24)
    assert metrics['altitude_loss_m'] == 15
    assert metrics['hover_peak_abs_pitch_deg'] == pytest.approx(18)
    assert metrics['longest_forward_without_lift_s'] is None
    view['transient_metrics']['altitude_loss_m'] = 999
    assert observer.snapshot()['transient_metrics']['altitude_loss_m'] == 15


def test_diagnostic_plots_use_physical_units_and_hide_low_speed_alpha():
    row = truth()
    row['state'][3:6] = [20,0,2]
    row['state'][10] = math.radians(-12)
    low = truth(2)
    low['state'][3:6] = [0,0,0]
    figures = sitl_figures({'history':[row,low], 'manifest':{'trim_control':{'de':0},
        'control_model':{'assumed_travel_deg':{'elevator':20,'aileron':20,'rudder':20}}}})
    assert figures['Alpha'].data[0].y == pytest.approx([math.degrees(math.atan2(2,20)),None])
    assert figures['Pitch rate'].data[0].y[0] == pytest.approx(-12)
    assert figures['Pusher'].data[0].y[0] == 30
    assert sorted(s.y0 for s in figures['Surface travel'].layout.shapes) == [-20,20]
