import json
from collections import deque
from types import SimpleNamespace

import pytest

from spearhead.sitl import ground_control
from spearhead.sitl.ground_control import ModeCommands
from spearhead.gui.ground_status import StatusTexts, firmware_health
from spearhead.gui.sitl_session import RunObserver, SITLManager


def test_mode_command_uses_real_mavlink_ack_then_heartbeat_confirmation(tmp_path, monkeypatch):
    monkeypatch.setattr(ground_control.time, "monotonic", lambda: 100)
    calls, events = [], []
    link = SimpleNamespace(target_system=1, target_component=1,
                           mav=SimpleNamespace(command_long_send=lambda *args: calls.append(args)))
    (tmp_path / "one.json").write_text(json.dumps({"id": "one", "kind": "set_mode", "mode": 0, "created_at": 99}))
    commands = ModeCommands(tmp_path, lambda **event: events.append(event))
    commands.update(link, ready=True)
    assert calls == [(1, 1, 176, 0, 1, 0, 0, 0, 0, 0, 0)]
    commands.update(link, {"mavpackettype": "COMMAND_ACK", "command": 176, "result": 0}, ready=True)
    assert events[-1]["status"] == "accepted" and commands.pending is not None
    commands.update(link, {"mavpackettype": "HEARTBEAT", "custom_mode": 5}, ready=True)
    assert commands.pending is not None
    commands.update(link, {"mavpackettype": "HEARTBEAT", "custom_mode": 0}, ready=True)
    assert events[-1]["status"] == "confirmed" and commands.pending is None
    assert len(calls) == 1


@pytest.mark.parametrize("result,label", [(1, "Temporarily rejected"), (2, "Denied"), (3, "Unsupported"), (4, "Failed")])
def test_firmware_rejections_are_visible(tmp_path, monkeypatch, result, label):
    monkeypatch.setattr(ground_control.time, "monotonic", lambda: 100)
    events = []
    commands = ModeCommands(tmp_path, lambda **event: events.append(event))
    commands.pending = {"id": "one", "mode": 5, "sent_at": 99}
    commands.update(None, {"mavpackettype": "COMMAND_ACK", "command": 176, "result": result})
    assert events[-1]["status"] == "rejected" and events[-1]["text"] == label
    assert commands.pending is None


def test_timeout_and_expired_request_never_report_a_success(tmp_path, monkeypatch):
    monkeypatch.setattr(ground_control.time, "monotonic", lambda: 100)
    events = []
    commands = ModeCommands(tmp_path, lambda **event: events.append(event))
    commands.pending = {"id": "old", "mode": 5, "sent_at": 90}
    (tmp_path / "expired.json").write_text(json.dumps({"id": "expired", "kind": "set_mode", "mode": 5, "created_at": 1}))
    commands.update(None, ready=True)
    assert [e["status"] for e in events] == ["timeout", "rejected"]
    assert commands.pending is None


def test_firmware_text_chunks_keep_severity_timestamp_and_literal_text():
    feed = deque(maxlen=2)
    texts = StatusTexts(feed)
    texts.add({"text": "A"*50, "severity": 3, "id": 12, "chunk_seq": 0}, 2)
    assert not feed
    texts.add({"text": "<tail>", "severity": 3, "id": 12, "chunk_seq": 1}, 3)
    assert feed[-1] == {"time_s": 2, "severity": 3, "level": "Error", "source": "ArduPilot", "text": "A"*50+"<tail>"}
    texts.add({"text": "PreArm: waiting", "severity": 4}, 4)
    texts.add({"text": "Ready", "severity": 6}, 5)
    assert len(feed) == 2 and feed[0]["level"] == "Warning"


def test_unknown_stale_and_unhealthy_telemetry_are_not_green():
    assert firmware_health({}, 1)["ekf_ok"] is None
    messages = {"SYS_STATUS": {"fdm_time_s": 1, "message": {
        "onboard_control_sensors_present": 33, "onboard_control_sensors_enabled": 33,
        "onboard_control_sensors_health": 1, "voltage_battery": 65535, "current_battery": -1, "battery_remaining": -1}},
        "EKF_STATUS_REPORT": {"fdm_time_s": 1, "message": {"flags": 63}},
        "GPS_RAW_INT": {"fdm_time_s": 1, "message": {"fix_type": 3, "satellites_visible": 255}}}
    health = firmware_health(messages, 1, live=True)
    assert health["sensor_faults"] == ["GPS"] and health["ekf_ok"]
    assert health["battery_v"] is None and health["satellites"] is None
    stale = firmware_health(messages, 5, live=True)
    assert stale["ekf_ok"] is None and not stale["system_known"]
    messages["SYS_STATUS"]["message"].update(onboard_control_sensors_present=4194304,
        onboard_control_sensors_enabled=4194304, onboard_control_sensors_health=0, voltage_battery=0)
    unavailable = firmware_health(messages, 1)
    assert unavailable["sensor_faults"] == ["Terrain data"] and unavailable["battery_v"] is None


def test_manager_routes_commands_only_to_live_interactive_run(tmp_path):
    manager = SITLManager(root=tmp_path)
    with pytest.raises(RuntimeError, match="live interactive"):
        manager.set_mode(5)
    manager.observer = RunObserver(tmp_path)
    manager.process = SimpleNamespace(poll=lambda: None)
    manager.cached.update(active=True, interactive=True, status="FREE FLIGHT", stale_truth=False,
                          stale_mavlink=False, error=None)
    ident = manager.set_mode(0)
    request = json.loads((tmp_path / "commands" / (ident+".json")).read_text())
    assert request["mode"] == 0 and request["kind"] == "set_mode"
    assert not list((tmp_path / "commands").glob("*.tmp"))
    with pytest.raises(RuntimeError, match="current mode change"):
        manager.set_mode(5)
    with pytest.raises(ValueError, match="supported"):
        manager.set_mode(19)  # No QuadPlane setup in this reference simulation.
    manager.pending_command = None
    manager.cached["stale_mavlink"] = True
    with pytest.raises(RuntimeError, match="live interactive"):
        manager.set_mode(5)


def test_observer_preserves_multiple_warnings_and_command_outcomes(tmp_path):
    rows = [{"fdm_time_s": n, "message": {"mavpackettype": "STATUSTEXT", "severity": 4, "text": text}}
            for n, text in enumerate(["PreArm: GPS", "EKF waiting"])]
    (tmp_path / "telemetry.jsonl").write_text("".join(json.dumps(row)+"\n" for row in rows))
    (tmp_path / "control-events.jsonl").write_text(json.dumps({"command_id": "one", "mode": 5,
        "time_s": 2, "status": "rejected", "text": "Denied"})+"\n")
    observer = RunObserver(tmp_path)
    observer.poll()
    view = observer.snapshot()
    assert [m["text"] for m in view["notifications"]] == ["PreArm: GPS", "EKF waiting", "Denied"]
    assert view["command"]["status"] == "rejected"
    observer.poll()
    assert len(observer.snapshot()["notifications"]) == 3


def test_fast_display_snapshot_is_bounded_and_does_not_mutate_cached_history():
    manager = SITLManager()
    manager.cached["history"] = [{"time_s": i, "state": [i]} for i in range(1000)]
    manager.cached["estimates"] = [{"time_s": 1}]
    display = manager.snapshot(include_history=False)
    assert len(display["history"]) == 600 and not display["estimates"]
    display["history"][-1]["state"][0] = -1
    full = manager.snapshot()
    assert len(full["history"]) == 1000 and full["history"][-1]["state"] == [999]


def test_pinned_471_catalog_includes_every_user_mode_and_internal_status():
    # Published enum in pinned mode.h; 9 is unused and 16 is internal startup.
    assert set(ground_control.MODE_NAMES) == set(range(27)) - {9}
    assert set(ground_control.QUADPLANE_MODES) == set(range(26)) - {9, 16}
    assert ground_control.MODE_NAMES[16] == "INITIALISING"
    assert 26 in ground_control.PLANE_MODES and 26 not in ground_control.QUADPLANE_MODES
    assert set(ground_control.MODE_HELP) == set(ground_control.MODE_NAMES) - {16}
    assert ground_control.session_modes({"parameter_overlay": {"Q_ENABLE": 1}}) == ground_control.QUADPLANE_MODES
    assert 19 not in ground_control.session_modes({"parameter_overlay": {"Q_ENABLE": 0}})


@pytest.mark.parametrize("mode", list(ground_control.QUADPLANE_MODES))
def test_all_quadplane_modes_route_to_firmware_and_require_matching_heartbeat(tmp_path, monkeypatch, mode):
    monkeypatch.setattr(ground_control.time, "monotonic", lambda: 100)
    calls, events = [], []
    link = SimpleNamespace(target_system=1, target_component=1,
                           mav=SimpleNamespace(command_long_send=lambda *args: calls.append(args)))
    (tmp_path / "mode.json").write_text(json.dumps({"id": "test", "kind": "set_mode", "mode": mode, "created_at": 99}))
    commands = ModeCommands(tmp_path, lambda **event: events.append(event), ground_control.QUADPLANE_MODES)
    commands.update(link, ready=True)
    assert calls == [(1, 1, 176, 0, 1, mode, 0, 0, 0, 0, 0)]
    commands.update(link, {"mavpackettype": "COMMAND_ACK", "command": 176, "result": 0}, ready=True)
    assert commands.pending is not None  # ACK alone never confirms the selection.
    commands.update(link, {"mavpackettype": "HEARTBEAT", "custom_mode": mode}, ready=True)
    assert events[-1]["status"] == "confirmed" and ground_control.MODE_NAMES[mode] in events[-1]["text"]
    assert commands.pending is None


@pytest.mark.parametrize("mode", [9, 16, 26, 27])
def test_quadplane_rejects_non_quadplane_and_internal_modes(tmp_path, mode):
    events = []
    (tmp_path / "mode.json").write_text(json.dumps({"id": "test", "kind": "set_mode", "mode": mode, "created_at": ground_control.time.monotonic()}))
    commands = ModeCommands(tmp_path, lambda **event: events.append(event), ground_control.QUADPLANE_MODES)
    commands.update(None, ready=True)
    assert events[-1]["status"] == "rejected" and commands.pending is None


def test_quad_manager_waits_for_bootstrap_and_queues_qloiter_atomically(tmp_path):
    manager = SITLManager(root=tmp_path)
    manager.observer = RunObserver(tmp_path)
    manager.process = SimpleNamespace(poll=lambda: None)
    manager.cached.update(active=True, interactive=True, status="FREE FLIGHT", stale_truth=False,
                          stale_mavlink=False, error=None, mode_controls_ready=False,
                          manifest={"parameter_overlay": {"Q_ENABLE": 1}})
    with pytest.raises(RuntimeError, match="live interactive"):
        manager.set_mode(19)
    manager.cached["mode_controls_ready"] = True
    ident = manager.set_mode(19)
    assert json.loads((tmp_path / "commands" / (ident+".json")).read_text())["mode"] == 19
    assert manager.pending_command["text"] == "Queued QLOITER"
    with pytest.raises(ValueError, match="supported"):
        manager.set_mode(26)
