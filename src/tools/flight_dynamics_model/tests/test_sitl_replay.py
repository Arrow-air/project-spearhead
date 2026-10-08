"""Recorded PWM independently reproduces physics and detects damaged histories."""

import json

import pytest

from spearhead.sitl.backend import CONFIG_ROOT, Session, load_setup
from spearhead.sitl.replay import replay
from test_sitl import packet


@pytest.fixture
def recorded(tmp_path):
    physics, mapping, manifest = load_setup(CONFIG_ROOT / "fixed_wing_control.yaml")
    physics.hold = False
    session = Session(physics, mapping)
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    rows = []
    for count in range(20):
        session.process(packet(count), ("127.0.0.1", 1234))
        rows.append(session.last_record)
    (tmp_path / "truth.jsonl").write_text("".join(json.dumps(r)+"\n" for r in rows))
    return tmp_path, rows


def test_pwm_replay_matches_original_physics_without_firmware(recorded):
    path, _ = recorded
    result = replay(path)
    assert result["passed"] and result["physics_source_hashes_match"]
    assert result["replayed_frames"] == 19
    assert result["max_position_error_m"] < 1e-10
    assert result["max_attitude_error_deg"] < 1e-10


def test_replay_detects_modified_pwm(recorded):
    path, rows = recorded
    for r in rows[1:]:
        r["pwm"] = list(r["pwm"])
        r["pwm"][1] = 2000
    (path / "truth.jsonl").write_text("".join(json.dumps(r)+"\n" for r in rows))
    assert not replay(path)["passed"]


def test_replay_rejects_missing_frame_and_changed_source_data(recorded):
    path, rows = recorded
    (path / "truth.jsonl").write_text("".join(json.dumps(r)+"\n" for r in rows if r["frame"] != 2))
    with pytest.raises(ValueError, match="contiguous"):
        replay(path)
    manifest = json.loads((path / "manifest.json").read_text())
    manifest["database_sha256"] = "changed"
    (path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="source data mismatch"):
        replay(path)
