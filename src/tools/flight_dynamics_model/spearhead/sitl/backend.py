"""Deterministic lockstep session; repeated datagrams never advance physics."""

from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import re
import socket
import threading
import time

import numpy as np
import yaml

from . import DATASET, FIRMWARE_COMMIT, FIRMWARE_TAG
from .actuators import ActuatorMap, Channel
from .physics import Physics, full_trim
from .protocol import ServoFrame, encode_sensor
from ..aircraft import aircraft_summary, load_aircraft_config

CONFIG_ROOT = Path(__file__).resolve().parents[2] / "sitl_configs"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_firmware(path):
    identity = re.findall(rb"ArduPlane V[^\x00\r\n]{1,90}", Path(path).read_bytes())
    expected = b"ArduPlane V4.7.1 (dbe79216)"
    if expected not in identity or any(item not in (expected, b"ArduPlane V4.7.1") for item in identity):
        raise ValueError(f"SITL executable must identify as {expected.decode()}")
    return {"tag": FIRMWARE_TAG, "source_commit": FIRMWARE_COMMIT,
            "identity": expected.decode(), "binary_sha256": sha256(path)}


def load_setup(path=None):
    path = Path(path or CONFIG_ROOT / "fixed_wing_new.yaml").resolve()
    config = yaml.safe_load(path.read_text())
    required = {"firmware", "aircraft", "aerodynamics", "cg_body_xyz_m", "speed_m_s",
                "altitude_m", "max_step_s", "physics_rate_hz", "channels"}
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError("Invalid SITL configuration fields")
    if config["firmware"] != FIRMWARE_TAG or config["aerodynamics"] != DATASET:
        raise ValueError("This integration requires Plane-4.7.1 and the new AeroDB")
    cg = np.asarray(config["cg_body_xyz_m"], dtype=float)
    if cg.shape != (3,) or not np.isfinite(cg).all():
        raise ValueError("CG must contain three finite values")
    params = replace(load_aircraft_config(config["aircraft"], aerodynamics=DATASET), cg_body_xyz=cg)
    if (not params.use_aero_database or params.deflection_database_path is None
            or params.rho != 1.021 or not np.array_equal(params.aerodb_reference.moment_center_body_xyz, [0, 0, 0])):
        raise ValueError("Unexpected AeroDB, deflections, density, or source moment center")
    channels = dict(config["channels"])
    if set(channels) not in ({"aileron", "elevator", "rudder", "pusher"},
                             {"aileron", "elevator", "rudder", "pusher", "lift"}):
        raise ValueError("Configure all four conventional actuator mappings explicitly")
    mapping = ActuatorMap(**{k: tuple(Channel(**item) for item in v) if k == "lift" else Channel(**v)
                             for k, v in channels.items()})
    state, control = full_trim(params, config["speed_m_s"], config["altitude_m"])
    rate = config["physics_rate_hz"]
    if not isinstance(rate, int) or not 100 <= rate <= 1000:
        raise ValueError("Physics rate must be an integer in [100, 1000] Hz")
    physics = Physics(params, state, max_step=config["max_step_s"])
    physics.rate_hz = rate
    manifest = {"config": config, "config_sha256": sha256(path), "aircraft": aircraft_summary(params),
                "database_sha256": sha256(params.aero_database_path),
                "deflections_sha256": sha256(params.deflection_database_path),
                "initial_state": state.tolist(), "trim_control": asdict(control),
                "limitations": ["diagnostic CG; not a physical flight recommendation",
                                "synthetic servo travel, provisional rates/inertia/pusher",
                                "bootstrap hold is a protocol test, not controlled flight",
                                "AP atmosphere is standard; FDM rho is fixed at 1.021"]}
    source_root = Path(__file__).resolve().parents[1]
    manifest["fdm_source_sha256"] = {str(p.relative_to(source_root)): sha256(p)
                                     for p in sorted(source_root.rglob("*.py"))}
    return physics, mapping, manifest


class Session:
    def __init__(self, physics, mapping):
        self.physics, self.mapping = physics, mapping
        self.peer = self.last_frame = self.last_packet = self.last_reply = None
        self.accepted = self.duplicates = 0
        self.last_record = None
        self.release_requested = threading.Event()

    def release(self):
        """Explicit one-way transition; applied at the next physics boundary."""
        self.release_requested.set()

    def process(self, packet, peer):
        frame = ServoFrame.decode(packet)
        if self.peer is not None and peer != self.peer:
            raise ValueError("Unexpected UDP peer")
        if self.last_frame is not None:
            delta = (frame.count - self.last_frame.count) & 0xffffffff
            if delta == 0:
                if packet != self.last_packet:
                    raise ValueError("Duplicate counter contains different controls")
                self.duplicates += 1
                return self.last_reply
            if delta != 1:
                # Lockstep cannot skip a frame; restart requires a fresh session.
                raise ValueError("Out-of-order, missing frame, or firmware restart; start a new run")
        control = self.mapping.control(frame.pwm)
        if self.release_requested.is_set():
            self.physics.hold = False
        # JSON 4.7.1 advertises its adaptive rate (which tends toward the
        # controller rate). Sensor timestamps own the external physics clock.
        sensor = self.physics.advance(1 / getattr(self.physics, "rate_hz", 200), control)
        reply = encode_sensor(sensor)
        self.peer, self.last_frame = peer, frame
        self.last_packet, self.last_reply = packet, reply
        self.accepted += 1
        speed = np.linalg.norm(self.physics.state[3:6])
        alpha = float(np.rad2deg(np.arctan2(self.physics.state[5], self.physics.state[3])))
        beta = float(np.rad2deg(np.arcsin(np.clip(self.physics.state[4]/max(speed, 1e-12), -1, 1))))
        self.last_record = {"frame": frame.count, "advertised_rate_hz": frame.rate,
                            "physics_rate_hz": getattr(self.physics, "rate_hz", 200), "time_s": self.physics.time,
                            "hold": self.physics.hold, "pwm": frame.pwm, "control": asdict(control),
                            "state": self.physics.state.tolist(), "sensor": sensor,
                            "control_coverage_exceeded": not (-5 <= alpha <= 0 and abs(beta) < 1e-6),
                            "alpha_deg": alpha, "beta_deg": beta}
        return reply


def serve(sock, session, output, stop, *, timeout_s=10., duration_s=None, release_file=None):
    """Return on duration/stop; fail explicitly on invalid traffic or silence."""
    sock.settimeout(.2)
    last_rx = time.monotonic()
    with Path(output).open("x") as log:
        while not stop.is_set():
            try:
                packet, peer = sock.recvfrom(256)
            except socket.timeout:
                if time.monotonic() - last_rx > timeout_s:
                    raise TimeoutError("No valid SITL control packet within timeout")
                continue
            before = session.accepted
            if release_file is not None and Path(release_file).is_file():
                if Path(release_file).read_text().strip() != "release":
                    raise ValueError("Release file must contain the explicit word 'release'")
                session.release()
            reply = session.process(packet, peer)
            sock.sendto(reply, peer)
            last_rx = time.monotonic()
            if session.accepted != before:
                log.write(json.dumps(session.last_record, allow_nan=False) + "\n")
            if duration_s is not None and session.physics.time >= duration_s:
                return
