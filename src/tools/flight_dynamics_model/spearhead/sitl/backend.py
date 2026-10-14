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
    if not isinstance(config, dict) or not required <= set(config) or set(config)-required-{"control_deflection_policy", "flight_regime", "ground_contact", "hover_trim", "wind", "failure"}:
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
    policy = config.get("control_deflection_policy", "hold")
    if policy not in ("hold", "linear"):
        raise ValueError("Control deflection policy must be hold or linear")
    params = replace(params, control_deflection_policy=policy,
                     control_surface_limits_deg=tuple(getattr(mapping, name).travel_deg
                                                      for name in ("elevator", "aileron", "rudder")))
    regime = config.get("flight_regime", "fixed_wing")
    hover_trim = config.get("hover_trim", "collective")
    if hover_trim not in ("collective", "independent_rotors") or ("hover_trim" in config and regime not in ("hover", "ground")):
        raise ValueError("Hover trim must select collective or independent_rotors for a hover or ground start")
    ground = None
    if regime in ("hover", "ground"):
        from ..propulsion import hover_collective, hover_trim_control
        from ..controls import Control
        from ..dynamics import dynamics
        if len(mapping.lift) != 4 or config["speed_m_s"] != 0 or not np.isfinite(config["altitude_m"]) or config["altitude_m"] <= 0:
            raise ValueError("Hover requires four lift channels, zero initial speed and positive altitude")
        state = np.zeros(12)
        state[2] = -config["altitude_m"]
        control = (hover_trim_control(params) if hover_trim == "independent_rotors" else
                   Control(0, 0, 0, 0, lift_throttle=hover_collective(params)))
        if np.linalg.norm(dynamics(0, state, params, control)[3:]) > 1e-6:
            raise ValueError("Symmetric collective does not balance the configured hover CG")
        if regime == "ground":
            from .ground import GroundContact
            if config["aircraft"] != "pt15" or not 3 <= config["altitude_m"] <= 30:
                raise ValueError("Ground start supports PT1 takeoff targets of 3–30 m")
            ground = GroundContact(**config.get("ground_contact", {}))
            state = ground.rest_state(params)
    elif regime == "fixed_wing":
        state, control = full_trim(params, config["speed_m_s"], config["altitude_m"])
    else:
        raise ValueError("Unknown flight regime")
    rate = config["physics_rate_hz"]
    if not isinstance(rate, int) or not 100 <= rate <= 1000:
        raise ValueError("Physics rate must be an integer in [100, 1000] Hz")
    if regime != "ground" and "ground_contact" in config:
        raise ValueError("Ground-contact settings require a ground initial point")
    from .battery import BatteryModel
    from ..wind import WindField
    wind = WindField.from_config(config.get("wind"))
    physics = Physics(params, state, hold=regime != "ground", max_step=config["max_step_s"], ground=ground,
                      battery=BatteryModel.from_params(params) if params.lift_system is not None else None,
                      wind=wind)
    from .failures import Failure
    physics.failure = Failure.from_config(config.get("failure"))
    if physics.failure is not None and physics.failure.kind.startswith("lift_") and not mapping.lift:
        raise ValueError("Lift-motor failure requires mapped lift motors")
    physics.rate_hz = rate
    manifest = {"config": config, "config_sha256": sha256(path), "aircraft": aircraft_summary(params),
                "database_sha256": sha256(params.aero_database_path),
                "deflections_sha256": sha256(params.deflection_database_path),
                "initial_state": state.tolist(), "trim_control": asdict(control),
                "control_model": {"deflection_policy": policy,
                                  "assumed_travel_deg": {name: getattr(mapping, name).travel_deg
                                                         for name in ("elevator", "aileron", "rudder")},
                                  "travel_calibrated": False,
                                  "extrapolation": "continue measured endpoint slopes; alpha/beta held"},
                "failure": None if physics.failure is None else {**physics.failure.as_config(), "time_reference": "seconds after release (ground: after lift spin-up)"},
                "wind": None if wind is None or wind.calm else {**wind.as_config(), "time_reference": "seconds after release (ground: after lift spin-up)",
                                                                 "model": "uniform steady wind plus 1-cosine gusts; no turbulence or shear"},
                "limitations": ["simulation CG (target, not measured); not a physical flight recommendation",
                                "assumed unmeasured servo travel; provisional rates/inertia/pusher",
                                "bootstrap hold is a protocol test, not controlled flight",
                                "AP atmosphere is standard; FDM rho is fixed at 1.021"]}
    source_root = Path(__file__).resolve().parents[1]
    if ground is not None:
        manifest["ground_contact"] = {"provisional": True, "datum_ned_z_m": 0, **asdict(ground)}
        manifest["limitations"] = [item for item in manifest["limitations"] if "bootstrap" not in item and "simulation CG" not in item]
        manifest["limitations"].append("PT1 CG, inertia and motor coordinates are provisional; simulation does not validate hardware")
        manifest["limitations"].append("Assumed gear footprint/compliance/friction; flat ground only; no rotor ground effect")
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
        self.pending_air_start = None

    def air_start(self, state):
        """Queue a logged bootstrap teleport, permitted only before release."""
        state = np.asarray(state, dtype=float)
        if state.shape != (12,) or not np.isfinite(state).all() or state[2] >= 0:
            raise ValueError("Air-start needs a finite airborne 12-state vector")
        if not self.physics.hold or self.release_requested.is_set():
            raise ValueError("Cannot teleport after flight release")
        self.pending_air_start = state.copy()

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
        air_start_applied = self.pending_air_start is not None
        if air_start_applied:
            self.physics.state = self.pending_air_start
            self.pending_air_start = None
        released = self.release_requested.is_set() and self.physics.hold
        if released:
            self.physics.hold = False
        # JSON 4.7.1 advertises its adaptive rate (which tends toward the
        # controller rate). Sensor timestamps own the external physics clock.
        failure = getattr(self.physics, "failure", None)
        commanded = control
        if failure is not None:
            control = failure.apply(control, self.physics)
        sensor = self.physics.advance(1 / getattr(self.physics, "rate_hz", 200), control)
        reply = encode_sensor(sensor)
        self.peer, self.last_frame = peer, frame
        self.last_packet, self.last_reply = packet, reply
        self.accepted += 1
        wind_ned = self.physics.wind_ned(self.physics.time)
        v_air = self.physics.state[3:6].copy()
        if wind_ned is not None:
            from ..rotations import euler_body_to_ned
            v_air -= euler_body_to_ned(*self.physics.state[6:9]).T @ wind_ned
        speed = np.linalg.norm(v_air)
        alpha = float(np.rad2deg(np.arctan2(v_air[2], v_air[0])))
        beta = float(np.rad2deg(np.arcsin(np.clip(v_air[1]/max(speed, 1e-12), -1, 1))))
        from ..deflections import load_deflection_database_cached
        coverage = None
        params = self.physics.params
        if params.deflection_database_path is not None and speed >= 1e-6:
            _, coverage = load_deflection_database_cached(params.deflection_database_path, params.aerodb_reference).increments_at(
                alpha, beta, control, deflection_policy=params.control_deflection_policy,
                limits_deg=params.control_surface_limits_deg, velocity_m_s=float(speed))
        self.last_record = {"frame": frame.count, "advertised_rate_hz": frame.rate,
                            "physics_rate_hz": getattr(self.physics, "rate_hz", 200), "time_s": self.physics.time,
                            "hold": self.physics.hold, "pwm": frame.pwm, "control": asdict(control),
                            "bootstrap_air_start": air_start_applied,
                            "release_applied": released,
                            "state": self.physics.state.tolist(), "sensor": sensor,
                            "control_coverage": coverage,
                            "control_coverage_exceeded": coverage is not None and not coverage["within_measured_coverage"],
                            "alpha_deg": alpha, "beta_deg": beta, "true_airspeed_m_s": float(speed),
                            "wind_ned_m_s": None if wind_ned is None else wind_ned.tolist(),
                            "failure_active": failure is not None and failure.active(self.physics),
                            "battery_consumed_mah": None if self.physics.battery is None else self.physics.battery.consumed_ah*1000}
        if control is not commanded:
            self.last_record["commanded_control"] = asdict(commanded)
        if params.pusher_installed and params.pusher_system is not None:
            from ..propulsion import pusher_performance
            pusher = pusher_performance(control, params, float(v_air[0]))
            self.last_record["pusher"] = {
                "thrust_n": pusher["thrust_n"], "current_a": pusher["current_a"],
                "sustainable_current_a": params.pusher_system.sustainable_current_a,
                "above_sustainable_current": bool(pusher.get("above_sustainable_current", False))}
        if self.physics.ground is not None:
            self.last_record["ground"] = self.physics.ground.loads(self.physics.state)[2]
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
