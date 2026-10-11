"""Normally armed ArduPlane 4.7.1 ground takeoff, VTOL loiter and landing."""

import argparse
import json
from pathlib import Path
import socket
import subprocess
import threading
import time

import numpy as np

from .backend import CONFIG_ROOT, Session, load_setup, serve, sha256, verify_firmware
from .__main__ import telemetry_record
from .flight import parameters_match
from .quad import GROUND_HOME, quad_parameters, ground_home_matches


class GroundMission:
    """Request-driven MAVLink upload; acceptance is required before arming."""
    def __init__(self, altitude, duration):
        self.items = [(16, 0., 0.), (84, 0., altitude), (3000, 3., 0.),
                      (19, duration, altitude), (85, 0., 0.)]
        self.sent = set()
        self.accepted = False
        self.started = None

    def update(self, link, message, now):
        if self.started is None:
            link.mav.mission_count_send(link.target_system, link.target_component, len(self.items))
            self.started = now
        if now-self.started > 15 and not self.accepted:
            raise TimeoutError("VTOL mission upload timed out")
        if message is None:
            return
        kind = message.get_type()
        if kind in ("MISSION_REQUEST", "MISSION_REQUEST_INT"):
            seq = message.seq
            if not 0 <= seq < len(self.items):
                raise ValueError("Firmware requested an invalid mission item")
            command, dwell, altitude = self.items[seq]
            common = (link.target_system, link.target_component, seq)
            # Plane accepts INT items even when its legacy request message is used.
            link.mav.mission_item_int_send(*common, 6, command, 0, 1, dwell, 0, 0, 0,
                                          round(GROUND_HOME[0]*1e7), round(GROUND_HOME[1]*1e7), altitude)
            self.sent.add(seq)
        elif kind == "MISSION_ACK":
            if message.type != 0 or self.sent != set(range(len(self.items))):
                raise RuntimeError("Firmware rejected or incompletely requested the VTOL mission")
            self.accepted = True


def analyse_ground(path, target, duration, normal_arm, mission_accepted, disarmed, events):
    rows = [json.loads(line) for line in Path(path).open()]
    x = np.array([row["state"] for row in rows])
    t = np.array([row["time_s"] for row in rows])
    altitude = -x[:, 2]
    contact = np.array([row["ground"]["contact"] for row in rows])
    hover = (altitude > target-1) & (altitude < target+1) & (np.abs(x[:, 5]) < .5)
    # Bound settled hover by contiguous samples, not scattered passages in climb/descent.
    spans = np.split(np.flatnonzero(hover), np.flatnonzero(np.diff(np.flatnonzero(hover)) != 1)+1)
    hover_time = max((t[s[-1]]-t[s[0]] for s in spans if len(s)), default=0.)
    phases = {e["phase"]: e["time_s"] for e in events if e["event"] == "phase"}
    mission_hover_time = phases.get("LANDING", 0)-phases.get("HOVER", 0)
    checks = {"normal_arm": normal_arm, "mission_accepted": mission_accepted,
              "initial_ground_contact": bool(contact[0]), "actual_liftoff": bool(np.any(~contact & (altitude > 1))),
              "target_reached": bool(max(altitude) >= target-.5),
              "mission_hover_dwell": bool(mission_hover_time >= duration-.2),
              "settled_hover": bool(hover_time >= max(5, duration-3)),
              "landed_on_gear": bool(contact[-1] and abs(x[-1, 5]) < .1),
              "firmware_disarmed": disarmed, "attitude_envelope": bool(np.max(abs(np.rad2deg(x[:, 6:8]))) < 15),
              "pusher_off": all(row["control"]["throttle"] == 0 for row in rows)}
    return {"passed": all(checks.values()), "checks": checks, "hover_duration_s": float(hover_time),
            "mission_hover_duration_s": mission_hover_time,
            "altitude_range_m": [float(min(altitude)), float(max(altitude))],
            "max_tilt_deg": float(np.max(abs(np.rad2deg(x[:, 6:8])))),
            "final_speed_m_s": float(np.linalg.norm(x[-1, 3:6])),
            "final_ground_normal_force_n": rows[-1]["ground"]["normal_force_n"]}


def run(args):
    from pymavlink import mavutil
    physics, mapping, manifest = load_setup(args.config or CONFIG_ROOT / "quad_ground.yaml")
    if manifest["config"].get("flight_regime") != "ground":
        raise ValueError("Ground flight requires a ground initial point")
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=False)
    firmware = verify_firmware(args.binary)
    params = quad_parameters(manifest["trim_control"]["lift_throttle"], CONFIG_ROOT / "quad_hover_controller.json")
    params.update(SERVO3_FUNCTION=70, SERVO3_MIN=mapping.pusher.minimum,
                  SERVO3_MAX=mapping.pusher.maximum, SERVO3_TRIM=mapping.pusher.minimum,
                  Q_GUIDED_MODE=1, Q_OPTIONS=65536, Q_NAVALT_MIN=2,
                  Q_WP_SPD_UP=1, Q_WP_SPD_DN=1, Q_LAND_FINAL_SPD=.5, TERRAIN_ENABLE=0)
    overlay = output / "quad.parm"
    overlay.write_text((CONFIG_ROOT / "plane_4_7_1.parm").read_text()+"\n"+
                       "\n".join(f"{key} {value}" for key, value in params.items())+"\n")
    mission = GroundMission(manifest["config"]["altitude_m"], args.duration)
    manifest.update(firmware=firmware, test_kind="quad_ground", interactive=False,
                    scenario_options={"hover_duration_s": args.duration}, parameter_overlay=params,
                    parameter_sha256=sha256(overlay), ground_home_lat_lon_alt_m=list(GROUND_HOME),
                    mission=mission.items, forced_arm=False, bootstrap="Free ground physics from first frame; normal arm; AUTO VTOL mission",
                    limitations=manifest["limitations"]+["No battery sag, rotor transient, ground effect or measured gear calibration"])
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    session = Session(physics, mapping)
    stop, errors, events = threading.Event(), [], []
    telemetry, observed = {}, {}
    process = link = None
    armed_at = takeoff_at = None
    version = arm_ack = home_requested = home_confirmed = disarmed = False
    last_rc = last_mode = last_arm = -100.
    phase = None
    report = {"passed": False, "test_kind": "quad_ground", "scenario": "quad_ground", "forced_arm": False}
    started = time.monotonic()

    def event(name, **data):
        row = {"time_s": physics.time, "event": name, **data}
        events.append(row); print(json.dumps(row), flush=True)

    def set_phase(value):
        nonlocal phase
        if value == phase:
            return
        phase = value
        temporary = output / "phase.tmp"
        temporary.write_text(json.dumps({"phase": value, "time_s": physics.time, "flight_start_s": takeoff_at})+"\n")
        temporary.replace(output / "phase.json")
        event("phase", phase=value)

    set_phase("ON GROUND")
    with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
        sock.bind(("127.0.0.1", args.port))
        def backend():
            try:
                serve(sock, session, output / "truth.jsonl", stop, timeout_s=15)
            except Exception as exc:
                errors.append(repr(exc))
        worker = threading.Thread(target=backend, daemon=True); worker.start()
        try:
            with (output / "sitl.log").open("x") as slog, (output / "telemetry.jsonl").open("x") as tlog:
                command = [str(Path(args.binary).resolve()), "--model", "JSON:127.0.0.1", "--home", "47.397742,8.545594,500,0",
                           "--rate", "200", "--speedup", "1", "--sim-port-out", str(args.port),
                           "--base-port", str(args.mavlink_port), "--defaults", str(overlay), "--wipe"]
                manifest["sitl_command"] = command
                (output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
                process = subprocess.Popen(command, cwd=output, stdout=slog, stderr=subprocess.STDOUT)
                while time.monotonic()-started < max(240, (args.duration+3*mission.items[1][2]+70)*5):
                    if errors or not worker.is_alive():
                        raise RuntimeError("Physics backend stopped: "+str(errors))
                    if process.poll() is not None:
                        raise RuntimeError("Firmware exited")
                    if link is None:
                        try:
                            link = mavutil.mavlink_connection(f"tcp:127.0.0.1:{args.mavlink_port}", source_system=255, retries=1)
                        except OSError:
                            stop.wait(.1); continue
                    msg = link.recv_match(blocking=True, timeout=.02)
                    now = physics.time
                    if msg is not None and msg.get_type() != "BAD_DATA":
                        kind = msg.get_type(); telemetry[kind] = msg.to_dict()
                        tlog.write(json.dumps({"fdm_time_s": now, **telemetry_record(msg.to_dict())})+"\n")
                        if kind == "HEARTBEAT":
                            link.target_system, link.target_component = msg.get_srcSystem(), msg.get_srcComponent()
                            link.mav.heartbeat_send(6,8,0,0,4)
                            link.mav.request_data_stream_send(link.target_system, link.target_component, 0, 10, 1)
                            if not version:
                                link.mav.command_long_send(link.target_system,link.target_component,512,0,148,0,0,0,0,0,0)
                            for key in params:
                                if key not in observed:
                                    link.mav.param_request_read_send(link.target_system,link.target_component,key.encode(),-1)
                            if home_requested and not home_confirmed:
                                link.mav.command_long_send(link.target_system,link.target_component,512,0,242,0,0,0,0,0,0)
                            if msg.base_mode & 128 and armed_at is None:
                                armed_at = now; event("armed", mode=msg.custom_mode)
                            if armed_at is not None and not msg.base_mode & 128:
                                if phase != "LANDING" or not session.last_record["ground"]["contact"]:
                                    raise RuntimeError("Unexpected disarm before ground landing")
                                disarmed = True; set_phase("LANDED")
                                break
                        elif kind == "AUTOPILOT_VERSION":
                            if msg.flight_sw_version != 0x040701ff:
                                raise RuntimeError("Unexpected firmware version")
                            version = True
                        elif kind == "PARAM_VALUE" and msg.param_id in params:
                            observed[msg.param_id] = msg.param_value
                        elif kind == "HOME_POSITION" and home_requested:
                            home_confirmed = ground_home_matches(msg.to_dict())
                        elif kind == "COMMAND_ACK" and msg.command == 400:
                            arm_ack = msg.result == 0; event("arm_ack", result=msg.result)
                        elif kind == "STATUSTEXT":
                            event("firmware_message", severity=msg.severity, text=msg.text)
                        elif kind == "MISSION_CURRENT" and takeoff_at is not None:
                            set_phase({1: "TAKEOFF", 3: "HOVER", 4: "LANDING"}.get(msg.seq, phase))
                        elif kind == "EKF_STATUS_REPORT" and takeoff_at is not None:
                            if msg.flags & 63 != 63 or msg.flags & (1024|32768):
                                raise RuntimeError("EKF lost navigation validity during ground flight")
                    if "HEARTBEAT" not in telemetry:
                        continue
                    if now > 22 and not home_requested:
                        link.mav.command_int_send(link.target_system,link.target_component,0,179,0,0,0,0,0,0,
                                                  round(GROUND_HOME[0]*1e7),round(GROUND_HOME[1]*1e7),GROUND_HOME[2])
                        home_requested = True; event("ground_home_requested")
                    if home_confirmed and not mission.accepted:
                        mission.update(link, msg, now)
                    desired_mode = 10 if takeoff_at is not None else 17
                    if now-last_mode >= 1 and telemetry["HEARTBEAT"]["custom_mode"] != desired_mode:
                        link.mav.set_mode_send(link.target_system,1,desired_mode); last_mode = now
                        event("mode_requested", mode=desired_mode)
                    if now >= 25 and armed_at is None and now-last_arm >= 3:
                        if (version and home_confirmed and mission.accepted and parameters_match(observed, params)
                                and telemetry.get("EKF_STATUS_REPORT",{}).get("flags",0)&63 == 63
                                and telemetry.get("AHRS2") is not None):
                            link.mav.command_long_send(link.target_system,link.target_component,400,0,1,0,0,0,0,0,0)
                            last_arm = now; event("normal_arm_requested")
                    if now > 65 and armed_at is None:
                        raise RuntimeError(f"Normal ground arming failed; missing parameters: {set(params)-set(observed)}")
                    if armed_at is not None and takeoff_at is None and arm_ack and now-armed_at >= 2:
                        if telemetry.get("EKF_STATUS_REPORT",{}).get("flags",0)&63 == 63:
                            link.mav.mission_set_current_send(link.target_system,link.target_component,1)
                            takeoff_at = now; set_phase("TAKEOFF")
                    if now-last_rc >= .05:
                        link.mav.rc_channels_override_send(link.target_system,link.target_component,
                                                           1500,1500,1500 if takeoff_at is not None else 1000,1500,1500,1500,1500,1500)
                        last_rc = now
                    if takeoff_at is not None:
                        state = physics.state.copy()
                        if np.max(abs(np.rad2deg(state[6:8]))) > 20 or not -.05 < -state[2] < mission.items[1][2]+5:
                            raise RuntimeError("Ground flight attitude/altitude envelope exceeded")
                        if now-takeoff_at > args.duration+3*mission.items[1][2]+65:
                            raise TimeoutError("Ground mission did not finish within flight deadline")
                else:
                    raise TimeoutError("Ground mission wall deadline exceeded")
            stop.set(); worker.join(timeout=3)
            if worker.is_alive() or errors:
                raise RuntimeError("Physics log did not close cleanly")
            report.update(analyse_ground(output / "truth.jsonl", mission.items[1][2], args.duration, arm_ack, mission.accepted, disarmed, events))
        except Exception as exc:
            report["error"] = repr(exc); event("failed", error=repr(exc)); set_phase("FAILED")
        finally:
            stop.set()
            if process is not None and process.poll() is None:
                process.terminate()
                try: process.wait(timeout=5)
                except subprocess.TimeoutExpired: process.kill(); process.wait()
            worker.join(timeout=3)
            if link is not None: link.close()
            report.update(normal_arm=arm_ack, mission_accepted=mission.accepted, firmware_disarmed=disarmed,
                          observed_parameters=observed, ground_home_confirmed=home_confirmed,
                          flight_start_s=takeoff_at, wall_s=time.monotonic()-started, simulated_s=physics.time,
                          errors=errors, final_mode=telemetry.get("HEARTBEAT",{}).get("custom_mode"))
            (output / "events.json").write_text(json.dumps(events,indent=2)+"\n")
            (output / "report.json").write_text(json.dumps(report,indent=2)+"\n")
            print(json.dumps(report,indent=2),flush=True)
    return report["passed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config")
    parser.add_argument("--duration", type=float, default=60, help="Hover dwell; climb/landing are additional")
    parser.add_argument("--port", type=int, default=19602)
    parser.add_argument("--mavlink-port", type=int, default=16360)
    args = parser.parse_args()
    if not np.isfinite(args.duration) or not 15 <= args.duration <= 600:
        parser.error("Hover dwell must be 15–600 seconds")
    if args.duration != int(args.duration):
        parser.error("Mission hover dwell must be a whole number of seconds")
    raise SystemExit(0 if run(args) else 1)


if __name__ == "__main__":
    main()
