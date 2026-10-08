"""CLI: validate, serve a backend, or launch a pinned real SITL smoke test."""

import argparse
from collections import Counter
import json
from pathlib import Path
import socket
import subprocess
import threading
import time

import numpy as np

from .backend import CONFIG_ROOT, Session, load_setup, serve, sha256, verify_firmware


def telemetry_record(data):
    """MAVLink uses NaN for unavailable optional values; preserve that as null."""
    unavailable = []

    def convert(value, path):
        if isinstance(value, dict):
            return {k: convert(v, f"{path}.{k}") for k, v in value.items()}
        if isinstance(value, (list, tuple, bytes, bytearray)):
            return [convert(v, f"{path}[{i}]") for i, v in enumerate(value)]
        if isinstance(value, float) and not np.isfinite(value):
            unavailable.append(path)
            return None
        return value

    return {"message": convert(data, "message"), "unavailable_fields": unavailable}


def run(args):
    physics, mapping, manifest = load_setup(args.config)
    if args.command == "validate":
        print(json.dumps(manifest, indent=2, allow_nan=False))
        return
    if args.release_file and Path(args.release_file).exists():
        raise ValueError("Release file already exists; choose a new run-specific path")
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest["test_kind"] = "bootstrap_protocol" if args.command == "smoke" else "backend"
    manifest["bootstrap_hold"] = True
    manifest["fdm_commit"] = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=CONFIG_ROOT, text=True).strip()
    manifest["firmware"] = verify_firmware(args.binary)
    if args.command == "smoke":
        manifest["parameter_sha256"] = sha256(CONFIG_ROOT / "plane_4_7_1.parm")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    stop = threading.Event()
    session = Session(physics, mapping)
    errors = []
    process = None
    link = None
    worker = None
    with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
        sock.bind(("127.0.0.1", args.port))
        if args.command == "serve":
            serve(sock, session, output / "truth.jsonl", stop, duration_s=args.duration,
                  release_file=args.release_file)
            return

        from pymavlink import mavutil

        def backend():
            try:
                serve(sock, session, output / "truth.jsonl", stop, duration_s=args.duration)
            except Exception as exc:
                errors.append(repr(exc))

        worker = threading.Thread(target=backend, daemon=True)
        worker.start()
        started = time.monotonic()
        messages = Counter()
        version = None
        last_messages = {}
        try:
            with (output / "sitl.log").open("x") as sitl_log, (output / "telemetry.jsonl").open("x") as telemetry:
                command = [str(Path(args.binary).resolve()), "--model", "JSON:127.0.0.1",
                           "--home", "47.397742,8.545594,500,0", "--speedup", "1", "--rate", "200",
                           "--sim-port-out", str(args.port), "--base-port", str(args.mavlink_port),
                           "--defaults", str(CONFIG_ROOT / "plane_4_7_1.parm"), "--wipe"]
                process = subprocess.Popen(command, cwd=output, stdout=sitl_log, stderr=subprocess.STDOUT)
                manifest["sitl_command"] = command
                (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
                deadline = started + max(90., args.duration * 4)
                while time.monotonic() < deadline and worker.is_alive():
                    if process.poll() is not None:
                        raise RuntimeError(f"SITL exited with {process.returncode}; inspect sitl.log")
                    if link is None:
                        try:
                            link = mavutil.mavlink_connection(f"tcp:127.0.0.1:{args.mavlink_port}",
                                                             source_system=255, retries=1)
                        except OSError:
                            stop.wait(.2)
                            continue
                    msg = link.recv_match(blocking=True, timeout=.2)
                    if msg is None:
                        continue
                    kind = msg.get_type()
                    if kind == "BAD_DATA":
                        continue
                    messages[kind] += 1
                    data = msg.to_dict()
                    # Byte arrays in custom firmware identity need JSON conversion.
                    telemetry.write(json.dumps({"fdm_time_s": physics.time, **telemetry_record(data)},
                                               default=list, allow_nan=False) + "\n")
                    last_messages[kind] = data
                    if kind == "HEARTBEAT":
                        link.mav.heartbeat_send(6, 8, 0, 0, 4)
                        if messages[kind] == 1:
                            link.mav.request_data_stream_send(msg.get_srcSystem(), msg.get_srcComponent(),
                                                              mavutil.mavlink.MAV_DATA_STREAM_ALL, 10, 1)
                        if version is None:
                            link.mav.command_long_send(msg.get_srcSystem(), msg.get_srcComponent(),
                                mavutil.mavlink.MAV_CMD_REQUEST_MESSAGE, 0,
                                mavutil.mavlink.MAVLINK_MSG_ID_AUTOPILOT_VERSION, 0, 0, 0, 0, 0, 0)
                    if kind == "AUTOPILOT_VERSION":
                        version = msg.flight_sw_version
                if worker.is_alive():
                    raise TimeoutError("SITL smoke did not complete within wall-clock deadline")
                if errors:
                    raise RuntimeError("; ".join(errors))
                if physics.time < args.duration:
                    raise RuntimeError("Incomplete sensor exchange")
                if version is None or (version >> 8) != 0x040701 or (version & 255) != 255:
                    raise RuntimeError(f"MAVLink firmware version mismatch: {version}")
                required = ["HEARTBEAT", "ATTITUDE", "GLOBAL_POSITION_INT", "SERVO_OUTPUT_RAW"]
                if any(messages[k] == 0 for k in required):
                    raise RuntimeError(f"Missing required telemetry: {dict(messages)}")
                attitude = last_messages["ATTITUDE"]
                error = np.rad2deg(np.array([attitude[k] for k in ["roll", "pitch", "yaw"]]) - physics.state[6:9])
                error = (error + 180) % 360 - 180
                if not np.isfinite(error).all() or max(abs(error)) > 2:
                    raise RuntimeError(f"Attitude estimate/truth mismatch: {error}")
                servo = last_messages["SERVO_OUTPUT_RAW"]
                reported_pwm = [servo[f"servo{i}_raw"] for i in range(1, 17)]
                if reported_pwm != list(session.last_frame.pwm[:16]):
                    raise RuntimeError("MAVLink servo outputs differ from backend PWM")
                position = last_messages["GLOBAL_POSITION_INT"]
                from ..rotations import euler_body_to_ned
                velocity_truth = euler_body_to_ned(*physics.state[6:9]) @ physics.state[3:6]
                velocity_error = np.array([position[k]/100 for k in ["vx", "vy", "vz"]]) - velocity_truth
                altitude_error = position["alt"]/1000 - (500 - physics.state[2])
                if np.linalg.norm(velocity_error) > .5 or abs(altitude_error) > 2:
                    raise RuntimeError("GPS velocity/absolute altitude and FDM truth disagree")
                report = {"passed": True, "test_kind": "bootstrap_protocol", "controlled_flight": False,
                          "simulated_s": physics.time, "wall_s": time.monotonic() - started,
                          "accepted_frames": session.accepted, "duplicate_frames": session.duplicates,
                          "message_counts": dict(messages), "flight_sw_version": version,
                          "last_attitude_error_deg": error.tolist(),
                          "last_velocity_error_m_s": velocity_error.tolist(),
                          "last_absolute_altitude_error_m": altitude_error,
                          "servo_pwm_matches": True,
                          "last_ekf_status": last_messages.get("EKF_STATUS_REPORT"),
                          "last_system_status": last_messages.get("SYS_STATUS"),
                          "errors": errors}
                (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
                print(json.dumps(report, indent=2))
        except Exception as exc:
            (output / "report.json").write_text(json.dumps({"passed": False, "error": repr(exc),
                "simulated_s": physics.time, "accepted_frames": session.accepted,
                "controlled_flight": False}, indent=2) + "\n")
            raise
        finally:
            stop.set()
            if process is not None and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            if worker is not None:
                worker.join(timeout=2)
            if link is not None:
                link.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["validate", "serve", "smoke"])
    parser.add_argument("--config")
    parser.add_argument("--output", default="sitl-run")
    parser.add_argument("--port", type=int, default=9002)
    parser.add_argument("--mavlink-port", type=int, default=5760)
    parser.add_argument("--duration", type=float, default=30.)
    parser.add_argument("--binary")
    parser.add_argument("--release-file", help="serve only: new file containing 'release' ends bootstrap hold")
    args = parser.parse_args()
    if not np.isfinite(args.duration) or args.duration <= 0:
        parser.error("duration must be positive and finite")
    if args.command != "validate" and not args.binary:
        parser.error("serve/smoke requires an explicitly selected 4.7.1 --binary")
    if args.command != "serve" and args.release_file:
        parser.error("release-file is supported only by serve; smoke always holds")
    run(args)


if __name__ == "__main__":
    main()
