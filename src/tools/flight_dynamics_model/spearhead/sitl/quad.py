"""Official ArduPlane 4.7.1 QuadPlane hover with four independent KDE motors."""

import argparse
from dataclasses import asdict
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
from .ground_control import ModeCommands, QUADPLANE_MODES
from .transition import transition_parameters, transition_command, analyse_transition, forward_trim_parameters
from ..propulsion import load_static_performance, lift_propulsion

from .hover_doublets import HOVER_DOUBLETS, hover_doublet_plan, hover_doublet_command, analyse_hover_doublet

from . import SITL_HOME, ardupilot_density_ratio, sitl_home_argument

GROUND_HOME = SITL_HOME

PARAMETER_RENAMES = {"Q_P_JERK_NE": "Q_P_NE_JERK", "Q_P_JERK_D": "Q_P_D_JERK"}


def hover_thrust_setting(actuator, spin_min=.15, spin_max=.95, expo=.75, air_density_ratio=1.):
    """Inverse 4.7.1 motor linearization at fixed voltage, not PWM percentage.

    4.7.1 divides VTOL throttle by the AHRS air-density ratio before linearization, so
    the parameter that reproduces a physical actuator command is thrust*ratio.
    """
    ratio = (actuator-spin_min)/(spin_max-spin_min)
    if not 0 < ratio < 1:
        raise ValueError("Hover command outside usable motor range")
    return ((1-expo)*ratio + expo*ratio**2)*air_density_ratio


def quad_parameters(collective, baseline, air_density_ratio=1.):
    # Preserve supplied log's VTOL controller gains; do not copy board/CAN wiring.
    logged = json.loads(Path(baseline).read_text())["parameters"]
    params = {PARAMETER_RENAMES.get(key, key): value for key, value in logged.items()
              if key.startswith(("Q_A_", "Q_P_")) and isinstance(value, (int, float))}
    params.update(Q_ENABLE=1, Q_FRAME_CLASS=1, Q_FRAME_TYPE=1,
                  SERVO5_FUNCTION=33, SERVO6_FUNCTION=34, SERVO7_FUNCTION=35, SERVO8_FUNCTION=36,
                  Q_M_PWM_MIN=1100, Q_M_PWM_MAX=1940, Q_M_SPIN_MIN=.15, Q_M_SPIN_MAX=.95,
                  Q_M_SPIN_ARM=.1, Q_M_THST_EXPO=.75, Q_M_HOVER_LEARN=0,
                  Q_M_THST_HOVER=hover_thrust_setting(collective, air_density_ratio=air_density_ratio), Q_THROTTLE_EXPO=0,
                  Q_M_BAT_VOLT_MAX=0, Q_M_BAT_VOLT_MIN=0,
                  ARMING_REQUIRE=1, ARMING_SKIPCHK=0, AHRS_EKF_TYPE=3, EK3_ENABLE=1,
                  FLTMODE_CH=0, RC1_DZ=0, RC2_DZ=0, RC4_DZ=0, RC3_DZ=60,
                  RC3_MIN=1000, RC3_MAX=2000, RC3_TRIM=1500, RC3_REVERSED=0,
                  ARSPD_TYPE=0, ARSPD_USE=0, Q_ASSIST_SPEED=-1, Q_VFWD_GAIN=0,
                  SCHED_LOOP_RATE=100, PTCH_TRIM_DEG=0)
    return params


def quad_throttle_pwm(mode, hover_thrust, *, minimum=1000, maximum=2000, dead_zone=60):
    """Keep neutral climb input in altitude modes, calculated throttle in manual Q modes.

    Fixed-wing modes retain centre throttle input; PT1.5 now has an estimated
    electric pusher. This helper does not validate a transition or tune throttle.
    This is a pilot input only. Firmware remains the owner of all motor outputs.
    """
    # 4.7.1 RC_Channel::pwm_to_range_dz removes the low-end deadzone.
    # Q_THROTTLE_EXPO=0, so get_pilot_throttle returns that normalized range.
    low = minimum + dead_zone
    return round(low+(maximum-low)*hover_thrust) if mode in (17, 23) else 1500


def ground_home_matches(message):
    """Check HOME_POSITION against the ground datum, allowing float32 lat/lon rounding."""
    return (abs(message.get("latitude", 0)/1e7-GROUND_HOME[0]) < 2e-6
            and abs(message.get("longitude", 0)/1e7-GROUND_HOME[1]) < 2e-6
            and abs(message.get("altitude", 0)-GROUND_HOME[2]*1000) < 10)


def analyse_hover(path, initial_altitude, duration):
    # Keep only relevant physics fields; raw sensor records remain on disk.
    rows = []
    with Path(path).open() as stream:
        for line in stream:
            row = json.loads(line)
            if not row["hold"]:
                rows.append((row["time_s"], row["state"], row["pwm"][4:8]))
    if len(rows) < 2:
        raise ValueError("No free hover physics frames")
    t = np.array([r[0] for r in rows]); t -= t[0]
    x = np.array([r[1] for r in rows]); pwm = np.array([r[2] for r in rows])
    settled = t >= min(10, duration/3)
    altitude = -x[:, 2]
    tilt = np.rad2deg(x[:, 6:8])
    rms = np.sqrt(np.mean(tilt[settled]**2, axis=0))
    altitude_error = float(max(abs(altitude[settled]-initial_altitude)))
    checks = {"duration": bool(t[-1] >= duration-.02),
              "settled_tilt": bool(max(rms) < 2), "altitude_hold": altitude_error < 2,
              "attitude_envelope": bool(np.max(abs(tilt)) < 15),
              "motor_range": bool(np.all((pwm >= 1100) & (pwm <= 1940)))}
    return {"passed": all(checks.values()), "checks": checks, "test_kind": "quad_hover",
            "free_duration_s": float(t[-1]), "attitude_rms_error_deg": rms.tolist(),
            "max_tilt_deg": float(np.max(abs(tilt))), "altitude_range_m": [float(min(altitude)), float(max(altitude))],
            "settled_max_altitude_error_m": altitude_error,
            "horizontal_displacement_m": float(np.linalg.norm(x[-1,:2]-x[0,:2])),
            "lift_pwm_range_us": [int(pwm.min()), int(pwm.max())],
            "final_lift_pwm_us": pwm[-1].tolist()}


def run(args):
    from pymavlink import mavutil
    interactive = getattr(args, "interactive", False)
    transition = getattr(args, "transition", False)
    doublet = getattr(args, "doublet", None)
    if doublet and (interactive or transition or args.duration < 45):
        raise ValueError("Hover doublets require 45 s and cannot be interactive or transitions")
    parameter_set = getattr(args, "transition_parameters", "first_transition_2026_10_09")
    if interactive and transition:
        raise ValueError("Choose an interactive session or an automated transition")
    default = "quad_transition.yaml" if transition else "quad_hover.yaml"
    physics, mapping, manifest = load_setup(getattr(args, "config", None) or CONFIG_ROOT / default)
    if manifest["config"].get("flight_regime") != "hover":
        raise ValueError("Quad runner requires a hover initial point")
    if transition and (physics.params.aircraft_id != "pt15" or -physics.state[2] < 40):
        raise ValueError("Transition diagnostic requires PT1.5 and at least 40 m initial altitude")
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=False)
    firmware = verify_firmware(args.binary)
    collective = manifest["trim_control"]["lift_throttle"]
    baseline = CONFIG_ROOT / "quad_hover_controller.json"
    params = quad_parameters(collective, baseline, ardupilot_density_ratio(-physics.state[2]))
    plan = hover_doublet_plan(doublet, params) if doublet else None
    forward_params = {}
    if transition:
        params.update(transition_parameters(parameter_set))
        forward_params, forward_state, forward_control = forward_trim_parameters(
            physics.params, mapping, -physics.state[2], include_simulation_gains=False)
        params.update(forward_params)
    params.update(SERVO3_FUNCTION=70, SERVO3_MIN=mapping.pusher.minimum,
                  SERVO3_MAX=mapping.pusher.maximum, SERVO3_TRIM=mapping.pusher.minimum)
    overlay = output / "quad.parm"
    overlay.write_text((CONFIG_ROOT / "plane_4_7_1.parm").read_text()+"\n"+
                       "\n".join(f"{key} {value}" for key, value in params.items())+"\n")
    manifest.update(firmware=firmware, test_kind="interactive_session" if interactive else "quad_hover", interactive=interactive,
                    scenario_options={"duration_s": args.duration},
                    requires_controls_ready=interactive, commandable_modes=QUADPLANE_MODES,
                    ground_home_lat_lon_alt_m=list(GROUND_HOME),
                    parameter_overlay=params, parameter_sha256=sha256(overlay),
                    controller_log_path=str(baseline), controller_log_sha256=sha256(baseline),
                    parameter_renames=PARAMETER_RENAMES,
                    forced_arm=False, hover_collective_pwm_fraction=collective,
                    hover_firmware_thrust_fraction=params["Q_M_THST_HOVER"],
                    bootstrap="stationary airborne hold; normal arm in QSTABILIZE; spool at hover; explicit release; QHOVER",
                    fdm_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=CONFIG_ROOT, text=True).strip())
    if plan:
        manifest.update(test_kind="control_doublet", doublet_plan=plan)
    if transition:
        manifest.update(test_kind="quad_transition", transition_plan={
            "hover_settle_s": 10, "forward_mode": "FBWA", "rc3_ramp_s": 5,
            "rc3_target_us": 1900, "parameter_set": parameter_set,
            "forward_gains": "frozen 9 October proposal and factory rate gains; no tuning",
            "cg_policy": "retain configured prototype CG; no diagnostic trim relocation"})
        manifest["forward_trim"] = {"speed_m_s": 24, "state": forward_state.tolist(),
                                    "control": asdict(forward_control), "firmware_parameters": forward_params}
    curve = load_static_performance(physics.params.lift_system.performance_path)
    manifest["static_lift_to_weight"] = 4*curve.lookup(44.4, 1)["thrust_n"]/(physics.params.mass*physics.params.g)
    manifest["limitations"] += ["No ground/contact, battery sag on thrust, rotor transient or propwash model; air-start only",
                                 "Logged VTOL gains are a simulation starting point, not a PT1.5 hardware tune"]
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    session = Session(physics, mapping)
    stop, errors = threading.Event(), []
    telemetry, observed, events = {}, {}, []
    link = process = None
    armed_at = release_at = None
    forward_confirmed_at = None
    version = arm_ack = controls_ready = False
    home_requested = home_confirmed = False
    last_rc = last_arm = last_mode = -100.
    critical = {key: params[key] for key in ("Q_ENABLE", "Q_FRAME_CLASS", "Q_FRAME_TYPE", "Q_M_PWM_MIN", "Q_M_PWM_MAX", "Q_M_THST_HOVER", "ARMING_REQUIRE", "ARMING_SKIPCHK", "AHRS_EKF_TYPE", "RC3_MIN", "RC3_MAX", "RC3_REVERSED", "RC3_DZ", "Q_THROTTLE_EXPO")}
    critical.update({key: value for key, value in params.items() if key.startswith(("Q_A_", "Q_P_"))})
    critical.update(SCHED_LOOP_RATE=100, Q_ASSIST_SPEED=params["Q_ASSIST_SPEED"], SERVO5_FUNCTION=33, SERVO6_FUNCTION=34, SERVO7_FUNCTION=35, SERVO8_FUNCTION=36)
    critical.update({key: params[key] for key in ("SERVO3_FUNCTION", "SERVO3_MIN", "SERVO3_MAX", "SERVO3_TRIM")})
    if plan:
        critical.update({f"RC{channel}_{field}": value for channel in (1, 2, 4)
                         for field, value in (("MIN", 1000), ("MAX", 2000), ("TRIM", 1500), ("REVERSED", 0), ("DZ", 0))})
    if transition:
        critical.update(transition_parameters(parameter_set))
        critical.update(forward_params)
    report = {"passed": False, "scenario": "quad_interactive" if interactive else "quad_hover", "test_kind": "interactive_session" if interactive else "quad_hover",
              "interactive": interactive, "validation_run": not interactive, "forced_arm": False}
    if plan:
        report.update(scenario="quad_"+doublet, test_kind="control_doublet")
    if transition:
        report.update(scenario="quad_transition", test_kind="quad_transition")
    started = time.monotonic()

    def event(name, **data):
        row = {"time_s": physics.time, "event": name, **data}
        events.append(row); print(json.dumps(row), flush=True)

    def command_event(**data):
        with (output / "control-events.jsonl").open("a") as stream:
            stream.write(json.dumps({"time_s": physics.time, **data}) + "\n")
        event("ground_control", **data)

    mode_commands = ModeCommands(output / "commands", command_event, QUADPLANE_MODES)

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
                command = [str(Path(args.binary).resolve()), "--model", "JSON:127.0.0.1",
                           "--home", sitl_home_argument(), "--rate", "200", "--speedup", "1",
                           "--sim-port-out", str(args.port), "--base-port", str(args.mavlink_port),
                           "--defaults", str(overlay), "--wipe"]
                manifest["sitl_command"] = command
                (output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
                process = subprocess.Popen(command, cwd=output, stdout=slog, stderr=subprocess.STDOUT)
                while time.monotonic()-started < (getattr(args, "wall_timeout", None) or max(150, args.duration*5)):
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
                    if msg is not None and msg.get_type() != "BAD_DATA":
                        kind = msg.get_type(); telemetry[kind] = msg.to_dict()
                        tlog.write(json.dumps({"fdm_time_s": physics.time, **telemetry_record(msg.to_dict())})+"\n")
                        if kind == "HEARTBEAT":
                            link.target_system, link.target_component = msg.get_srcSystem(), msg.get_srcComponent()
                            link.mav.heartbeat_send(6,8,0,0,4)
                            link.mav.request_data_stream_send(link.target_system, link.target_component, 0, 10, 1)
                            if transition and release_at is not None and physics.time-release_at >= 10:
                                if msg.custom_mode == 5 and forward_confirmed_at is None:
                                    forward_confirmed_at = physics.time-release_at
                                    event("forward_mode_confirmed", mode=5,
                                          free_time_s=forward_confirmed_at)
                                elif msg.custom_mode != 5:
                                    # A fresh confirmation starts a fresh ramp after mode loss.
                                    forward_confirmed_at = None
                            if not version:
                                link.mav.command_long_send(link.target_system, link.target_component, 512, 0, 148,0,0,0,0,0,0)
                            for key in critical:
                                if key not in observed:
                                    link.mav.param_request_read_send(link.target_system, link.target_component, key.encode(), -1)
                            if home_requested and (not home_confirmed or release_at is None):
                                link.mav.command_long_send(link.target_system, link.target_component, 512, 0, 242,0,0,0,0,0,0)
                            if msg.base_mode & 128 and armed_at is None:
                                armed_at = physics.time; event("armed", mode=msg.custom_mode)
                            if armed_at is not None and not msg.base_mode & 128:
                                raise RuntimeError("Firmware disarmed")
                        elif kind == "AUTOPILOT_VERSION":
                            if msg.flight_sw_version != 0x040701ff:
                                raise RuntimeError("Unexpected firmware version")
                            version = True
                        elif kind == "PARAM_VALUE" and msg.param_id in critical:
                            observed[msg.param_id] = msg.param_value
                        elif kind == "HOME_POSITION" and home_requested:
                            home_confirmed = ground_home_matches(msg.to_dict())
                            if release_at is not None and not home_confirmed:
                                raise RuntimeError("Firmware ground home changed unexpectedly")
                        elif kind == "COMMAND_ACK" and msg.command == 400:
                            arm_ack = msg.result == 0; event("arm_ack", result=msg.result)
                        elif kind == "STATUSTEXT":
                            event("firmware_message", severity=msg.severity, text=msg.text)
                        elif kind == "EKF_STATUS_REPORT" and release_at is not None:
                            if msg.flags & 63 != 63 or msg.flags & (1024|32768):
                                raise RuntimeError("EKF lost navigation validity")
                    if "HEARTBEAT" not in telemetry:
                        continue
                    now = physics.time
                    if now > 22 and not home_requested:
                        link.mav.command_int_send(link.target_system, link.target_component, 0, 179, 0, 0,
                                                  0,0,0,0,round(GROUND_HOME[0]*1e7),round(GROUND_HOME[1]*1e7),GROUND_HOME[2])
                        home_requested = True
                        event("ground_home_requested", latitude=GROUND_HOME[0], longitude=GROUND_HOME[1], altitude_m=GROUND_HOME[2])
                    desired_mode = 18 if release_at is not None and now-release_at >= 3 else 17
                    if transition and release_at is not None and now-release_at >= 3:
                        desired_mode, _ = transition_command(now-release_at, telemetry["HEARTBEAT"]["custom_mode"], forward_confirmed_at)
                    if plan and release_at is not None and now-release_at >= 3:
                        desired_mode, _ = hover_doublet_command(plan, now-release_at, telemetry["HEARTBEAT"]["custom_mode"])
                    if (interactive and not controls_ready and release_at is not None
                            and now-release_at >= 3 and telemetry["HEARTBEAT"]["custom_mode"] == 18):
                        controls_ready = True
                        (output / "controls-ready.json").write_text(json.dumps({"time_s": now})+"\n")
                        event("mode_controls_ready")
                    if interactive:
                        mode_commands.update(link, None if msg is None else msg.to_dict(), ready=controls_ready)
                    if (not interactive or not controls_ready) and now-last_mode >= 1 and telemetry["HEARTBEAT"]["custom_mode"] != desired_mode:
                        link.mav.set_mode_send(link.target_system,1,desired_mode); last_mode = now
                        event("mode_requested", mode=desired_mode)
                    if now >= 25 and armed_at is None and now-last_arm >= 3:
                        if not version or not parameters_match(observed, critical):
                            raise RuntimeError(f"Unverified QuadPlane parameters: {observed}")
                        if not home_confirmed:
                            raise RuntimeError("Explicit ground home not confirmed before arming")
                        link.mav.command_long_send(link.target_system,link.target_component,400,0,1,0,0,0,0,0,0)
                        last_arm = now; event("normal_arm_requested")
                    if now > 55 and armed_at is None:
                        raise RuntimeError("Normal QuadPlane arming did not succeed")
                    rc = [1500,1500,1000,1500,1500,1500,1500,1500]
                    if armed_at is not None:
                        rc[2] = quad_throttle_pwm(telemetry["HEARTBEAT"]["custom_mode"] if controls_ready else desired_mode, params["Q_M_THST_HOVER"],
                                                  minimum=params["RC3_MIN"], maximum=params["RC3_MAX"], dead_zone=params["RC3_DZ"])
                    if transition and release_at is not None and now-release_at >= 3:
                        _, rc[2] = transition_command(now-release_at, telemetry["HEARTBEAT"]["custom_mode"], forward_confirmed_at)
                    if plan and release_at is not None and now-release_at >= 3:
                        _, rc = hover_doublet_command(plan, now-release_at, telemetry["HEARTBEAT"]["custom_mode"])
                    if now-last_rc >= .05:
                        link.mav.rc_channels_override_send(link.target_system,link.target_component,*rc); last_rc = now
                    if armed_at is not None and release_at is None and now-armed_at >= 4:
                        if not arm_ack or not home_confirmed or telemetry["HEARTBEAT"]["custom_mode"] != 17 or telemetry.get("EKF_STATUS_REPORT",{}).get("flags",0)&63 != 63:
                            continue
                        force, moment, _ = lift_propulsion(mapping.control(session.last_frame.pwm), physics.params)
                        if abs(-force[2]/(physics.params.mass*physics.params.g)-1) > .05:
                            if now-armed_at > 15:
                                raise RuntimeError("Quad motor outputs did not reach static hover support")
                            continue
                        session.release(); release_at = now
                        event("release", pwm=list(session.last_frame.pwm), lift_force_n=force.tolist(), lift_moment_nm=moment.tolist())
                    if release_at is not None:
                        x = physics.state.copy()
                        if np.max(abs(np.rad2deg(x[6:8]))) > (60 if interactive or transition else 20) or not 5 < -x[2] < (120 if interactive or transition or plan else 40):
                            raise RuntimeError("QuadPlane diagnostic envelope exceeded (attitude or altitude; ground contact is not modelled)")
                        if transition and np.linalg.norm(x[3:6]) > 30:
                            raise RuntimeError("Transition diagnostic exceeded the 30 m/s test-speed limit")
                        if now-release_at >= args.duration:
                            break
                else:
                    raise TimeoutError("Quad hover wall-clock deadline exceeded")
            stop.set(); worker.join(timeout=3)
            if worker.is_alive() or errors:
                raise RuntimeError("Physics log did not close cleanly")
            report.update({"passed": True, "test_kind": "interactive_session", "validation_run": False,
                           "free_duration_s": physics.time-release_at} if interactive else
                          {} if transition or plan else analyse_hover(output / "truth.jsonl", -manifest["initial_state"][2], args.duration))
            if not interactive and not transition and not plan and telemetry.get("HEARTBEAT",{}).get("custom_mode") != 18:
                report["passed"] = False; report["mode_error"] = "QHOVER not confirmed"
        except Exception as exc:
            report["error"] = repr(exc); event("failed", error=repr(exc))
        finally:
            stop.set()
            if process is not None and process.poll() is None:
                process.terminate()
                try: process.wait(timeout=5)
                except subprocess.TimeoutExpired: process.kill(); process.wait()
            worker.join(timeout=3)
            if link is not None: link.close()
            report.update(normal_arm=arm_ack, ground_home_confirmed=home_confirmed, observed_parameters=observed, release_at_s=release_at,
                          final_mode=telemetry.get("HEARTBEAT",{}).get("custom_mode"), errors=errors,
                          wall_s=time.monotonic()-started, simulated_s=physics.time)
            if plan:
                report.update(analyse_hover_doublet(output / "truth.jsonl", release_at, args.duration, manifest["initial_state"], plan))
                if report.get("error") or errors:
                    report["passed"] = False
            if transition:
                idle_pwm = int(np.ceil(params["Q_M_PWM_MIN"] + params["Q_M_SPIN_ARM"] *
                                       (params["Q_M_PWM_MAX"] - params["Q_M_PWM_MIN"]))) + 1
                analysis = analyse_transition(output / "truth.jsonl", release_at, args.duration,
                                              -manifest["initial_state"][2], report["final_mode"], idle_pwm)
                report.update(analysis)
                report["firmware_transition_done"] = any(e.get("text") == "Transition done" for e in events)
                if report.get("error") or errors:
                    report.update(passed=False, transition_completed=False)
            (output / "events.json").write_text(json.dumps(events,indent=2)+"\n")
            (output / "report.json").write_text(json.dumps(report,indent=2)+"\n")
            print(json.dumps(report,indent=2),flush=True)
    return report["passed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--doublet", choices=tuple(HOVER_DOUBLETS), help="Hover pilot-input doublet; frozen VTOL gains")
    parser.add_argument("--interactive", action="store_true", help="Enable Flight Deck mode commands after QHOVER bootstrap")
    parser.add_argument("--transition", action="store_true", help="PT1.5 hover-to-FBWA diagnostic; no CG substitution")
    parser.add_argument("--transition-parameters", choices=("first_transition_2026_10_09",), default="first_transition_2026_10_09",
                        help="Frozen transition proposal; this option does not tune gains")
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config", help="Validated aircraft and initial-condition configuration")
    parser.add_argument("--duration", type=float, default=60)
    parser.add_argument("--wall-timeout", type=float, help="Host runtime allowance; leaves physics and controller timing unchanged")
    parser.add_argument("--port", type=int, default=19602)
    parser.add_argument("--mavlink-port", type=int, default=16360)
    args = parser.parse_args()
    if not np.isfinite(args.duration) or args.duration < 15:
        parser.error("hover duration must be finite and at least 15 seconds")
    if args.wall_timeout is not None and (not np.isfinite(args.wall_timeout) or args.wall_timeout <= 0):
        parser.error("wall timeout must be finite and positive")
    if args.transition and (args.interactive or args.duration < 60):
        parser.error("transition requires at least 60 seconds and cannot be interactive")
    if args.doublet and (args.interactive or args.transition or args.duration < 45):
        parser.error("hover doublets require at least 45 seconds and cannot be interactive or transitions")
    raise SystemExit(0 if run(args) else 1)


if __name__ == "__main__":
    main()
