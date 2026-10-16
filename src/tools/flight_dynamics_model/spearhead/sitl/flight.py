"""Real 4.7.1 arming, surface direction and hold-free fixed-wing scenarios."""

import argparse
import json
from pathlib import Path
import socket
import subprocess
import threading
import time

import numpy as np

from ..controls import Control
from ..dynamics import dynamics
from . import SITL_HOME, ardupilot_density_ratio, sitl_home_argument
from .backend import CONFIG_ROOT, Session, default_wall_timeout, load_setup, serve, sha256, verify_firmware
from .__main__ import telemetry_record
from .metrics import analyse_truth, analyse_elevator_truth, command_at, elevator_at
from .ground_control import ModeCommands
from .doublets import DOUBLETS, doublet_plan, doublet_delta, doublet_mode, analyse_doublet
from .back_transition import back_transition_command, analyse_back_transition, forward_release_ready, HOVER_MODE


def surface_pwm(channel, angle):
    normalized = np.rad2deg(angle) / channel.travel_deg * (-1 if channel.reverse else 1)
    if not -1 <= normalized <= 1:
        raise ValueError("Trim outside configured surface travel")
    span = channel.maximum-channel.trim if normalized >= 0 else channel.trim-channel.minimum
    return round(channel.trim + normalized*span)


def manual_surface_rc(channel, trim_angle, firmware_trim_pwm, delta_deg):
    """Invert Plane's asymmetric servo-trim mapping to command physical degrees.

    MANUAL RC2 is centred at 1500; firmware output is centred at SERVO2_TRIM.
    Neither direct PWM injection nor an FDM control override is used.
    """
    desired_pwm = surface_pwm(channel, trim_angle+np.deg2rad(delta_deg))
    span = (channel.maximum-firmware_trim_pwm if desired_pwm >= firmware_trim_pwm
            else firmware_trim_pwm-channel.minimum)
    if span <= 0:
        raise ValueError("Elevator servo trim must leave authority in both directions")
    ideal_rc = round(1500+500*(desired_pwm-firmware_trim_pwm)/span)
    # The final PWM conversion truncates. A rounded RC inverse can miss by
    # another output count; choose the nearest realizable command locally.
    def output_pwm(rc):
        normalized = (rc-1500)/500
        output_span = (channel.maximum-firmware_trim_pwm if normalized >= 0
                       else firmware_trim_pwm-channel.minimum)
        return firmware_trim_pwm+int(normalized*output_span)
    candidates = range(max(1000,ideal_rc-2),min(2000,ideal_rc+2)+1)
    return min(candidates,key=lambda rc:(abs(output_pwm(rc)-desired_pwm),abs(rc-ideal_rc)))


def manual_elevator_rc(mapping, trim, firmware_trim_pwm, delta_deg):
    return manual_surface_rc(mapping.elevator, trim.de, firmware_trim_pwm, delta_deg)


def startup_states(initial, roll_deg, pitch_deg):
    if not np.isfinite([roll_deg, pitch_deg]).all() or max(abs(roll_deg), abs(pitch_deg)) > .5:
        raise ValueError("Initial attitude perturbations must be within +/-0.5 degrees")
    airborne = initial.copy()
    airborne[6:8] += np.deg2rad([roll_deg, pitch_deg])
    stationary = np.zeros(12)
    stationary[6:9] = airborne[6:9]
    return stationary, airborne


def bootstrap_mode(scenario, armed_at, now):
    if scenario not in ("manual", "elevator", "fbwa", "steps", "back_transition", *DOUBLETS):
        raise ValueError("Unknown flight scenario")
    if scenario == "back_transition":
        # Frozen Q_OPTIONS requires normal arming in a Q mode. Remain held,
        # then confirm FBWA before releasing the precomputed forward trim.
        return 17 if armed_at is None or now-armed_at < 3 else 5
    return 5 if scenario not in ("manual", "elevator") and armed_at is not None and now-armed_at >= 3 else 0


def parameters_match(observed, expected):
    return set(expected) <= set(observed) and all(
        np.isclose(observed[k], v, rtol=0, atol=1e-6) for k, v in expected.items())


def flight_parameters(mapping, state, trim):
    """Firmware servo trims, not hidden forces, establish the physical trim."""
    pitch_scale = 5 / mapping.elevator.travel_deg
    roll_scale = 5 / mapping.aileron.travel_deg
    return {"SERVO1_TRIM": surface_pwm(mapping.aileron, trim.da),
            "SERVO2_TRIM": surface_pwm(mapping.elevator, trim.de),
            "SERVO4_TRIM": surface_pwm(mapping.rudder, trim.dr),
            "RC1_DZ": 0, "RC2_DZ": 0, "RC4_DZ": 0,
            "RC3_TRIM": 1000, "FLTMODE_CH": 0,
            "ARMING_REQUIRE": 1, "ARMING_SKIPCHK": 0,
            "AHRS_EKF_TYPE": 3, "EK3_ENABLE": 1,
            "PTCH_TRIM_DEG": float(np.rad2deg(state[7])),
            "STAB_PITCH_DOWN": 0, "KFF_THR2PTCH": 0,
            "TRIM_THROTTLE": float(trim.throttle*100),
            "AIRSPEED_CRUISE": float(np.linalg.norm(state[3:6])),
            "SCALING_SPEED": float(np.linalg.norm(state[3:6])),
            "AIRSPEED_MIN": 15, "AIRSPEED_MAX": 35,
            # Gains are tuned for provisional +/-20 travel and scale with travel.
            # Pitch FF uses local rate damping / elevator authority; P/I and
            # roll I were rechecked against quantized-output step/recovery tests.
            "PTCH_RATE_FF": 1.8*pitch_scale, "PTCH_RATE_P": .32*pitch_scale,
            "PTCH_RATE_I": .6*pitch_scale, "PTCH2SRV_TCONST": .2,
            "RLL_RATE_FF": .345*roll_scale, "RLL_RATE_P": .08*roll_scale,
            "RLL_RATE_I": .3*roll_scale,
            "ROLL_LIMIT_DEG": 30, "PTCH_LIM_MIN_DEG": -10, "PTCH_LIM_MAX_DEG": 10}


def prototype_forward_parameters(mapping, state, trim, parameter_set="first_transition_2026_10_09"):
    """Reuse the saved forward baseline; derive only the physical trim mapping."""
    from .transition import transition_parameters
    params = flight_parameters(mapping, state, trim)
    baseline = transition_parameters(parameter_set)
    params.update({key: value for key, value in baseline.items()
                   if not key.startswith(("Q_", "FLTMODE", "FS_", "BATT_"))})
    params.update(SERVO3_MIN=mapping.pusher.minimum, SERVO3_MAX=mapping.pusher.maximum,
                  SERVO3_TRIM=mapping.pusher.minimum)
    params["TRIM_THROTTLE"] = round(params["TRIM_THROTTLE"])
    return params


def firmware_pitot_speed(telemetry, timestamps, observed, now, max_age_s=2.):
    """Use fresh firmware airspeed only with the simulated pitot verified healthy.

    Older pymavlink dialects do not decode AIRSPEED (295). VFR_HUD provides
    firmware airspeed too, but is accepted only with ARSPD_TYPE/USE confirmed
    and fresh present/enabled/healthy differential-pressure sensor status.
    """
    if not parameters_match(observed, {"ARSPD_TYPE": 100, "ARSPD_USE": 1}):
        return None
    status = telemetry.get("SYS_STATUS", {})
    if not 0 <= now-timestamps.get("SYS_STATUS", -np.inf) <= max_age_s:
        return None
    # MAV_SYS_STATUS_SENSOR_DIFFERENTIAL_PRESSURE = 16 in the bundled dialect.
    if not all(status.get(key, 0) & 16 for key in
               ("onboard_control_sensors_present", "onboard_control_sensors_enabled", "onboard_control_sensors_health")):
        return None
    for kind in ("AIRSPEED", "VFR_HUD"):
        value = telemetry.get(kind, {}).get("airspeed")
        if (value is not None and np.isfinite(value) and value >= 0 and
                0 <= now-timestamps.get(kind, -np.inf) <= max_age_s):
            return float(value)
    return None


def run(args):
    from pymavlink import mavutil
    default_config = "elevator_demo.yaml" if args.scenario == "elevator" else "fixed_wing_control.yaml"
    physics, mapping, manifest = load_setup(args.config or CONFIG_ROOT / default_config)
    if manifest["config"].get("flight_regime", "fixed_wing") != "fixed_wing":
        raise ValueError("Fixed-wing runner requires a forward-trim initial point")
    initial = physics.state.copy()
    stationary, air_start_state = startup_states(initial, args.perturb_roll, args.perturb_pitch)
    trim = Control(**manifest["trim_control"])
    firmware = verify_firmware(args.binary)
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    interactive = getattr(args, "interactive", False)
    if interactive:
        (output / "commands").mkdir()
    prototype = physics.params.aircraft_id == "pt15"
    back_transition = args.scenario == "back_transition"
    is_doublet = args.scenario in DOUBLETS
    plan = doublet_plan(args.scenario, trim, mapping) if is_doublet else None
    if back_transition and (not prototype or len(mapping.lift) != 4):
        raise ValueError("Forward-to-hover requires PT1.5 with four mapped lift motors")
    parameter_set = getattr(args, "forward_parameters", None) or "first_transition_2026_10_09"
    params = (prototype_forward_parameters(mapping, initial, trim, parameter_set) if prototype else
              flight_parameters(mapping, initial, trim))
    if back_transition:
        from .quad import quad_parameters
        from ..propulsion import hover_collective
        from .transition import transition_parameters
        forward_params = params
        params = quad_parameters(hover_collective(physics.params), CONFIG_ROOT / "quad_hover_controller.json",
                                ardupilot_density_ratio(-physics.state[2]))
        params.update(transition_parameters(parameter_set))
        params.update({key: value for key, value in forward_params.items() if not key.startswith(
            ("RC", "FLTMODE", "ARMING", "AHRS", "EK3"))})
    parameter_path = output / "flight.parm"
    parameter_path.write_text((CONFIG_ROOT / "plane_4_7_1.parm").read_text() + "\n" +
                              "\n".join(f"{k} {v}" for k, v in params.items()) + "\n")
    manifest.update(firmware=firmware, test_kind=args.scenario, interactive=interactive, parameter_overlay=params,
                    parameter_sha256=sha256(parameter_path), forced_arm=False,
                    scenario_options={"duration_s": args.duration, "instance": args.instance,
                                      "perturb_roll_deg": args.perturb_roll,
                                      "perturb_pitch_deg": args.perturb_pitch},
                    bootstrap="stationary pitot calibration, logged air-start, normal arm, explicit release",
                    bootstrap_air_start_state=air_start_state.tolist(),
                    fdm_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=CONFIG_ROOT, text=True).strip())
    if prototype:
        manifest["forward_parameter_set"] = parameter_set
        manifest["forward_trim"] = {"speed_m_s": float(np.linalg.norm(initial[3:6])),
                                    "state": initial.tolist(), "control": manifest["trim_control"]}
        manifest["limitations"].append("Forward trim air start with lift off; only fixed-wing mode commands available")
    if plan:
        manifest.update(test_kind="control_doublet", doublet_plan=plan)
    if back_transition:
        manifest.update(test_kind="back_transition", transition_plan={"direction": "forward_to_hover", "forward_settle_s": 10,
                        "hover_mode": "QLOITER", "parameter_set": parameter_set})
        manifest["limitations"][-1] = "FW trim air start, then firmware-managed QLOITER braking and hover; no tuning"
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    session = Session(physics, mapping)
    # Calibrate the simulated pitot at rest. Starting at 28 m/s would subtract
    # the flight dynamic pressure as sensor zero and corrupt airspeed scaling.
    physics.state = stationary
    stop, failures = threading.Event(), []
    process = link = None
    events, telemetry, telemetry_times = [], {}, {}
    version, armed_at, release_at, air_started_at = None, None, None, None
    airspeed_calibrated = False
    home_requested = home_confirmed = False
    arm_ack = False
    last_rc = last_arm = last_mode = -100.
    last_rc_values = None
    critical = {"ARMING_REQUIRE": 1., "ARMING_SKIPCHK": 0., "AHRS_EKF_TYPE": 3.,
                "PTCH_RATE_FF": params["PTCH_RATE_FF"], "PTCH_RATE_P": params["PTCH_RATE_P"],
                "PTCH_RATE_I": params["PTCH_RATE_I"], "PTCH2SRV_TCONST": params["PTCH2SRV_TCONST"]}
    if prototype:
        critical.update({key: value for key, value in params.items() if key.startswith(
            ("RLL_RATE_", "PTCH_RATE_", "SERVO", "ARSPD_", "SCALING_SPEED"))})
    if back_transition:
        critical.update({key: value for key, value in params.items() if key.startswith(("Q_", "RC3_", "SCHED_LOOP_RATE"))})
    observed = {}
    signs = {}
    neutral_pwm = [params["SERVO1_TRIM"], params["SERVO2_TRIM"],
                   round(mapping.pusher.minimum+(mapping.pusher.maximum-mapping.pusher.minimum)*trim.throttle),
                   params["SERVO4_TRIM"], *([0]*12)]
    neutral_control = mapping.control(neutral_pwm)
    base_accel = dynamics(0, initial, physics.params, neutral_control)[9:12]
    mode = 0 if args.scenario in ("manual", "elevator") else 5
    started = time.monotonic()

    def event(name, **data):
        row = {"time_s": physics.time, "event": name, **data}
        events.append(row)
        print(json.dumps(row), flush=True)

    def command_event(**data):
        with (output / "control-events.jsonl").open("a") as stream:
            stream.write(json.dumps({"time_s": physics.time, **data}) + "\n")
        event("ground_control", **data)

    mode_commands = ModeCommands(output / "commands", command_event)

    with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
        sock.bind(("127.0.0.1", args.port))

        def backend():
            try:
                serve(sock, session, output / "truth.jsonl", stop, timeout_s=15.)
            except Exception as exc:
                failures.append(repr(exc))

        worker = threading.Thread(target=backend, daemon=True)
        worker.start()
        report = {"passed": False, "scenario": args.scenario, "forced_arm": False,
                  "interactive": interactive, "validation_run": not interactive}
        try:
            with (output / "sitl.log").open("x") as sitl_log, (output / "telemetry.jsonl").open("x") as tlog:
                command = [str(Path(args.binary).resolve()), "--model", "JSON:127.0.0.1",
                           "--home", sitl_home_argument(), "--rate", "200",
                           "--speedup", "1", "--instance", str(args.instance), "--sim-port-out", str(args.port),
                           "--base-port", str(args.mavlink_port), "--defaults", str(parameter_path), "--wipe"]
                manifest["sitl_command"] = command
                (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
                process = subprocess.Popen(command, cwd=output, stdout=sitl_log, stderr=subprocess.STDOUT)
                while time.monotonic()-started < (getattr(args, "wall_timeout", None) or default_wall_timeout(args.duration)):
                    if failures or not worker.is_alive():
                        raise RuntimeError("Physics backend stopped: " + str(failures))
                    if process.poll() is not None:
                        raise RuntimeError(f"Firmware exited {process.returncode}")
                    if link is None:
                        try:
                            link = mavutil.mavlink_connection(f"tcp:127.0.0.1:{args.mavlink_port}",
                                                             source_system=255, retries=1)
                        except OSError:
                            stop.wait(.1)
                            continue
                    msg = link.recv_match(blocking=True, timeout=.02)
                    if msg is not None and msg.get_type() != "BAD_DATA":
                        kind = msg.get_type()
                        telemetry[kind] = msg.to_dict()
                        telemetry_times[kind] = physics.time
                        tlog.write(json.dumps({"fdm_time_s": physics.time,
                                               **telemetry_record(msg.to_dict())}) + "\n")
                        if kind == "HEARTBEAT":
                            link.target_system, link.target_component = msg.get_srcSystem(), msg.get_srcComponent()
                            link.mav.heartbeat_send(6, 8, 0, 0, 4)
                            allowed_modes = {5, HOVER_MODE} if back_transition else {0, 5} if is_doublet else {mode}
                            if release_at is not None and msg.custom_mode not in allowed_modes and not interactive:
                                raise RuntimeError("Firmware changed flight mode after release")
                            if version is None:
                                link.mav.request_data_stream_send(link.target_system, link.target_component, 0, 10, 1)
                                link.mav.command_long_send(link.target_system, link.target_component, 512, 0,
                                                          148, 0, 0, 0, 0, 0, 0)
                            for key in critical:
                                if key not in observed:
                                    link.mav.param_request_read_send(link.target_system, link.target_component, key.encode(), -1)
                            if home_requested and not home_confirmed:
                                link.mav.command_long_send(link.target_system, link.target_component, 512, 0,
                                                          242, 0, 0, 0, 0, 0, 0)
                            if msg.base_mode & 128 and armed_at is None:
                                armed_at = physics.time
                                event("armed", mode=msg.custom_mode)
                            elif armed_at is not None and not (msg.base_mode & 128):
                                raise RuntimeError("Firmware disarmed during test")
                        elif kind == "AUTOPILOT_VERSION":
                            version = msg.flight_sw_version
                            if version != 0x040701ff:
                                raise RuntimeError("Unexpected firmware version")
                        elif kind == "PARAM_VALUE" and msg.param_id in critical:
                            observed[msg.param_id] = msg.param_value
                        elif kind == "EKF_STATUS_REPORT" and release_at is not None:
                            if msg.flags & 63 != 63 or msg.flags & (1024 | 32768):
                                raise RuntimeError("EKF lost navigation validity during flight")
                        elif kind == "HOME_POSITION" and home_requested:
                            home_confirmed = abs(msg.altitude-SITL_HOME[2]*1000) < 10
                        elif kind == "COMMAND_ACK" and msg.command == 400:
                            arm_ack = msg.result == 0
                            event("arm_ack", result=msg.result)
                        elif kind == "STATUSTEXT":
                            event("firmware_message", text=msg.text)
                            if msg.text == "Airspeed 1 calibrated":
                                airspeed_calibrated = True
                    if "HEARTBEAT" not in telemetry:
                        continue
                    now = physics.time
                    if interactive:
                        mode_commands.update(link, None if msg is None else msg.to_dict(),
                                             ready=release_at is not None)
                    if airspeed_calibrated and now >= 5 and air_started_at is None:
                        # Establish the disturbed attitude before EKF alignment;
                        # do not teleport an aligned estimator without gyro motion.
                        session.air_start(air_start_state)
                        air_started_at = now
                        event("bootstrap_air_start_queued", stationary_pitot_calibration=True)
                    if now > 15 and air_started_at is None:
                        raise RuntimeError("Stationary pitot calibration not confirmed")
                    if now > 22 and not home_requested:
                        link.mav.command_long_send(link.target_system, link.target_component, 179, 0,
                                                  0, 0, 0, 0, SITL_HOME[0], SITL_HOME[1], SITL_HOME[2])
                        home_requested = True
                        event("ground_home_requested", altitude_m=SITL_HOME[2])
                    # Keep controllers reset in MANUAL during artificial hold;
                    # holding a nonzero FBWA error for 25 s would wind up I.
                    desired_mode = bootstrap_mode(args.scenario, armed_at, now)
                    if release_at is not None:
                        elapsed = now-release_at
                        if is_doublet:
                            desired_mode = doublet_mode(plan, elapsed)
                        elif back_transition:
                            desired_mode, _ = back_transition_command(elapsed, telemetry["HEARTBEAT"]["custom_mode"], trim.throttle)
                    if (not interactive or release_at is None) and now-last_mode >= 1 and telemetry["HEARTBEAT"]["custom_mode"] != desired_mode:
                        link.mav.set_mode_send(link.target_system, 1, desired_mode)
                        last_mode = now
                    rc = [1500, 1500, 1000, 1500, 1500, 1500, 1500, 1500]
                    if armed_at is not None:
                        rc[2] = round(1000+1000*trim.throttle)
                        if back_transition:
                            _, rc[2] = back_transition_command(0 if release_at is None else now-release_at,
                                                               telemetry["HEARTBEAT"]["custom_mode"], trim.throttle)
                            if release_at is None and telemetry["HEARTBEAT"]["custom_mode"] == 17:
                                from .quad import quad_throttle_pwm
                                rc[2] = quad_throttle_pwm(17, params["Q_M_THST_HOVER"])
                        if args.scenario == "manual" and release_at is None:
                            stage_time = now-armed_at-3
                            phase = int(stage_time)
                            pulses = [(0, 1), (0, -1), (1, 1), (1, -1), (3, 1), (3, -1)]
                            if 1 <= phase <= 6:
                                channel, sign = pulses[phase-1]
                                rc[channel] += sign*20
                                key = f"{channel+1}:{sign:+d}"
                                if stage_time-phase > .6 and key not in signs:
                                    actual = mapping.control(session.last_frame.pwm)
                                    acceleration = dynamics(0, initial, physics.params, actual)[9:12] - base_accel
                                    axis = {0: 0, 1: 1, 3: 2}[channel]
                                    signs[key] = {"passed": bool(acceleration[axis]*sign > 0),
                                                  "rc_channel": channel+1, "rc_delta_us": sign*20,
                                                  "pwm": list(session.last_frame.pwm),
                                                  "angular_accel_delta_deg_s2": np.rad2deg(acceleration).tolist()}
                                    event("direction_check", **signs[key])
                    if now >= 25 and armed_at is None and now-last_arm >= 3:
                        if not parameters_match(observed, critical) or version is None:
                            raise RuntimeError(f"Required firmware/arming parameters unverified: {observed}")
                        # Normal arm: param2 is zero; no magic force value, no skipped checks.
                        link.mav.command_long_send(link.target_system, link.target_component, 400, 0, 1, 0, 0, 0, 0, 0, 0)
                        last_arm = now
                        event("normal_arm_requested")
                    if now > 55 and armed_at is None:
                        raise RuntimeError("Normal arming did not succeed")
                    pre_release_s = 12 if args.scenario == "manual" else 3
                    if (armed_at is not None and release_at is None and now-armed_at >= pre_release_s
                            and telemetry["HEARTBEAT"]["custom_mode"] == mode):
                        if not arm_ack:
                            raise RuntimeError("Arming ACK not confirmed before release")
                        pitot_speed = (firmware_pitot_speed(telemetry, telemetry_times, observed, now) if prototype else
                                       telemetry.get("AIRSPEED", {}).get("airspeed"))
                        if air_started_at is None or pitot_speed is None or pitot_speed < 20:
                            raise RuntimeError("Pitot airspeed invalid before release")
                        if args.scenario == "manual" and (len(signs) != 6 or not all(r["passed"] for r in signs.values())):
                            raise RuntimeError("Manual surface direction check failed")
                        if not home_confirmed:
                            raise RuntimeError("Explicit ground home altitude not confirmed")
                        if telemetry.get("EKF_STATUS_REPORT", {}).get("flags", 0) & 63 != 63:
                            raise RuntimeError("EKF attitude/velocity/position not ready")
                        ready = not back_transition or forward_release_ready(session.last_frame.pwm, mapping, params, trim.throttle,
                                                                            any(e.get("text") == "Transition done" for e in events))
                        if not ready:
                            if now-armed_at > 25:
                                raise RuntimeError("Forward trim startup did not finish firmware handover")
                        if ready:
                            release_at = now
                            session.release()
                            event("release", pwm=list(session.last_frame.pwm), mode=mode, confirmed_pitot_speed_m_s=pitot_speed)
                    if release_at is not None and args.scenario == "steps":
                        axis, delta_pwm, _ = command_at(now-release_at)
                        if axis is not None:
                            rc[axis] += delta_pwm
                    if release_at is not None and args.scenario == "elevator":
                        rc[1] = manual_elevator_rc(mapping, trim, params["SERVO2_TRIM"],
                                                   elevator_at(now-release_at))
                    if release_at is not None and is_doublet:
                        elapsed = now-release_at
                        delta = doublet_delta(plan, elapsed)
                        if delta and telemetry["HEARTBEAT"]["custom_mode"] != 0:
                            raise RuntimeError("MANUAL mode not confirmed for the doublet input")
                        if plan["channel"] == "pusher":
                            rc[2] = round(1000+1000*(trim.throttle+delta))
                        else:
                            channel = getattr(mapping, plan["channel"])
                            rc[channel.number-1] = manual_surface_rc(channel, getattr(trim, plan["surface"]),
                                                                   params[f"SERVO{channel.number}_TRIM"], delta)
                    if now-last_rc >= .05:
                        link.mav.rc_channels_override_send(link.target_system, link.target_component, *rc)
                        if tuple(rc) != last_rc_values:
                            event("pilot_command", rc_pwm=rc.copy())
                            last_rc_values = tuple(rc)
                        last_rc = now
                    if release_at is not None:
                        x = physics.state.copy()
                        if abs(np.rad2deg(x[6])) > (60 if interactive or back_transition else 20) or abs(np.rad2deg(x[7])) > (60 if back_transition else 45 if interactive else 15):
                            raise RuntimeError("Diagnostic attitude envelope exceeded")
                        if back_transition and (not 5 < -x[2] < 120 or np.linalg.norm(x[3:6]) > 30):
                            raise RuntimeError("Forward-to-hover diagnostic altitude/speed envelope exceeded")
                        if now-release_at >= args.duration:
                            break
                else:
                    raise TimeoutError("Flight wall-clock deadline exceeded")
            stop.set()
            worker.join(timeout=2)
            if worker.is_alive() or failures:
                raise RuntimeError("Physics log did not close cleanly")
            metrics = ({"passed": True, "test_kind": "interactive_session", "validation_run": False}
                       if interactive else {} if is_doublet or back_transition else analyse_elevator_truth(output / "truth.jsonl", initial, trim, args.duration,
                                             surface_limit_deg=mapping.elevator.travel_deg)
                       if args.scenario == "elevator" else
                       analyse_truth(output / "truth.jsonl", initial, args.scenario, args.duration,
                                     perturbation=[args.perturb_roll, args.perturb_pitch]))
            report.update(metrics,
                          controlled_flight=not interactive and mode == 5, normal_arm=arm_ack, arm_parameters=observed,
                          release_at_s=release_at, direction_checks=signs,
                          last_airspeed=telemetry.get("AIRSPEED"), last_vfr_hud=telemetry.get("VFR_HUD"),
                          stationary_pitot_calibration=airspeed_calibrated,
                          air_start_at_s=air_started_at, ground_home_confirmed=home_confirmed,
                          perturbation_deg=[args.perturb_roll, args.perturb_pitch],
                          last_ekf_status=telemetry.get("EKF_STATUS_REPORT"), errors=failures)
        except Exception as exc:
            report.update(error=repr(exc), simulated_s=physics.time, errors=failures,
                          controlled_flight=not interactive and mode == 5 and release_at is not None,
                          normal_arm=arm_ack, arm_parameters=observed,
                          release_at_s=release_at, direction_checks=signs,
                          stationary_pitot_calibration=airspeed_calibrated,
                          ground_home_confirmed=home_confirmed)
            event("failed", error=repr(exc))
        finally:
            stop.set()
            if process is not None and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            worker.join(timeout=2)
            if link is not None:
                link.close()
            report["wall_s"] = time.monotonic()-started
            report["simulated_s"] = physics.time
            if is_doublet or back_transition:
                analysis = (analyse_doublet(output / "truth.jsonl", initial, plan, args.duration, release_at) if is_doublet else
                            analyse_back_transition(output / "truth.jsonl", release_at, args.duration, -initial[2],
                                                    telemetry.get("HEARTBEAT", {}).get("custom_mode")))
                report.update(analysis)
                if report.get("error") or failures:
                    report["passed"] = False
            (output / "events.json").write_text(json.dumps(events, indent=2) + "\n")
            (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
            print(json.dumps(report, indent=2), flush=True)
    return report["passed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=["manual", "elevator", "fbwa", "steps", "back_transition", *DOUBLETS])
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config")
    parser.add_argument("--forward-parameters", choices=("first_transition_2026_10_09",),
                        help="Saved forward baseline for a PT1.5 trim start; no gain tuning")
    parser.add_argument("--port", type=int, default=19002)
    parser.add_argument("--mavlink-port", type=int, default=15760)
    parser.add_argument("--instance", type=int, default=0)
    parser.add_argument("--duration", type=float)
    parser.add_argument("--interactive", action="store_true", help="Allow Flight Deck mode commands after bootstrap")
    parser.add_argument("--perturb-roll", type=float, default=0.)
    parser.add_argument("--perturb-pitch", type=float, default=0.)
    parser.add_argument("--wall-timeout", type=float, help="Host runtime allowance; leaves physics and controller timing unchanged")
    args = parser.parse_args()
    if args.wall_timeout is not None and (not np.isfinite(args.wall_timeout) or args.wall_timeout <= 0):
        parser.error("wall timeout must be finite and positive")
    if args.interactive and args.scenario != "fbwa":
        parser.error("interactive sessions bootstrap in FBWA")
    if args.duration is None:
        args.duration = {"steps": 118., "elevator": 45.}.get(args.scenario, 60.)
    if not np.isfinite(args.duration) or args.duration <= 0:
        parser.error("duration must be finite and positive")
    if args.scenario == "steps" and args.duration < 118:
        parser.error("six step tests require at least 118 seconds")
    if args.scenario == "elevator" and args.duration < 45:
        parser.error("elevator demonstration requires at least 45 seconds")
    if args.scenario in DOUBLETS and args.duration < 45:
        parser.error("doublets require at least 45 seconds including recovery")
    if args.scenario == "back_transition" and args.duration < 60:
        parser.error("forward-to-hover requires at least 60 seconds")
    if not 0 <= args.instance <= 20:
        parser.error("instance must be between 0 and 20")
    if not np.isfinite([args.perturb_roll, args.perturb_pitch]).all() or max(abs(args.perturb_roll), abs(args.perturb_pitch)) > .5:
        parser.error("release perturbations must be within +/-0.5 degrees")
    raise SystemExit(0 if run(args) else 1)


if __name__ == "__main__":
    main()
