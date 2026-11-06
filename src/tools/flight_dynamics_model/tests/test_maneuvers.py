import json
from dataclasses import asdict

import numpy as np
import pytest

from spearhead.controls import Control
from spearhead.gui.flight_setup import FlightSetup, maneuver_options
from spearhead.sitl.actuators import ActuatorMap, Channel
from spearhead.sitl.back_transition import analyse_back_transition, back_transition_command, forward_release_ready
from spearhead.sitl.doublets import DOUBLETS, analyse_doublet, doublet_plan, doublet_delta, doublet_mode
from spearhead.sitl.flight import bootstrap_mode, manual_surface_rc, surface_pwm


MAPPING = ActuatorMap(aileron=Channel(1, travel_deg=30), elevator=Channel(2, reverse=True, travel_deg=30),
                      rudder=Channel(4, travel_deg=30), pusher=Channel(3, minimum=1050, maximum=1940))
TRIM = Control(de=np.deg2rad(-8), da=.001, dr=.005, throttle=.88)


@pytest.mark.parametrize("name", list(DOUBLETS))
def test_doublet_timing_and_actual_servo_inverse(name):
    plan = doublet_plan(name, TRIM, MAPPING)
    assert doublet_delta(plan, 9.99) == 0
    assert doublet_delta(plan, 10) == plan["amplitude"]
    assert doublet_delta(plan, 10+plan["pulse_s"]) == -plan["amplitude"]
    assert doublet_delta(plan, 10+2*plan["pulse_s"]) == 0
    assert doublet_mode(plan, 8.99) == 5
    assert doublet_mode(plan, 9) == 0
    assert doublet_mode(plan, plan["recovery_at_s"]) == 5
    if name != "throttle_doublet":
        c = getattr(MAPPING, plan["channel"]); neutral = surface_pwm(c, getattr(TRIM, plan["surface"]))
        for delta in [-plan["amplitude"], 0, plan["amplitude"]]:
            rc = manual_surface_rc(c, getattr(TRIM, plan["surface"]), neutral, delta)
            span = 2000-neutral if rc >= 1500 else neutral-1000
            pwm = [1500]*4; pwm[c.number-1] = neutral+int((rc-1500)/500*span)
            assert np.rad2deg(c.surface(pwm)-getattr(TRIM, plan["surface"])) == pytest.approx(delta, abs=.07)


def test_throttle_headroom_is_explicit_and_inadequate_excitation_rejected():
    plan = doublet_plan("throttle_doublet", Control(0, 0, 0, .95), MAPPING)
    assert plan["amplitude"] == pytest.approx(.045)
    with pytest.raises(ValueError, match="authority"):
        doublet_plan("throttle_doublet", Control(0, 0, 0, .98), MAPPING)


@pytest.mark.parametrize("name", ["back_transition", *DOUBLETS])
def test_forward_maneuver_routes_correct_start_and_motor_mapping(name):
    setup = FlightSetup(initial_point="trim24", speed_m_s=24, altitude_m=100, maneuver=name, cg_case="forward25")
    assert setup.scenario == name
    config = setup.config()
    assert config["flight_regime"] == "fixed_wing" and config["speed_m_s"] == 24
    assert ("lift" in config["channels"]) == (name == "back_transition")
    assert name in maneuver_options("pt15", "trim24")
    assert (name in maneuver_options("pt15", "hover20")) == (name in DOUBLETS)


def test_back_transition_arms_normally_in_q_mode_then_releases_fw_and_retains_throttle_until_confirmation():
    assert bootstrap_mode("back_transition", None, 25) == 17
    assert bootstrap_mode("back_transition", 25, 27.9) == 17
    assert bootstrap_mode("back_transition", 25, 28) == 5
    assert back_transition_command(9, 5, .88) == (5, 1887)
    assert back_transition_command(10, 5, .88) == (19, 1887)
    assert back_transition_command(10, 19, .88) == (19, 1500)
    params = {"Q_M_PWM_MIN": 1100, "Q_M_PWM_MAX": 1940, "Q_M_SPIN_ARM": .1}
    pwm = [1500, 1637, 1833, 1505]+[1184]*4
    assert forward_release_ready(pwm, MAPPING, params, .88, True)
    assert not forward_release_ready(pwm, MAPPING, params, .88, False)
    assert not forward_release_ready(pwm[:4]+[1700]*4, MAPPING, params, .88, True)
    assert not forward_release_ready([1500, 1637, 1050, 1505]+[1184]*4, MAPPING, params, .88, True)


def telemetry(path, mode_at, duration=60):
    path.write_text("".join(json.dumps({"fdm_time_s": 30+t, "message": {"mavpackettype": "HEARTBEAT", "custom_mode": mode_at(t)}})+"\n"
                            for t in np.arange(0, duration+1, .5)))


@pytest.mark.parametrize("name", list(DOUBLETS))
def test_doublet_analysis_checks_both_actual_lobes_and_mode_evidence(tmp_path, name):
    plan = doublet_plan(name, TRIM, MAPPING); initial = np.array([0, 0, -100, 24, 0, 0, 0, .01, 0, 0, 0, 0])
    rows = []
    for t in np.arange(0, 45.02, .02):
        control = asdict(TRIM); delta = doublet_delta(plan, t)
        control[plan["surface"]] += np.deg2rad(delta) if plan["unit"] == "deg" else delta
        rows.append({"time_s": 30+t, "hold": False, "state": initial.tolist(), "control": control})
    path = tmp_path/"truth.jsonl"
    def save(): path.write_text("".join(json.dumps(r)+"\n" for r in rows))
    save(); telemetry(tmp_path/"telemetry.jsonl", lambda t: doublet_mode(plan, t), 45)
    report = analyse_doublet(path, initial, plan, 45, 30)
    assert report["passed"] and report["returned_near_trim"]
    for r in rows:
        if doublet_delta(plan, r["time_s"]-30) < 0: r["control"][plan["surface"]] = getattr(TRIM, plan["surface"])
    save()
    assert not analyse_doublet(path, initial, plan, 45, 30)["passed"]
    telemetry(tmp_path/"telemetry.jsonl", lambda t: 5, 45)
    assert not analyse_doublet(path, initial, plan, 45, 30)["passed"]


def test_hover_acceptance_requires_pusher_off_lift_active_and_fresh_mode_throughout(tmp_path):
    rows = [{"time_s": 30+t, "hold": False, "state": [0, 0, -100, 24 if t < 20 else .1, 0, 0, 0, 0, 0, 0, 0, 0],
             "pwm": [1500, 1500, 1050, 1500]+[1670]*4, "control": {"throttle": .88 if t < 10 else 0}}
            for t in np.arange(0, 60.02, .02)]
    path = tmp_path/"truth.jsonl"
    def save(): path.write_text("".join(json.dumps(r)+"\n" for r in rows))
    save(); telemetry(tmp_path/"telemetry.jsonl", lambda t: 5 if t < 10 else 19)
    report = analyse_back_transition(path, 30, 60, 100, 19)
    assert report["passed"] and report["longest_hover_s"] > 39
    for r in rows: r["control"]["throttle"] = .1
    save()
    assert not analyse_back_transition(path, 30, 60, 100, 19)["passed"]
    for r in rows: r["control"]["throttle"] = 0
    save(); telemetry(tmp_path/"telemetry.jsonl", lambda t: 5 if t < 59 else 19)
    assert not analyse_back_transition(path, 30, 60, 100, 19)["passed"]


def test_doublet_plots_show_physical_requested_input_and_yaw_response():
    pytest.importorskip("plotly")
    from spearhead.gui.sitl_page import sitl_figures
    plan = doublet_plan("pitch_doublet", TRIM, MAPPING)
    state = [0, 0, -100, 24, .2, 0, 0, .01, .1, .02, .03, .04]
    control = asdict(TRIM); control["de"] += np.deg2rad(3)
    row = {"time_s": 40.3, "hold": False, "state": state, "control": control, "pwm": [1500]*8}
    figures = sitl_figures({"history": [row], "manifest": {"initial_state": [0]*12, "trim_control": asdict(TRIM), "doublet_plan": plan},
                            "report": {"release_at_s": 30}})
    assert figures["Doublet input"].data[0].y == (3,)
    assert figures["Doublet input"].data[1].y[0] == pytest.approx(3)
    assert figures["Body rates"].data[2].y[0] == pytest.approx(np.rad2deg(.04))
    assert figures["Yaw / sideslip"].data[0].y[0] == pytest.approx(np.rad2deg(.1))


def test_flight_plots_count_time_from_release():
    pytest.importorskip("plotly")
    from spearhead.gui.sitl_page import flight_time_origin, sitl_figures
    state = [0, 0, -20, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    control = asdict(TRIM)
    rows = [{"time_s": t, "hold": t < 47, "state": state, "control": control, "pwm": [1500]*8} for t in (10., 47., 107.)]
    snapshot = {"history": rows, "report": {"release_at_s": 47.}}
    assert flight_time_origin(snapshot) == 47.
    figures = sitl_figures(snapshot)
    assert tuple(figures["Roll"].data[0].x) == (-37., 0., 60.)
    assert flight_time_origin({"history": rows, "time_s": 107., "free_time_s": 60.}) == 47.
    assert flight_time_origin({"history": rows}) == 47.
