import json

import numpy as np
import pytest

from spearhead.gui.flight_setup import FlightSetup, maneuver_options
from spearhead.sitl.hover_doublets import HOVER_DOUBLETS, analyse_hover_doublet, hover_doublet_command, hover_doublet_plan

PARAMS = {"Q_M_THST_HOVER": .474, "Q_A_ANGLE_MAX": 20}


@pytest.mark.parametrize("name", HOVER_DOUBLETS)
def test_hover_selection_keeps_lift_mapping_and_uses_quad_runner(name):
    setup = FlightSetup(maneuver=name, duration_s=45).validate()
    assert setup.scenario == "quad_"+name
    assert setup.config()["flight_regime"] == "hover"
    assert len(setup.config()["channels"]["lift"]) == 4
    assert "elevator" not in maneuver_options("pt1_as_built", "hover20")["pitch_doublet"]
    with pytest.raises(ValueError, match="45"):
        FlightSetup(maneuver=name, duration_s=44).validate()


@pytest.mark.parametrize("name", HOVER_DOUBLETS)
def test_input_inverse_and_unconfirmed_mode_guard(name):
    plan = hover_doublet_plan(name, PARAMS)
    mode = plan["excitation_mode"]
    first = hover_doublet_command(plan, 10.1, mode)[1][plan["channel"]-1]
    second = hover_doublet_command(plan, 12.1, mode)[1][plan["channel"]-1]
    assert first > plan["neutral_pwm"] > second
    assert (first-second)/2/plan["pwm_per_unit"] == pytest.approx(plan["amplitude"], abs=.001)
    assert hover_doublet_command(plan, 15, mode)[0] == 18
    # A delayed QHOVER confirmation keeps calculated collective, avoiding a drop.
    assert hover_doublet_command(plan, 3.1, 17)[1][2] == plan["hover_pwm"]
    if mode == 17:
        assert hover_doublet_command(plan, 10.1, 18)[1][2] == 1500
    else:
        assert hover_doublet_command(plan, 10.1, 17)[1][plan["channel"]-1] == 1500


@pytest.mark.parametrize("name", HOVER_DOUBLETS)
def test_report_requires_both_readback_lobes_fresh_modes_and_neutral(tmp_path, name):
    plan = hover_doublet_plan(name, PARAMS)
    initial = [0, 0, -20]+[0]*9
    path = tmp_path/'truth.jsonl'
    path.write_text(''.join(json.dumps({"time_s": 30+t, "hold": False, "state": initial,
                           "pwm": [1500]*8, "control": {"throttle": 0}})+'\n' for t in np.arange(0, 45.01, .05)))
    records = []
    for t in np.arange(0, 45.01, .1):
        mode = plan["excitation_mode"] if 9 <= t < 15 else 18
        rc = hover_doublet_command(plan, t, mode)[1]
        records.extend([{"fdm_time_s": 30+t, "message": {"mavpackettype": "HEARTBEAT", "custom_mode": mode}},
                        {"fdm_time_s": 30+t, "message": {"mavpackettype": "RC_CHANNELS", f"chan{plan['channel']}_raw": rc[plan['channel']-1]}}])
    log = tmp_path/'telemetry.jsonl'
    def save(): log.write_text(''.join(json.dumps(r)+'\n' for r in records))
    save()
    assert analyse_hover_doublet(path, 30, 45, initial, plan)["passed"]
    for r in records:
        if 42.3 <= r['fdm_time_s'] < 43.9 and r['message']['mavpackettype'] == 'RC_CHANNELS':
            r['message'][f"chan{plan['channel']}_raw"] = plan['neutral_pwm']
    save()
    assert not analyse_hover_doublet(path, 30, 45, initial, plan)["passed"]
    log.write_text('')
    assert not analyse_hover_doublet(path, 30, 45, initial, plan)["passed"]


def test_hover_plan_and_plot_use_hover_units_and_actual_attitude():
    pytest.importorskip('plotly')
    from spearhead.gui.sitl_page import flight_plan_text, sitl_figures
    plan = hover_doublet_plan('pitch_doublet', PARAMS)
    row = {"time_s": 40.5, "hold": False, "state": [0, 0, -20, 0, 0, 0, 0, .1, 0, 0, 0, 0],
           "pwm": [1500]*8, "control": {"de": 0, "da": 0, "dr": 0, "throttle": 0}}
    figures = sitl_figures({"history": [row], "manifest": {"initial_state": [0]*12, "trim_control": row['control'], "doublet_plan": plan},
                           "report": {"release_at_s": 30}})
    assert figures['Doublet input'].data[0].y == (8,)
    assert figures['Doublet input'].data[2].y == pytest.approx([np.rad2deg(.1)])
    text = flight_plan_text(FlightSetup(maneuver='pitch_doublet', duration_s=45))
    assert 'QHOVER' in text and '8°' in text and 'FBWA' not in text


@pytest.mark.parametrize("name", HOVER_DOUBLETS)
def test_manager_launches_hover_doublet_with_correct_runner(tmp_path, monkeypatch, name):
    from types import SimpleNamespace
    from spearhead.gui import sitl_session
    calls = []
    process = SimpleNamespace(pid=987654, returncode=None, poll=lambda: None)
    monkeypatch.setattr(sitl_session, "verify_firmware", lambda path: {})
    monkeypatch.setattr(sitl_session.subprocess, "Popen", lambda command, **kwargs: calls.append(command) or process)
    monkeypatch.setattr(sitl_session.threading.Thread, "start", lambda self: None)
    manager = sitl_session.SITLManager(root=tmp_path/"runs")
    manager.start_configured(FlightSetup(maneuver=name, duration_s=45))
    command = calls[0]
    assert command[2] == "spearhead.sitl.quad"
    assert command[command.index("--doublet")+1] == name
    assert command[command.index("--duration")+1] == "45"


def test_observer_preserves_actual_hover_rc_readback_for_chart(tmp_path):
    from spearhead.gui.sitl_session import RunObserver
    plan = hover_doublet_plan("pitch_doublet", PARAMS)
    (tmp_path/"manifest.json").write_text(json.dumps({"doublet_plan": plan}))
    (tmp_path/"telemetry.jsonl").write_text(json.dumps({"fdm_time_s": 40.5, "message": {
        "mavpackettype": "RC_CHANNELS", "chan1_raw": 1500, "chan2_raw": 1700, "chan3_raw": 1500, "chan4_raw": 1500}})+'\n')
    observer = RunObserver(tmp_path); observer.poll()
    assert observer.snapshot()["pilot_inputs"] == [{"time_s": 40.5, "channels": [1500, 1700, 1500, 1500]}]
