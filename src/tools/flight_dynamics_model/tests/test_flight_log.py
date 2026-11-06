import json
import math
from pathlib import Path

import numpy as np
import pytest

from spearhead.flight_log import (hover_summary, log_takeoffs, log_traces, parameter_diff, read_parameter_file,
                                  run_parameters, run_traces)

SAMPLE_LOG = Path("/mnt/user-data/uploads/00000016.BIN")


def test_parameter_files_in_mission_planner_and_qgc_formats(tmp_path):
    mp = tmp_path / "a.param"; mp.write_text("# comment\nQ_M_SPIN_MIN,0.15\nSCHED_LOOP_RATE 300\n")
    qgc = tmp_path / "b.params"; qgc.write_text("1\t1\tQ_M_SPIN_MIN\t0.17\t9\n")
    assert read_parameter_file(mp) == {"Q_M_SPIN_MIN": .15, "SCHED_LOOP_RATE": 300.}
    assert read_parameter_file(qgc) == {"Q_M_SPIN_MIN": .17}


def test_parameter_diff_keeps_flight_settings_and_drops_hardware_noise():
    aircraft = {"SCHED_LOOP_RATE": 300, "COMPASS_OFS_X": 12, "MAV1_EXTRA1": 4, "RC7_MIN": 1000, "RC3_TRIM": 1000,
                "Q_M_SPIN_MIN": .15, "ONLY_ON_AIRCRAFT": 1}
    sitl = {"SCHED_LOOP_RATE": 100, "COMPASS_OFS_X": 0, "MAV1_EXTRA1": 10, "RC7_MIN": 1100, "RC3_TRIM": 1500,
            "Q_M_SPIN_MIN": .15000001, "Q_P_D_JERK": 5}
    rows = {row["name"]: row for row in parameter_diff(aircraft, sitl)}
    assert set(rows) == {"SCHED_LOOP_RATE", "RC3_TRIM", "Q_P_D_JERK"}
    assert rows["Q_P_D_JERK"]["status"] == "SITL only" and rows["SCHED_LOOP_RATE"]["aircraft"] == "300"


def _rows(times, **columns):
    return [{"TimeUS": int(t*1e6), **{k: v[i] for k, v in columns.items()}} for i, t in enumerate(times)]


def test_takeoffs_split_by_landing_and_traces_align_to_takeoff():
    t = np.arange(0, 100, .1)
    flying = ((t > 10) & (t < 30)) | ((t > 60) & (t < 90))
    pwm = np.where(flying, 1600, 1100)
    series = {"RCOU": _rows(t, C5=pwm, C6=pwm, C7=pwm-50, C8=pwm-50),
              "QTUN": _rows(t, Alt=np.where(flying, 10., 0.)),
              "ATT": _rows(t, Roll=np.full(t.size, 2.), Pitch=np.zeros(t.size))}
    flights = log_takeoffs(series)
    assert [round(s) for s, _ in flights] == [10, 60] and round(flights[1][1]) == 30
    traces = log_traces(series, flights[1][0], 20)
    assert traces["lift_FR"][50] == 1600 and traces["lift_FL"][50] == 1550 and traces["roll_deg"][0] == 2
    summary = hover_summary(traces)
    assert summary["front_share_pwm"] == pytest.approx((500+450)/(2*500+2*450))


def test_run_traces_and_parameters_from_a_recorded_run(tmp_path):
    rows = []
    for i in range(300):
        held = i < 100
        rows.append({"time_s": i*.05, "hold": held, "pwm": [1500]*4 + [1650, 1640, 1630, 1620] + [0]*8,
                     "state": [0, 0, -20, 0, 0, 0, math.radians(1), 0, 0, 0, 0, 0],
                     "sensor": {"battery": {"voltage": 46., "current": 80.}}})
    (tmp_path / "truth.jsonl").write_text("".join(json.dumps(r)+"\n" for r in rows))
    (tmp_path / "quad.parm").write_text("Q_M_SPIN_MIN 0.15\nSCHED_LOOP_RATE 50\nSCHED_LOOP_RATE 100\n")
    (tmp_path / "manifest.json").write_text(json.dumps({"parameter_overlay": {"Q_P_D_JERK": 5}}))
    traces = run_traces(tmp_path)
    assert traces["time_s"][0] == 0 and traces["lift_FR"][0] == 1650 and traces["lift_RR"][0] == 1620
    assert traces["roll_deg"][3] == pytest.approx(1) and traces["battery_v"][0] == 46
    parameters, source = run_parameters(tmp_path)
    assert parameters == {"Q_M_SPIN_MIN": .15, "SCHED_LOOP_RATE": 100., "Q_P_D_JERK": 5.}
    assert source == "parameter files and overlay"


@pytest.mark.skipif(not SAMPLE_LOG.is_file(), reason="PT1 DataFlash sample not available")
def test_reads_a_real_pt1_log():
    from spearhead.flight_log import read_dataflash
    parameters, series = read_dataflash(SAMPLE_LOG)
    assert parameters["Q_M_PWM_MIN"] == 1100 and len(log_takeoffs(series)) >= 1
