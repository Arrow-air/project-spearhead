"""Fit the KDE5215XF-220 electrical model used for the KDE pusher option.

The KDE static test (KDE5215XF-220 with the KDE-CF215-DP 21.5x7.3 propeller,
`data/propulsion/kde5215xf_220_cf215_dp.csv`) gives voltage, throttle, battery
current, RPM and torque. With Kv = 220 RPM/V and Kt = 60/(2*pi*Kv):

    I_winding = Q/Kt + I0
    V_motor   = RPM/Kv + R_equivalent * I_winding
    I_battery = V_motor * I_winding / (eta_ESC * V_battery)

R_equivalent and I0 are fitted to the battery current of all 28 points
(eta_ESC = 0.97, as for the MAD estimate). The effective voltage fraction
V_motor/V_battery at each test point gives a throttle-to-duty map over
throttle and battery voltage; the KDE ESC does not apply full voltage at
full throttle, so commanded throttle is not used as duty directly.

The propeller is not part of this fit: the pusher model combines these motor
constants with the 20x12 propeller map. Writes
`data/propulsion/kde5215xf_220_pusher_fit.json` (or --output).
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

FDM_ROOT = Path(__file__).resolve().parents[1]
SOURCE = FDM_ROOT / "data/propulsion/kde5215xf_220_cf215_dp.csv"
KV = 220.
ESC_EFFICIENCY = .97


def load_points(path=SOURCE):
    rows = list(csv.DictReader(Path(path).open()))
    return {key: np.array([float(row[key]) for row in rows])
            for key in ("voltage_v", "throttle", "rpm", "torque_nm", "current_a")}


def fit(points):
    kt = 60/(2*np.pi*KV)
    v, rpm, torque, current = (points[k] for k in ("voltage_v", "rpm", "torque_nm", "current_a"))

    def battery_current(p):
        resistance, idle = p
        winding = torque/kt+idle
        return (rpm/KV+resistance*winding)*winding/(ESC_EFFICIENCY*v)

    solution = least_squares(lambda p: (battery_current(p)-current)/current, [.15, 1.], bounds=([0, 0], [2, 5]))
    resistance, idle = (float(x) for x in solution.x)
    duty = (rpm/KV+resistance*(torque/kt+idle))/v
    voltages = sorted(set(points["voltage_v"].tolist()))
    throttles = sorted(set(points["throttle"].tolist()))
    table = [[float(duty[(v == volt) & (points["throttle"] == thr)][0]) for thr in throttles] for volt in voltages]
    relative = (battery_current(solution.x)-current)/current
    return dict(resistance=resistance, idle=idle, voltages=voltages, throttles=throttles, duty=table,
                current_rms=float(np.sqrt(np.mean(relative**2))), current_max=float(np.max(abs(relative))))


def duty_from_map(throttle, voltage, voltages, throttles, duty):
    """Bilinear in throttle and voltage; linear ramp from zero below the first throttle."""
    voltage = float(np.clip(voltage, voltages[0], voltages[-1]))
    column = [np.interp(throttle, [0., *throttles], [0., *row]) for row in duty]
    return float(np.interp(voltage, voltages, column))


def check_against_test(points, result):
    """Re-predict the static test with a quadratic propeller law fitted to the same data."""
    kt = 60/(2*np.pi*KV)
    rpm, torque, current = points["rpm"], points["torque_nm"], points["current_a"]
    c = float(np.sum(torque*rpm**2)/np.sum(rpm**4))
    predicted_rpm, predicted_current = [], []
    for voltage, throttle in zip(points["voltage_v"], points["throttle"]):
        d = duty_from_map(throttle, voltage, result["voltages"], result["throttles"], result["duty"])
        n = brentq(lambda n: n/KV+result["resistance"]*(c*n*n/kt+result["idle"])-d*voltage, 1, 20000)
        winding = c*n*n/kt+result["idle"]
        predicted_rpm.append(n)
        predicted_current.append((n/KV+result["resistance"]*winding)*winding/(ESC_EFFICIENCY*voltage))
    rpm_error = (np.array(predicted_rpm)-rpm)/rpm
    current_error = (np.array(predicted_current)-current)/current
    return dict(rpm_rms=float(np.sqrt(np.mean(rpm_error**2))), rpm_max=float(np.max(abs(rpm_error))),
                current_rms=float(np.sqrt(np.mean(current_error**2))), current_max=float(np.max(abs(current_error))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=FDM_ROOT / "data/propulsion/kde5215xf_220_pusher_fit.json")
    args = parser.parse_args()
    points = load_points()
    result = fit(points)
    check = check_against_test(points, result)
    record = {
        "status": "provisional estimate; no installed KDE pusher / 20x12 bench test",
        "source_csv": "kde5215xf_220_cf215_dp.csv",
        "source_csv_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_workbook": "sources/KDE5215-220kV_rev2.xlsx (see kde5215xf_220_source.json)",
        "test_esc": "KDE-UAS125UVC (the aircraft uses KDE-UAS75HVC+)",
        "equation": "I_w = Q/Kt + I0; V_m = RPM/Kv + R_eq*I_w; I_batt = V_m*I_w/(eta*V); Kt = 60/(2*pi*Kv)",
        "kv_rpm_per_v": KV,
        "esc_efficiency_assumed": ESC_EFFICIENCY,
        "equivalent_resistance_ohm": result["resistance"],
        "no_load_current_a": result["idle"],
        "battery_current_fit_rms_relative": result["current_rms"],
        "battery_current_fit_max_relative": result["current_max"],
        "throttle_duty_map": {"voltages_v": result["voltages"], "throttles": result["throttles"],
                              "duty": result["duty"]},
        "static_test_check": check,
        "published_phase_resistance_ohm": 0.078,
    }
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps({k: record[k] for k in ("equivalent_resistance_ohm", "no_load_current_a",
                                             "battery_current_fit_rms_relative", "static_test_check")}, indent=2))


if __name__ == "__main__":
    main()
