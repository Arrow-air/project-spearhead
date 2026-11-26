"""Fit the HOBBYWING PLATINUM 5230SL 380KV electrical model for the PT1.5 pusher.

The manufacturer sheet (`data/propulsion/hobbywing_5230sl_380_source.json`)
gives three 29.6 V throttle sweeps (APC 17x10E, 18x10E, 18x12E) with battery
current, RPM and torque. With Kv = 380 RPM/V and Kt = 60/(2*pi*Kv):

    I_winding = Q/Kt + I0
    V_motor   = RPM/Kv + R_equivalent * I_winding  (= throttle * V_battery)
    I_battery = V_motor * I_winding / (eta_ESC * V_battery)

R_equivalent and I0 are fitted to the commanded duty and the battery current
of the points at or below the AMPX 80A continuous rating, the range the pusher
can use (53 of 66 points; eta_ESC = 0.97, as for the other pusher estimates).
The sheet states nominal pack voltage, so the high-current rows also carry
unrecorded battery sag. The
sheet's duty follows throttle within a few percent, so the pusher model uses
duty = throttle (as for the MAD AMPX sheet). The propeller is not part of this
fit: the pusher model combines these constants with the APC 20x12 map.

Writes `data/propulsion/hobbywing_5230sl_380_pusher_fit.json` (or --output).
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

FDM_ROOT = Path(__file__).resolve().parents[1]
SOURCE = FDM_ROOT/"data/propulsion/hobbywing_5230sl_380_source.json"
ESC_EFFICIENCY = .97
MAX_FIT_CURRENT_A = 80.  # MAD AMPX 80A continuous rating


def load_points(path=SOURCE):
    data = json.loads(Path(path).read_text())
    rows = np.vstack([np.asarray(table["samples"], float) for table in data["tables"].values()])
    points = {name: rows[:, i] for i, name in enumerate(data["columns"])}
    points["voltage_v"] = np.full(len(rows), float(data["test_voltage_v"]))
    return data, points


def fit(data, points):
    kv = float(data["specifications"]["kv_rpm_per_v"])
    kt = 60/(2*np.pi*kv)
    v, throttle, rpm, torque, current = (points[k] for k in ("voltage_v", "throttle", "rpm", "torque_nm", "current_a"))

    def model(p):
        resistance, idle = p
        winding = torque/kt+idle
        motor_voltage = rpm/kv+resistance*winding
        return motor_voltage/v, motor_voltage*winding/(ESC_EFFICIENCY*v)

    def residual(p):
        duty, battery = model(p)
        return np.r_[duty-throttle, (battery-current)/current]

    solution = least_squares(residual, [.04, 4.], bounds=([0, 0], [1, 20]))
    duty, battery = model(solution.x)
    relative = (battery-current)/current
    return dict(resistance=float(solution.x[0]), idle=float(solution.x[1]),
                duty_error_max=float(np.max(abs(duty-throttle))),
                current_rms=float(np.sqrt(np.mean(relative**2))), current_max=float(np.max(abs(relative))))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=FDM_ROOT/"data/propulsion/hobbywing_5230sl_380_pusher_fit.json")
    args = parser.parse_args()
    data, points = load_points()
    used = points["current_a"] <= MAX_FIT_CURRENT_A
    points = {key: value[used] for key, value in points.items()}
    result = fit(data, points)
    record = {
        "source": SOURCE.name, "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "equations": __doc__.split("\n\n")[1].strip(),
        "kv_rpm_per_v": data["specifications"]["kv_rpm_per_v"], "esc_efficiency_assumed": ESC_EFFICIENCY,
        "equivalent_resistance_ohm": result["resistance"], "no_load_current_a": result["idle"],
        "published_phase_resistance_ohm": data["specifications"]["phase_resistance_ohm"],
        "published_no_load_current_a": data["specifications"]["no_load_current_a"],
        "points": int(len(points["throttle"])), "max_fit_current_a": MAX_FIT_CURRENT_A,
        "duty_minus_throttle_max": result["duty_error_max"],
        "battery_current_relative_rms": result["current_rms"],
        "battery_current_relative_max": result["current_max"],
        "limitations": [
            "Fitted at 29.6 V only; PT1.5 runs it at 12S (36-50.4 V), outside the 6-8S rating, behind an output limit.",
            "Test ESC is not stated; the AMPX 80A is assumed to apply duty = throttle as in the MAD sheet.",
            "Constant I0: iron losses rising with RPM are absorbed in the fit only within the tested RPM range.",
            "No thermal model; the sheet's full-throttle runs reached 80 C within 1.4-6 minutes.",
        ],
    }
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps({key: record[key] for key in ("equivalent_resistance_ohm", "no_load_current_a",
                                                    "duty_minus_throttle_max", "battery_current_relative_rms",
                                                    "battery_current_relative_max")}, indent=2))


if __name__ == "__main__":
    main()
