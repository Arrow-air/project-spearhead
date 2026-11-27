"""Build the PT1.5 lift curve: MAD V8013 150KV + AMPX 80A + FLUXER PRO 26x7.8.

The manufacturer sheet (`data/propulsion/mad_v8013_150_fluxer26_source.json`)
gives one 48 V throttle sweep at sea level and 25 C. The FDM lift lookup takes
battery voltage and throttle and applies no density correction, so this script
rebuilds the curve at the profile density (1.021 kg/m^3, Ankara hot day) over a
battery-voltage grid:

    Propeller (static):  T = Ct(RPM) rho n^2 D^4,  Q = Cq(RPM) rho n^2 D^5
                         Ct, Cq linear in RPM, fitted to the sheet and held
                         constant outside the measured RPM range.
    Motor and ESC:       throttle * V = RPM/Kv + R_eq (I0 + Q/Kt)
                         I_battery = (throttle * V) (I0 + Q/Kt) / (eta_ESC V)
                         R_eq and I0 fitted to the sheet thrust and current
                         through this closed-loop model,
                         eta_ESC = 0.97 as for the pusher estimates.

Rows below the sheet's first throttle (30 %) come from the same model and are
marked `model_below_sheet`. All rows are static: no inflow, ground effect,
thermal or battery-sag model.

Writes `data/propulsion/mad_v8013_150_fluxer26_rho1021.csv` and
`data/propulsion/mad_v8013_150_fluxer26_fit.json` (or --output-dir).
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

FDM_ROOT = Path(__file__).resolve().parents[1]
PROPULSION = FDM_ROOT / "data/propulsion"
SOURCE = PROPULSION / "mad_v8013_150_fluxer26_source.json"
GRAVITY = 9.80665
ESC_EFFICIENCY = .97
TARGET_DENSITY = 1.021
VOLTAGES = (36.0, 40.0, 42.0, 44.4, 48.0, 50.4)
THROTTLES = tuple(round(.15 + .05*i, 2) for i in range(18))  # 0.15 ... 1.00


def load_sheet(path=SOURCE):
    data = json.loads(Path(path).read_text())
    table = np.asarray(data["samples"], dtype=float)
    columns = {name: table[:, i] for i, name in enumerate(data["columns"])}
    return data, columns


def fit(data, sheet):
    spec, prop = data["specifications"], data["propeller_specifications"]
    kv = float(spec["kv_rpm_per_v"])
    kt = 60/(2*np.pi*kv)
    diameter = prop["diameter_in"]*.0254
    rho = float(data["test_density_kg_m3"])
    rpm, torque = sheet["rpm"], sheet["torque_nm"]
    n = rpm/60
    ct = sheet["thrust_gf"]*GRAVITY/1000/(rho*n**2*diameter**4)
    cq = torque/(rho*n**2*diameter**5)
    ct_line, cq_line = np.polyfit(rpm, ct, 1), np.polyfit(rpm, cq, 1)
    throttle, voltage, current = sheet["throttle"], sheet["voltage_v"], sheet["current_a"]

    def residual(p):
        resistance, idle = p
        winding = idle+torque/kt
        motor_voltage = rpm/kv+resistance*winding
        return np.r_[(motor_voltage-throttle*voltage)/voltage,
                     (motor_voltage*winding/(ESC_EFFICIENCY*voltage)-current)/current]

    # Start from the open-loop fit (measured torque and RPM), then refine R_eq and I0
    # on the closed-loop thrust and current the curve actually reproduces.
    start = least_squares(residual, [.1, 2.], bounds=([0, 0], [1, 10])).x
    model = dict(kv=kv, kt=kt, diameter_m=diameter, ct_line=[float(x) for x in ct_line],
                 cq_line=[float(x) for x in cq_line], rpm_range=[float(rpm.min()), float(rpm.max())],
                 test_density=rho)

    def closed_loop(p):
        trial = dict(model, resistance=p[0], idle=p[1])
        errors = []
        for i, command in enumerate(throttle):
            point = operating_point(trial, voltage[i], command, rho)
            errors += [point["thrust_n"]/(sheet["thrust_gf"][i]*GRAVITY/1000)-1, point["current_a"]/current[i]-1]
        return np.asarray(errors)

    resistance, idle = (float(x) for x in least_squares(closed_loop, start, bounds=([0, 0], [1, 10])).x)
    return dict(model, resistance=resistance, idle=idle)


def operating_point(model, voltage, throttle, rho):
    """Static steady state for a commanded throttle (duty) and battery voltage."""
    low, high = model["rpm_range"]

    def coefficients(rpm):
        clipped = float(np.clip(rpm, low, high))
        return np.polyval(model["ct_line"], clipped), np.polyval(model["cq_line"], clipped)

    def torque(rpm):
        return coefficients(rpm)[1]*rho*(rpm/60)**2*model["diameter_m"]**5

    def balance(rpm):
        return rpm/model["kv"]+model["resistance"]*(model["idle"]+torque(rpm)/model["kt"])-throttle*voltage

    if balance(0.) >= 0:
        raise ValueError("No rotating solution below the idle voltage drop")
    rpm = brentq(balance, 0., throttle*voltage*model["kv"], xtol=1e-6)
    ct, cq = coefficients(rpm)
    n = rpm/60
    shaft_torque = cq*rho*n**2*model["diameter_m"]**5
    winding = model["idle"]+shaft_torque/model["kt"]
    current = throttle*voltage*winding/(ESC_EFFICIENCY*voltage)
    return dict(thrust_n=float(ct*rho*n**2*model["diameter_m"]**4), current_a=float(current),
                power_w=float(current*voltage), rpm=float(rpm), torque_nm=float(shaft_torque))


def validation(model, sheet):
    rows = []
    for i, throttle in enumerate(sheet["throttle"]):
        point = operating_point(model, sheet["voltage_v"][i], throttle, model["test_density"])
        rows.append(dict(throttle=float(throttle),
                         thrust_error=point["thrust_n"]/(sheet["thrust_gf"][i]*GRAVITY/1000)-1,
                         current_error=point["current_a"]/sheet["current_a"][i]-1,
                         rpm_error=point["rpm"]/sheet["rpm"][i]-1,
                         torque_error=point["torque_nm"]/sheet["torque_nm"][i]-1))
    summary = {key: float(np.max(np.abs([row[key] for row in rows])))
               for key in ("thrust_error", "current_error", "rpm_error", "torque_error")}
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, default=PROPULSION)
    args = parser.parse_args()
    data, sheet = load_sheet()
    model = fit(data, sheet)
    rows, summary = validation(model, sheet)
    first_sheet_throttle = float(sheet["throttle"].min())
    curve_path = args.output_dir/"mad_v8013_150_fluxer26_rho1021.csv"
    with curve_path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["voltage_v", "throttle", "thrust_n", "current_a", "power_w", "rpm", "torque_nm", "basis"])
        for voltage in VOLTAGES:
            for throttle in THROTTLES:
                point = operating_point(model, voltage, throttle, TARGET_DENSITY)
                basis = "model_below_sheet" if throttle < first_sheet_throttle-1e-9 else "sheet_fit"
                writer.writerow([f"{voltage:.1f}", f"{throttle:.2f}", f"{point['thrust_n']:.3f}",
                                 f"{point['current_a']:.3f}", f"{point['power_w']:.2f}", f"{point['rpm']:.1f}",
                                 f"{point['torque_nm']:.4f}", basis])
    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    record = {
        "source": SOURCE.name, "source_sha256": source_hash,
        "curve": curve_path.name, "curve_sha256": hashlib.sha256(curve_path.read_bytes()).hexdigest(),
        "target_density_kg_m3": TARGET_DENSITY, "test_density_kg_m3": model["test_density"],
        "voltages_v": list(VOLTAGES), "throttles": list(THROTTLES),
        "equations": __doc__.split("\n\n")[1].strip(),
        "kv_rpm_per_v": model["kv"], "esc_efficiency_assumed": ESC_EFFICIENCY,
        "equivalent_resistance_ohm": model["resistance"], "no_load_current_a": model["idle"],
        "ct_vs_rpm_polyfit": model["ct_line"], "cq_vs_rpm_polyfit": model["cq_line"],
        "coefficient_rpm_range": model["rpm_range"],
        "max_abs_error_vs_sheet": summary, "validation_rows": rows,
        "limitations": [
            "One manufacturer sweep at 48 V; other voltages come from the motor model.",
            "Static thrust only: no inflow, ground effect, rotor-speed dynamics or thermal limits.",
            "Rows below 30 % throttle are model extrapolation (basis model_below_sheet).",
            "Sheet is for the V8013 IPE edition; installed PRO edition to be confirmed.",
            "No installed bench test yet; the curve is a provisional surrogate.",
        ],
    }
    fit_path = args.output_dir/"mad_v8013_150_fluxer26_fit.json"
    fit_path.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps({key: record[key] for key in ("equivalent_resistance_ohm", "no_load_current_a",
                                                    "max_abs_error_vs_sheet")}, indent=2))


if __name__ == "__main__":
    main()
