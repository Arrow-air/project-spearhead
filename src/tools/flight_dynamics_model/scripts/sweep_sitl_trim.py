"""Sweep six-axis trims using the same new AeroDB and CG as SITL.

Run from the repository with PYTHONPATH=src/tools/flight_dynamics_model.
Mass variations hold inertia and CG fixed: they are sensitivity cases only.
"""

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from spearhead.dynamics import dynamics
from spearhead.force_moment import compute_force_moment_breakdown
from spearhead.sitl.backend import CONFIG_ROOT, load_setup
from spearhead.sitl.physics import full_trim
from spearhead.trim import trim_fixed_wing_longitudinal


def local_modes(state, control, params, step):
    # Omit position and heading symmetries; retain coupled body motion.
    indices = [3, 4, 5, 6, 7, 9, 10, 11]
    columns = []
    for index in indices:
        plus, minus = state.copy(), state.copy()
        plus[index] += step
        minus[index] -= step
        columns.append((dynamics(0, plus, params, control) -
                        dynamics(0, minus, params, control))[indices] / (2*step))
    return np.linalg.eigvals(np.column_stack(columns))


def run(output):
    physics, _, manifest = load_setup(CONFIG_ROOT / "fixed_wing_control.yaml")
    rows, details = [], []
    for mass in (22.5, 25., 27.5):
        params = replace(physics.params, mass=mass)
        for speed in (22., 24., 26., 28., 30., 32., 34., 36.):
            row = {"mass_kg": mass, "speed_m_s": speed, "trim_found": False}
            try:
                state, control = full_trim(params, speed, 100.)
            except ValueError as exc:
                row["reason"] = str(exc)
                longitudinal = trim_fixed_wing_longitudinal(params, speed, 100.)
                if longitudinal["success"]:
                    # Elevator data extend to +5 deg alpha, but lateral control
                    # data stop at zero. This is diagnostic, not a six-axis trim.
                    row["longitudinal_only_alpha_deg"] = float(longitudinal["alpha_deg"])
                    row["longitudinal_only_elevator_deg"] = float(longitudinal["de_deg"])
            else:
                loads = compute_force_moment_breakdown(0, state, params, control)["total"]
                eigenvalues = local_modes(state, control, params, 1e-5)
                finer = local_modes(state, control, params, 5e-6)
                maximum = float(max(eigenvalues.real))
                row.update(trim_found=True, reason="",
                           alpha_deg=float(np.rad2deg(np.arctan2(state[5], state[3]))),
                           pitch_deg=float(np.rad2deg(state[7])),
                           roll_deg=float(np.rad2deg(state[6])),
                           elevator_deg=float(np.rad2deg(control.de)),
                           aileron_deg=float(np.rad2deg(control.da)),
                           rudder_deg=float(np.rad2deg(control.dr)),
                           throttle_pct=float(100*control.throttle),
                           surface_margin_deg=float(5-max(abs(np.rad2deg(
                               [control.de, control.da, control.dr])))),
                           force_residual_N=float(np.linalg.norm(loads["force_b"])),
                           moment_residual_Nm=float(np.linalg.norm(loads["moment_b"])),
                           vertical_speed_m_s=float(dynamics(0, state, params, control)[2]),
                           max_real_eigenvalue_per_s=maximum,
                           eigenvalue_step_check_per_s=abs(maximum-float(max(finer.real))),
                           local_model_stable=maximum < -1e-6)
                details.append({"mass_kg": mass, "speed_m_s": speed,
                                "state": state.tolist(), "control": asdict(control),
                                "eigenvalues": [[float(v.real), float(v.imag)] for v in eigenvalues]})
            rows.append(row)
            print(json.dumps(row), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (output / "trim_sweep.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (output / "trim_sweep.json").write_text(json.dumps({
        "baseline": manifest, "altitude_m": 100,
        "limitations": ["Diagnostic CG, provisional inertia/rate derivatives/pusher",
                        "Mass sensitivity holds CG and inertia fixed",
                        "Measured control alpha range [-5,0] deg; beta=0 only",
                        "Lateral linearization assumes control increments constant versus beta",
                        "Trim feasibility and local model eigenvalues are not flight validation",
                        "No ArduPilot gains or sensor noise parameters changed"],
        "rows": rows, "trim_states_and_modes": details}, indent=2)+"\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for mass in (22.5, 25., 27.5):
        valid = [r for r in rows if r["mass_kg"] == mass and r["trim_found"]]
        speeds = [r["speed_m_s"] for r in valid]
        for axis, key, label in zip(axes.flat,
                ["alpha_deg", "elevator_deg", "throttle_pct", "max_real_eigenvalue_per_s"],
                ["Alpha (deg)", "Elevator (deg)", "Throttle (%)", "Largest Re(eigenvalue) (1/s)"]):
            axis.plot(speeds, [r[key] for r in valid], "o-", label=f"{mass:g} kg")
            axis.set(xlabel="TAS (m/s)", ylabel=label)
            axis.grid(alpha=.25)
    axes[0, 0].axhline(0, color="gray", linestyle=":")
    axes[0, 1].axhline(-5, color="gray", linestyle=":")
    axes[1, 1].axhline(0, color="gray", linestyle=":")
    axes[0, 0].legend()
    figure.suptitle("New AeroDB · six-axis level trims · diagnostic CG fixed\n"
                    "Only feasible trims plotted; local modes use provisional dynamics")
    figure.savefig(output / "trim_sweep.png", dpi=160)
    plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
