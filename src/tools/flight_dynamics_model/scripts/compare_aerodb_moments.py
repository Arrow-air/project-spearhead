"""Overlay old/new body moments at a common density, speed and moment center.

Run from the repository root with PYTHONPATH=src/tools/flight_dynamics_model.
These are zero-control database slices; no provisional controls or damping.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from spearhead.aircraft import load_aircraft_config
from spearhead.aerodb import dimensional_body_force_moment, load_aerodynamic_database


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rho, speed = 1.021, 25.0
    center = np.zeros(3)
    qbar = 0.5 * rho * speed**2
    databases = {
        "old": load_aerodynamic_database(load_aircraft_config("reference_design", aerodynamics="adb_v1_1").aero_database_path),
        "new": load_aerodynamic_database(load_aircraft_config("reference_design", aerodynamics="tail_update_2026_10_07").aero_database_path),
    }

    def moments(database, alpha, beta):
        _, moment = dimensional_body_force_moment(
            database.coefficients_at(float(alpha), float(beta)), qbar,
            database.reference, cg_body_xyz=center,
        )
        return moment

    angles = np.arange(-90.0, 90.1, 1.0)
    slices = {}
    rows = []
    for axis in ("alpha", "beta"):
        slices[axis] = {}
        for version, database in databases.items():
            values = np.array([moments(database, a if axis == "alpha" else 0,
                                       a if axis == "beta" else 0) for a in angles])
            slices[axis][version] = values
            for angle, vector in zip(angles, values):
                rows.append({"sweep": axis, "version": version, "angle_deg": angle,
                             "alpha_deg": angle if axis == "alpha" else 0,
                             "beta_deg": angle if axis == "beta" else 0,
                             "Mx_Nm": vector[0], "My_Nm": vector[1], "Mz_Nm": vector[2]})
    with (args.output_dir / "moment_slices.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # Finite differences use the same origin and pressure as the plots.
    summary = {"rho_kg_m3": rho, "speed_m_s": speed, "moment_center_body_xyz_m": center.tolist(),
               "controls": "zero", "sources": {key: str(load_aircraft_config("reference_design", aerodynamics=dataset).aero_database_path)
                  for key, dataset in [("old", "adb_v1_1"), ("new", "tail_update_2026_10_07")]}, "datasets": {}}
    for version, database in databases.items():
        pitch_points = {str(alpha): float(moments(database, alpha, 0)[1]) for alpha in [-5, 0, 5]}
        alpha_derivative = (moments(database, 1, 0) - moments(database, -1, 0)) / np.deg2rad(2)
        beta_derivative = (moments(database, 0, 1) - moments(database, 0, -1)) / np.deg2rad(2)
        summary["datasets"][version] = {
            "pitch_moment_at_alpha_Nm": pitch_points,
            "dMy_dalpha_Nm_per_rad_at_origin": float(alpha_derivative[1]),
            "dMx_dbeta_Nm_per_rad_at_origin": float(beta_derivative[0]),
            "dMz_dbeta_Nm_per_rad_at_origin": float(beta_derivative[2]),
        }
    (args.output_dir / "moment_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "axes.labelcolor": "#334155",
                         "xtick.color": "#475569", "ytick.color": "#475569"})
    panels = [("beta", 0, "Roll moment · Mx", "β (°) · α = 0°", (-15, 15)),
              ("alpha", 1, "Pitch moment · My", "α (°) · β = 0°", (-10, 20)),
              ("beta", 2, "Yaw moment · Mz", "β (°) · α = 0°", (-15, 15))]
    colors = {"old": "#2563eb", "new": "#ea580c"}
    labels = {"old": "Previous · ADB v1.1", "new": "New · October 2026"}
    for name, full in [("moment_overlay", False), ("moment_overlay_full_range", True)]:
        figure, axes = plt.subplots(1, 3, figsize=(15, 5.6))
        handles = []
        for ax, (axis, component, title, xlabel, limits) in zip(axes, panels):
            bounds = (-90, 90) if full else limits
            mask = (angles >= bounds[0]) & (angles <= bounds[1])
            for version in ("old", "new"):
                line, = ax.plot(angles[mask], slices[axis][version][mask, component],
                                color=colors[version], linewidth=2.3, label=labels[version],
                                linestyle="--" if version == "old" else "-")
                if len(handles) < 2:
                    handles.append(line)
            ax.axhline(0, color="#64748b", linewidth=0.8)
            ax.axvline(0, color="#cbd5e1", linewidth=0.8)
            ax.set(title=title, xlabel=xlabel, ylabel="Moment (N·m)", xlim=bounds)
            ax.grid(True, color="#e2e8f0", linewidth=0.6)
            ax.set_axisbelow(True)
            if component == 1 and not full:
                for version in ("old", "new"):
                    value = summary["datasets"][version]["pitch_moment_at_alpha_Nm"]["0"]
                    ax.scatter([0], [value], color=colors[version], s=35, zorder=5)
                    ax.annotate(f"{value:.1f} N·m", (0, value), xytext=(10, 8),
                                textcoords="offset points", color=colors[version], fontsize=10)
        figure.suptitle("Previous and new AeroDB — moments (−90°..+90°)" if full else "Previous and new AeroDB — base moment comparison",
                        fontsize=18, fontweight="bold", x=0.055, ha="left", y=0.98, color="#0f172a")
        figure.text(0.055, 0.895, "Common moment center (0, 0, 0) m  |  ρ = 1.021 kg/m³  |  V = 25 m/s  |  Controls = 0",
                    fontsize=11, color="#475569")
        figure.legend(handles, [labels[key] for key in ("old", "new")],
                      loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=2, frameon=False)
        figure.text(0.055, 0.02,
                     "Previous moments transferred from their own reference to the origin. New data use the confirmed origin; the raw header mismatch remains open.",
                     fontsize=9, color="#64748b")
        figure.subplots_adjust(left=0.07, right=0.98, top=0.78, bottom=0.23, wspace=0.30)
        figure.savefig(args.output_dir / f"{name}.png", dpi=180, facecolor="white")
        figure.savefig(args.output_dir / f"{name}.pdf", facecolor="white")
        plt.close(figure)
    print(json.dumps(summary["datasets"], indent=2))


if __name__ == "__main__":
    main()
