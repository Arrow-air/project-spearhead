"""Reproduce the provisional PT1 pusher estimate and forward-trim checks.

Run with PYTHONPATH=src/tools/flight_dynamics_model; outputs are standalone
Markdown, JSON and a PNG plot. No reference/old aerodynamic database is used.
"""

import argparse
from dataclasses import replace
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from spearhead import load_aircraft_config
from spearhead.aircraft import AIRCRAFT_CONFIGS, aircraft_summary
from spearhead.controls import Control
from spearhead.force_moment import compute_force_moment_breakdown
from spearhead.propulsion import pusher_performance
from spearhead.pusher import estimate_pusher
from spearhead.sitl.physics import full_trim


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("artifacts/pt1_pusher_2026_10_11"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    pt1 = load_aircraft_config("pt15")
    system = pt1.pusher_system
    speeds = [0, 10, 20, 22, 24, 26, 28]
    grid = []
    envelopes = []
    for voltage in (42., 44.4, 50.4):
        for speed in speeds:
            def point(duty):
                return estimate_pusher(system, voltage, duty, speed, pt1.rho)
            maximum = point(1)
            sustainable_duty = (brentq(lambda duty: point(duty)["current_a"]-system.sustainable_current_a, .01, 1)
                                if maximum["above_sustainable_current"] else 1.)
            envelopes.append(dict(voltage_v=voltage, speed_m_s=speed, full_throttle=maximum,
                                  sustainable_current_guideline=point(sustainable_duty)))
        for speed in np.linspace(0, 28, 29):
            for throttle in (.5, .65, .8, .9, 1):
                grid.append(dict(voltage_v=voltage, speed_m_s=float(speed), throttle=throttle,
                                 **estimate_pusher(system, voltage, throttle, float(speed), pt1.rho)))

    # The original PT1 CG is preserved. The second case is explicitly a
    # numerical diagnostic, using the previously selected FDM CG/travel only.
    diagnostic = replace(pt1, cg_body_xyz=np.array([.013595893369954063, 0., .45]),
                         control_deflection_policy="linear", control_surface_limits_deg=(20, 20, 20))
    trims = []
    thrust_budgets = []
    for name, params in (("PT1 profile CG / measured ±5° surface range", pt1),
                         ("Diagnostic CG / assumed ±20° travel + control extrapolation", diagnostic)):
        for speed in (22, 24, 26, 28):
            row = dict(case=name, speed_m_s=speed, cg_body_xyz_m=params.cg_body_xyz.tolist())
            try:
                state, control = full_trim(params, speed, 100)
                loads = compute_force_moment_breakdown(0, state, params, control)
                propulsion = pusher_performance(control, params, state[3])
                row.update(balanced=True, state=state.tolist(), elevator_deg=float(np.rad2deg(control.de)),
                           throttle=control.throttle, pusher=propulsion,
                           force_residual_n=float(np.linalg.norm(loads["total"]["force_b"])),
                           moment_residual_nm=float(np.linalg.norm(loads["total"]["moment_b"])))
            except ValueError as exc:
                row.update(balanced=False, reason=str(exc))
            trims.append(row)
            print(name, speed, "balanced" if row["balanced"] else row["reason"], flush=True)

    # Solve the new AeroDB's required forward force independently of motor
    # saturation. The 120 N scale only parameterizes the unknown thrust here.
    for speed in (22, 24, 26, 28):
        state, control = full_trim(replace(diagnostic, pusher_system=None, Tmax=120), speed, 100)
        envelope = next(row for row in envelopes if row["voltage_v"] == 44.4 and row["speed_m_s"] == speed)
        thrust_budgets.append(dict(speed_m_s=speed, required_thrust_n=control.throttle*120,
                                  available_full_thrust_n=estimate_pusher(system, 44.4, 1, state[3], pt1.rho)["thrust_n"],
                                  available_guideline_thrust_n=envelope["sustainable_current_guideline"]["thrust_n"]))

    sources = json.loads((AIRCRAFT_CONFIGS.parent / "data/propulsion/pt1_pusher_source.json").read_text())
    result = dict(aircraft=aircraft_summary(pt1), source=sources, envelopes=envelopes,
                  performance_grid=grid, trim_checks=trims, diagnostic_thrust_budgets=thrust_budgets,
                  validation_scope="steady propulsion estimate and numerical force/moment balance; not transition or flight stability validation")
    (args.output / "performance.json").write_text(json.dumps(result, indent=2)+"\n")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for voltage, colour in ((44.4, "#0943bf"), (50.4, "#008c86")):
        for throttle, style in ((.8, "--"), (1., "-")):
            rows = [row for row in grid if row["voltage_v"] == voltage and row["throttle"] == throttle]
            x = [r["speed_m_s"] for r in rows]
            label = f"{voltage:g} V · {100*throttle:g}% duty"
            axes[0].plot(x, [r["thrust_n"] for r in rows], style, color=colour, label=label)
            axes[1].plot(x, [r["current_a"] for r in rows], style, color=colour, label=label)
    axes[1].axhline(24, color="#b45f06", linestyle=":", label="MAD sustainable guideline · 24 A")
    for ax, ylabel in zip(axes, ["Pusher thrust (N)", "Battery current (A)"]):
        ax.set(xlabel="Axial inflow speed (m/s)", ylabel=ylabel, xlim=(0, 28))
        ax.grid(alpha=.2)
        ax.legend(fontsize=8, loc="upper left" if ax is axes[0] else "lower left")
    fig.suptitle("PT1 · MAD V8013 150 KV + AMPX 80A + APC 20×12\nProvisional estimate · ρ = 1.021 kg/m³ · no installation loss or battery sag")
    fig.savefig(args.output / "performance.png", dpi=180)
    plt.close(fig)

    lines = ["# PT1 pusher estimate — 11 October 2026", "",
             "**Status: provisional analytical estimate, integrated into the PT1 FDM.** The user confirmed 150 KV. PT1 remains 21 kg with four KDE lift motors and the October AeroDB (moment centre 0,0,0; density 1.021 kg/m³). Component masses are not added again. The 25 kg reference profile is unchanged.", "",
             "## Estimated performance", "",
             "Full throttle is a short-run operating point, not a continuous rating. `Duty` means the assumed ESC voltage fraction, not measured throttle-to-RPM calibration.", "",
             "| Battery | Axial speed | Full-duty thrust | Full-duty input | Current | Thrust at ≤24 A | Duty at ≤24 A |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for row in envelopes:
        if row["voltage_v"] == 42:
            continue
        p, s = row["full_throttle"], row["sustainable_current_guideline"]
        lines.append(f"| {row['voltage_v']:g} V | {row['speed_m_s']:g} m/s | {p['thrust_n']:.1f} N | {p['power_w']:.0f} W | {p['current_a']:.1f} A | {s['thrust_n']:.1f} N | {100*s['effective_duty']:.1f}% |")
    lines += ["", "![Estimated thrust and current](performance.png)", "", "## Model and sources", "",
              "- [APC LP20012 Sport 20×12](https://www.apcprop.com/product/20x12/), not 20×12E or a MAD lift propeller. [APC's analytical coefficient map](https://www.apcprop.com/files/PER3_20x12.dat) is preserved with its v2022-0915 / 2022-09-22 header. APC calls these computational results, not bench measurements.",
              "- [Official MAD 150KV test sheet](https://cdn.shopify.com/s/files/1/0697/2293/1501/files/KV150-V8013IPE.webp?v=1722922466): Kv=150 RPM/V, I0=2.35 A at 30 V, phase resistance 37 mΩ. The sheet is labelled V8013 IPE; the requested product says PRO IPE. Confirm the installed edition. Its 26×7.8 thrust is never reused for the APC propeller.",
              "- Effective motor resistance is fitted to all 15 test-sheet RPM/torque/duty points: 0.105689 Ω, voltage-fit RMS 0.205 V. This incorporates commutation/ESC effects and differs from phase resistance. `Kt=60/(2πKv)`; `I_winding=Q/Kt+I0`; `duty·V=RPM/Kv+R_equivalent·I_winding`.",
              "- For each RPM, interpolate Ct and Cp versus advance ratio and between bracketing RPM sections. `J=V_axial/(nD)`, `T=Ct·ρ·n²·D⁴`, `P_shaft=Cp·ρ·n³·D⁵`, `Q=P_shaft/(2πn)`. Coefficients are scaled at the FDM density, rather than copying sea-level thrust numbers.",
              "- ESC efficiency is an assumed 97%. Battery current is `duty·I_winding/η_ESC`, rather than equating winding current with battery current. The steady solution is restricted to 77.8 A, 3569 W and a conservative 9500 RPM Sport propeller ceiling. These are model limits, not simulated ESC protection.",
              "- The manufacturer sheet says below 24 A is sustainable; 24–78 A is short term, approximately 10–30 s with temperature/cooling dependence. The FDM flags excursions above 24 A; it does not invent a thermal derating timeline. The ≤24 A table is an operating guideline, not a calibrated continuous endurance guarantee.",
              "- [MAD AMPX manual](https://ueeshop.ly200-cdn.com/u_file/UPAZ/UPAZ594/2407/04/file/AMPX80A5-14SV2Manual.pdf): 5–14S, 80 A with cooling, default 1050–1940 μs. PT1's simulator channel 3 and `SERVO3_MIN/MAX/TRIM` use 1050/1940/1050; `SERVO3_FUNCTION=70` remains throttle. No ArduPilot controller gains were changed. Actual aircraft calibration is still required.",
              "", "## Integration and limitations", "",
              "PT1's installed pusher now uses this estimate throughout the shared force/moment engine, trim solver, realtime telemetry and ArduPilot JSON backend. Nominal static Tmax is derived from the estimate; forward thrust uses the current body's axial velocity. KDE lift curves and commands remain independent.", "",
              "Assumptions: 12S 36–50.4 V modelling envelope; fixed nominal battery voltage during a run; constant no-load current; quasi-steady motor/propeller equilibrium; no installation efficiency reduction, fuselage/wake interference, crossflow correction, battery sag, temperature, motor dynamics or propwash. The thrust line provisionally passes through CG and pusher reaction torque is omitted because mount coordinates/rotation direction were not supplied. These omissions matter for trim and stability.", "",
              "Below APC's first RPM section, use its endpoint coefficients and flag outside-data operation. Beyond the available positive-thrust advance-ratio envelope, set Ct=0 and hold endpoint Cp; report the condition. Negative axial inflow is clipped to zero and flagged. Zero command yields zero thrust/power: stopped-prop drag, windmilling and regeneration are not modelled. The solver starts electric-pusher trim with positive thrust authority to avoid the flat zero-Ct branch.", "",
              "The deck still offers PT1 hover starts. A mode change can exercise the estimated pusher, but the full transition mission and ground contact are not implemented or validated by this deliverable.", "",
              "## Forward-trim checks", "",
              "Only the new AeroDB is used. The diagnostic case is a numerical check of thrust availability; it does not relocate the physical PT1 CG or establish a stable aircraft.", "",
              "| Case | Speed | Balanced | Elevator | Duty | Pusher current |",
              "|---|---:|---|---:|---:|---:|"]
    for row in trims:
        if row["balanced"]:
            lines.append(f"| {row['case']} | {row['speed_m_s']:g} m/s | Yes | {row['elevator_deg']:.2f}° | {100*row['throttle']:.1f}% | {row['pusher']['current_a']:.1f} A |")
        else:
            lines.append(f"| {row['case']} | {row['speed_m_s']:g} m/s | No | — | — | — |")
    lines += ["", "The new AeroDB's independently solved diagnostic forward-force demand is compared below with the estimated pusher envelope at nominal voltage. No older aerodynamic dataset is involved.", "",
              "| Speed | Required forward thrust | Full-duty available | ≤24 A available |",
              "|---|---:|---:|---:|"]
    for row in thrust_budgets:
        lines.append(f"| {row['speed_m_s']:g} m/s | {row['required_thrust_n']:.2f} N | {row['available_full_thrust_n']:.2f} N | {row['available_guideline_thrust_n']:.2f} N |")
    lines += ["", "At 26 m/s, the ≤24 A thrust margin is less than 1%; at 28 m/s, the nominal full-duty thrust differs from the demand by only about 0.06 N. Those differences are far smaller than the unquantified installation/model uncertainty. Treat 26–28 m/s as marginal in this estimate, rather than asserting a precise maximum flight speed. No numerical confidence interval can be established without installed tests."]
    lines += ["", "A balanced trim is not evidence of dynamic stability, transition reliability or safe real-flight operation. An unsuccessful profile-CG trim can reflect CG/control coverage as well as thrust limits. Obtain measured CG, inertia, mount/rotation direction and control limits before treating these estimates as an aircraft tune.", "",
              "## Reproduction and deliverables", "",
              "Run `PYTHONPATH=src/tools/flight_dynamics_model python src/tools/flight_dynamics_model/scripts/estimate_pt1_pusher.py`. `performance.json` contains the full 42/44.4/50.4 V grid and all trim results; `performance.png` plots thrust and current. Source files, motor transcription, fitted constants and SHA256 hashes are in `data/propulsion/pt1_pusher_source.json` and its referenced files. Automated validation and any SITL run are recorded separately in `validation.md`.", ""]
    (args.output / "report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
