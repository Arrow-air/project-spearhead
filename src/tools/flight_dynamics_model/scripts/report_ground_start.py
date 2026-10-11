"""Regenerate the ground-start validation summary and flight traces."""

import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spearhead.sitl.ground_flight import analyse_ground


def main():
    root = Path(__file__).resolve().parents[1] / "artifacts/ground_start_2026_10_11"
    summary = {}
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for height, name in ((5, "verified_5m"), (10, "verified_10m"), (20, "verified_20m_ui")):
        run = root / name
        report = json.loads((run / "report.json").read_text())
        events = json.loads((run / "events.json").read_text())
        result = analyse_ground(run / "truth.jsonl", height, 15, report["normal_arm"],
                                report["mission_accepted"], report["firmware_disarmed"], events)
        result["run_directory"] = name
        summary[str(height)] = result
        t, states, contact = [], [], []
        with (run / "truth.jsonl").open() as stream:
            for index, line in enumerate(stream):
                if index % 10:
                    continue
                row = json.loads(line)
                t.append(row["time_s"]-report["flight_start_s"])
                states.append(row["state"])
                contact.append(row["ground"]["normal_force_n"])
        t, states = np.array(t), np.array(states)
        axes[0].plot(t, -states[:, 2], label=f"{height} m target")
        axes[1].plot(t, np.rad2deg(np.linalg.norm(states[:, 6:8], axis=1)))
        axes[2].plot(t, contact)
    axes[0].set_ylabel("CG height above ground (m)"); axes[0].legend()
    axes[1].set_ylabel("Roll/pitch norm (deg)")
    axes[2].set_ylabel("Gear reaction (N)"); axes[2].set_xlabel("Time since takeoff command (s)")
    for ax in axes:
        ax.grid(alpha=.25); ax.set_xlim(-5, 75)
    fig.suptitle("PT1 21 kg · ArduPlane 4.7.1 · Ground → takeoff → hover → landing")
    fig.tight_layout(); fig.savefig(root / "flight_traces.png", dpi=150)
    (root / "validation_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    lines = ["# PT1 ground-start integration — 11 October 2026", "",
             "Implemented in the Flight Deck: **PT1 → Ground · VTOL takeoff → Takeoff → hover → land**.",
             "Set the takeoff altitude and hover dwell independently; climb and landing add to that time.", "",
             "All three real ArduPlane 4.7.1 runs normally armed, flew an accepted AUTO mission, physically contacted the ground and disarmed.",
             "The pusher stayed off throughout these verified runs. The 20 m run was launched using the browser controls.", "",
             "| Target | Maximum CG height | Settled hover interval | Mission hover interval | Maximum tilt | Result |",
             "| --- | ---: | ---: | ---: | ---: | --- |"]
    for height, result in summary.items():
        lines.append(f"| {height} m | {result['altitude_range_m'][1]:.3f} m | {result['hover_duration_s']:.2f} s | "
                     f"{result['mission_hover_duration_s']:.2f} s | {result['max_tilt_deg']:.3f}° | {'PASS' if result['passed'] else 'FAIL'} |")
    lines += ["", "![Flight traces](flight_traces.png)", "",
              "Settled hover means a contiguous interval within ±1 m of the target and body vertical speed below 0.5 m/s. "
              "The mission hover interval comes from firmware mission-sequence telemetry. The first seconds after climb allow settling; "
              "acceptance requires at least the requested mission dwell and at least dwell minus 3 seconds of settled hover.", "",
              "The gear supports 205.94 N at rest. Its assumed four-foot footprint, compliance and friction feed both rigid-body motion "
              "and accelerometer specific force. Initial/resting CG height is approximately 0.297 m; there is no hold or imposed flight trajectory.", "",
              "**Scope and assumptions:** PT1 mass is 21 kg, density is 1.021, and the October aerodynamic data remain selected. "
              "Landing gear, CG/inertia, motor positions and vertical-flight aero are provisional. No terrain variation, rotor ground effect, "
              "propwash or battery sag is represented. This validates simulator integration, not the real aircraft or forward transition.", "",
              "The initial development trial (`trial_5m_b`) landed but exposed a short pusher pulse at loiter entry. "
              "The verified mission adds DO_VTOL_TRANSITION(MC) before loiter; all verified pusher outputs are zero. "
              "The 10 m run uses the final runner, including navigation-validity monitoring and terrain disabled for the flat-ground scenario.", "",
              "Startup messages about AHRS/pre-arm checks are preserved: normal arming waits for EKF3 navigation readiness. "
              "The runner does not bypass checks or force-disarm. Raw telemetry, actuator outputs, FDM truth, manifests and firmware logs "
              "are retained in each run folder. Original reports remain unchanged; `validation_summary.json` applies the final analysis consistently to all runs.", "",
              "See [ground-start setup and parameter rationale](../../docs/ground_start.md) for the configuration and official sources."]
    (root / "report.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({height: result["passed"] for height, result in summary.items()}))


if __name__ == "__main__":
    main()
