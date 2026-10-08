"""Replay recorded PWM through the common FDM without running firmware or GUI."""

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np
import yaml

from .backend import load_setup, sha256


def replay(directory, *, max_seconds=None):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    # Reuse recorded config, not a possibly edited current default.
    with tempfile.TemporaryDirectory(prefix="spearhead-replay-") as temporary:
        path = Path(temporary) / "setup.yaml"
        path.write_text(yaml.safe_dump(manifest["config"]))
        physics, mapping, current = load_setup(path)
    for name in ["database_sha256", "deflections_sha256"]:
        if current[name] != manifest[name]:
            raise ValueError(f"Replay source data mismatch: {name}")
    # Source mismatches must be visible even if the trajectory still agrees.
    source_match = current["fdm_source_sha256"] == manifest["fdm_source_sha256"]
    physical_modules = ["aero.py", "aerodb.py", "aircraft.py", "controls.py", "dynamics.py",
                        "force_moment.py", "params.py", "propulsion.py", "rotations.py",
                        "sitl/physics.py", "sitl/actuators.py"]
    physics_source_match = all(current["fdm_source_sha256"][name] ==
                               manifest["fdm_source_sha256"].get(name) for name in physical_modules)
    first = previous = None
    maxima = np.zeros(3)
    samples = 0
    with (directory / "truth.jsonl").open() as stream:
        for line in stream:
            row = json.loads(line)
            if first is None:
                # Seed only the first recorded state; all later states come
                # from recorded PWM, except the explicit bootstrap air-start.
                first = row
                physics.state = np.array(row["state"])
                physics.time = row["time_s"]
                previous = row
                continue
            elapsed = row["time_s"]-first["time_s"]
            if max_seconds is not None and elapsed > max_seconds+1e-9:
                break
            dt = row["time_s"]-previous["time_s"]
            if row["frame"] != (previous["frame"]+1) & 0xffffffff or not np.isclose(
                    dt, 1/row["physics_rate_hz"], rtol=0, atol=1e-9):
                raise ValueError("Replay requires contiguous lockstep records")
            if row.get("bootstrap_air_start"):
                if not previous["hold"]:
                    raise ValueError("Unexpected air-start in free flight")
                physics.state = np.array(manifest["bootstrap_air_start_state"])
            if not previous["hold"] and row["hold"]:
                raise ValueError("Hold cannot resume after release")
            physics.hold = row["hold"]
            # Log timestamps accumulate floating-point roundoff; the original
            # runner always advances by the exact reciprocal configured rate.
            physics.advance(1/row["physics_rate_hz"], mapping.control(row["pwm"]))
            delta = physics.state-np.array(row["state"])
            angle_delta = (delta[6:9]+np.pi) % (2*np.pi)-np.pi
            maxima = np.maximum(maxima, [np.max(abs(delta[:3])), np.max(abs(delta[3:6])),
                                          np.rad2deg(np.max(abs(angle_delta)))])
            samples += 1
            previous = row
    if samples == 0:
        raise ValueError("No replayable physics frames")
    limits = np.array([.01, .001, .01])
    return {"passed": bool(np.all(maxima < limits)), "source_hashes_match": source_match,
            "physics_source_hashes_match": physics_source_match,
            "replayed_frames": samples, "replayed_s": previous["time_s"]-first["time_s"],
            "seed": "first truth frame only; explicit bootstrap air-start from manifest",
            "max_position_error_m": maxima[0], "max_velocity_error_m_s": maxima[1],
            "max_attitude_error_deg": maxima[2],
            "truth_sha256": sha256(directory / "truth.jsonl"),
            "thresholds": dict(zip(["position_m", "velocity_m_s", "attitude_deg"], limits))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-seconds", type=float)
    args = parser.parse_args()
    if args.max_seconds is not None and (not np.isfinite(args.max_seconds) or args.max_seconds <= 0):
        parser.error("max-seconds must be positive and finite")
    if not (args.directory / "report.json").is_file():
        parser.error("Replay requires a finished run with report.json")
    result = replay(args.directory, max_seconds=args.max_seconds)
    content = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(content)
    print(content)
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
