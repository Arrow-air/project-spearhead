"""Compare run_open_loop numerical compatibility across two worktrees."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


CASES = {
    "default_trimmed": "default_trimmed",
    "untrimmed_schedule": "untrimmed_schedule",
    "constant_control": "constant_control",
}


RUNNER = r"""
import sys
from pathlib import Path

import numpy as np

root = Path(sys.argv[1])
case = sys.argv[2]
output = Path(sys.argv[3])
sys.path.insert(0, str(root / "src/tools/flight_dynamics_model"))

from spearhead.controls import Control
from spearhead.force_moment import compute_force_moment_history
from spearhead.params import AircraftParams
from spearhead.simulation import default_initial_state, run_open_loop

params = AircraftParams()
if case == "default_trimmed":
    sol = run_open_loop(t_final=1.0, n_points=11, params=params)
elif case == "untrimmed_schedule":
    sol = run_open_loop(t_final=1.0, n_points=11, params=params, start_from_trim=False)
elif case == "constant_control":
    sol = run_open_loop(
        t_final=0.5,
        n_points=6,
        params=params,
        initial_state=default_initial_state(),
        control_override=Control(de=0.01, da=-0.02, dr=0.03, throttle=0.5),
        start_from_trim=False,
    )
else:
    raise ValueError(f"Unknown case: {case}")

history = compute_force_moment_history(
    sol.t,
    sol.y,
    params,
    control_override=getattr(sol, "control_override", None),
)
np.savez(
    output,
    t=sol.t,
    y=sol.y,
    success=np.array([bool(sol.success)]),
    status=np.array([int(sol.status)]),
    message=np.array([str(sol.message)]),
    total_force_b=history["total"]["force_b"],
    total_moment_b=history["total"]["moment_b"],
)
"""


def _run_case(root: Path, case: str, output: Path) -> None:
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    subprocess.run(
        [sys.executable, "-c", RUNNER, str(root), case, str(output)],
        check=True,
        env=env,
    )


def _max_relative_difference(expected: np.ndarray, actual: np.ndarray) -> float:
    scale = np.maximum(np.abs(expected), 1e-300)
    return float(np.max(np.abs(actual - expected) / scale))


def _compare_array(case: str, name: str, expected: np.ndarray, actual: np.ndarray) -> None:
    try:
        if name == "t":
            np.testing.assert_array_equal(actual, expected)
        else:
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    except AssertionError as exc:
        max_abs = float(np.max(np.abs(actual - expected)))
        max_rel = _max_relative_difference(expected, actual)
        raise AssertionError(
            f"{case} {name} mismatch: max_abs={max_abs:.17e}, max_rel={max_rel:.17e}"
        ) from exc


def compare_worktrees(main_root: Path, feature_root: Path) -> None:
    """Run and compare representative run_open_loop cases."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for case in CASES:
            main_output = tmp_path / f"main_{case}.npz"
            feature_output = tmp_path / f"feature_{case}.npz"
            _run_case(main_root, case, main_output)
            _run_case(feature_root, case, feature_output)

            main = np.load(main_output)
            feature = np.load(feature_output)
            for scalar_name in ("success", "status", "message"):
                if feature[scalar_name][0] != main[scalar_name][0]:
                    raise AssertionError(
                        f"{case} {scalar_name} mismatch: "
                        f"main={main[scalar_name][0]!r}, feature={feature[scalar_name][0]!r}"
                    )
            for array_name in ("t", "y", "total_force_b", "total_moment_b"):
                _compare_array(case, array_name, main[array_name], feature[array_name])

            print(
                f"{case}: identical t/status/message; "
                "y and force/moment histories allclose at rtol=1e-12 atol=1e-12"
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-root", required=True, type=Path)
    parser.add_argument("--feature-root", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    compare_worktrees(args.main_root, args.feature_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
