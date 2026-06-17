"""CLI for preliminary stability analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ..scenarios import load_scenario
from .analysis import analyze_stability
from .export import export_json, export_modes_csv
from .report import export_markdown
from .sweep import run_cg_sweep


def _write_outputs(result, args: argparse.Namespace) -> None:
    if args.json:
        export_json(result, args.json)
    if args.csv:
        export_modes_csv(result, args.csv)
    if getattr(args, "markdown", None):
        if hasattr(result, "points"):
            successful = next((point.result for point in result.points if point.result is not None), None)
            if successful is not None:
                export_markdown(successful, args.markdown)
        else:
            export_markdown(result, args.markdown)


def _analyze(args: argparse.Namespace) -> int:
    config = load_scenario(args.scenario)
    result = analyze_stability(config)
    _write_outputs(result, args)
    print(
        f"{config.name}: stability analysis complete; "
        f"longitudinal_modes={len(result.longitudinal.modes)} "
        f"lateral_directional_modes={len(result.lateral_directional.modes)}"
    )
    return 0


def _cg_sweep(args: argparse.Namespace) -> int:
    config = load_scenario(args.scenario)
    y_ref = float(config.params.cg_body_xyz[1])
    z_ref = float(config.params.cg_body_xyz[2])
    cg_positions = [np.array([x_ref, y_ref, z_ref], dtype=float) for x_ref in args.cg_x]
    result = run_cg_sweep(config, cg_positions)
    _write_outputs(result, args)
    print(
        f"{config.name}: cg sweep complete; "
        f"successful={len(result.successful_points)} total={len(result.points)}"
    )
    return 0 if len(result.successful_points) == len(result.points) else 1


def _add_output_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--json", type=Path, help="write JSON summary")
    parser.add_argument("--csv", type=Path, help="write CSV eigenvalue table")
    parser.add_argument("--markdown", type=Path, help="write Markdown summary")


def build_parser() -> argparse.ArgumentParser:
    """Build the stability CLI parser."""
    parser = argparse.ArgumentParser(prog="python -m spearhead.stability")
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze_parser = subparsers.add_parser("analyze", help="analyze one scenario")
    analyze_parser.add_argument("scenario", type=Path)
    _add_output_args(analyze_parser)
    analyze_parser.set_defaults(func=_analyze)

    sweep_parser = subparsers.add_parser("cg-sweep", help="analyze runtime CG x positions")
    sweep_parser.add_argument("scenario", type=Path)
    sweep_parser.add_argument("--cg-x", nargs="+", type=float, required=True)
    _add_output_args(sweep_parser)
    sweep_parser.set_defaults(func=_cg_sweep)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the stability CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
