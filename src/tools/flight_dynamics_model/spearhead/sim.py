"""Command-line entry point for shared Spearhead simulations."""

from __future__ import annotations

import argparse
from pathlib import Path

from .export import export_csv, export_json
from .scenarios import load_scenario
from .simulation import run_simulation


def _run(args: argparse.Namespace) -> int:
    config = load_scenario(args.scenario)
    result = run_simulation(config)
    if args.csv:
        export_csv(result, args.csv)
    if args.json:
        export_json(result, args.json)

    print(f"{config.name}: success={result.success} samples={result.time.size}")
    if not result.success:
        print(result.message)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(prog="python -m spearhead.sim")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="run a YAML simulation scenario")
    run_parser.add_argument("scenario", type=Path)
    run_parser.add_argument("--csv", type=Path, help="write CSV time-history output")
    run_parser.add_argument("--json", type=Path, help="write JSON time-history output")
    run_parser.set_defaults(func=_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the Spearhead simulation CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
