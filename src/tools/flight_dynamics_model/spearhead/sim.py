"""Command-line entry point for shared Spearhead simulations."""

from __future__ import annotations

import argparse
from pathlib import Path

from .export import export_csv, export_json
from .scenarios import load_scenario
from .simulation import run_simulation
from .aircraft import aircraft_options, aerodynamic_options, load_aircraft_config


def _run(args: argparse.Namespace) -> int:
    params = None if args.aircraft is None else load_aircraft_config(args.aircraft)
    config = load_scenario(args.scenario, params=params, aerodynamics=args.aerodynamics)
    result = run_simulation(config)
    if args.csv:
        export_csv(result, args.csv)
    if args.json:
        export_json(result, args.json)

    print(f"{config.name} [{config.params.aircraft_id}]: success={result.success} samples={result.time.size}")
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
    run_parser.add_argument("--aircraft", choices=list(aircraft_options()), help="override the scenario aircraft profile")
    run_parser.add_argument("--aerodynamics", choices=list(aerodynamic_options()), help="override the aerodynamic dataset independently of hardware")
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
