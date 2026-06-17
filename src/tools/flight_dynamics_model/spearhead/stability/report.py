"""Markdown reporting for preliminary stability analysis."""

from __future__ import annotations

from pathlib import Path
from typing import TextIO

from .types import StabilityResult


def stability_markdown(result: StabilityResult) -> str:
    """Return a compact Markdown stability report."""
    trim = result.trim_result
    static = result.static_derivatives.derivatives
    lines = [
        f"# Stability Analysis: {result.config.name}",
        "",
        f"> {result.warning}",
        "",
        "## Method",
        "",
        "This analysis uses `SimulationConfig` and the shared FDM backend. Trim is obtained through",
        "`run_simulation(config)`, and linearization calls the existing `dynamics()` function with",
        "the existing force/moment, aero, propulsion, and control paths.",
        "",
        "## Limitations",
        "",
        "- This PR restores source for the stability workflow but does not add GUI/dashboard features.",
        "- The current tracked aerodynamic database is the nominal ADB only.",
        "- CG sweep support uses runtime `cg_body_xyz` moment shifting; it does not recreate missing",
        "  Nondimit per-CG CSV cases.",
        "- Inertia, rate damping, control derivatives, and propulsion remain preliminary placeholders.",
        "",
        "## Trim",
        "",
        f"- Success: `{bool(trim['success'])}`",
        f"- Message: `{trim['message']}`",
        f"- Alpha: `{float(trim['alpha_deg']):+.4f} deg`",
        f"- Theta: `{float(trim['theta_deg']):+.4f} deg`",
        f"- Elevator: `{float(trim['de_deg']):+.4f} deg`",
        f"- Throttle: `{float(trim['throttle']):.5f}`",
        "",
        "## Static Derivatives",
        "",
    ]
    for key in sorted(static):
        lines.append(f"- `{key}`: `{static[key]:+.6f}`")

    lines.extend(["", "## Modes", ""])
    for analysis in (result.longitudinal, result.lateral_directional):
        lines.extend(
            [
                f"### {analysis.subsystem.replace('_', ' ').title()}",
                "",
                "| Mode | Eigenvalue | wn rad/s | zeta | period s | half s | double s | Classification |",
                "|---|---:|---:|---:|---:|---:|---:|---|",
            ]
        )
        for mode in analysis.modes:
            eigenvalue = f"{mode.eigenvalue.real:+.6g}{mode.eigenvalue.imag:+.6g}j"
            damping = "n/a" if mode.damping_ratio is None else f"{mode.damping_ratio:.6g}"
            period = "n/a" if mode.oscillation_period_s is None else f"{mode.oscillation_period_s:.6g}"
            half = "n/a" if mode.time_to_half_s is None else f"{mode.time_to_half_s:.6g}"
            double = "n/a" if mode.time_to_double_s is None else f"{mode.time_to_double_s:.6g}"
            lines.append(
                f"| {mode.mode} | `{eigenvalue}` | {mode.natural_frequency_rad_s:.6g} | "
                f"{damping} | {period} | {half} | {double} | {mode.classification} |"
            )
        lines.append("")

    return "\n".join(lines)


def export_markdown(result: StabilityResult, destination: str | Path | TextIO) -> None:
    """Write a Markdown stability report."""
    text = stability_markdown(result)
    if hasattr(destination, "write"):
        destination.write(text)
        destination.write("\n")
        return
    Path(destination).write_text(text + "\n")
