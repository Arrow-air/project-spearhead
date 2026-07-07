"""Export adapters for GUI code."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import TextIO

from ...export import export_csv as _export_simulation_csv
from ...export import export_json as _export_simulation_json
from ...result import SimulationResult
from ...stability.export import export_json as _export_stability_json
from ...stability.export import export_modes_csv as _export_stability_modes_csv
from ...stability.report import export_markdown as _export_stability_markdown
from ...stability.types import CGSweepResult, StabilityResult


def export_simulation_csv(result: SimulationResult, destination: str | Path | TextIO) -> None:
    """Delegate simulation CSV export to the shared exporter."""
    _export_simulation_csv(result, destination)


def export_simulation_json(result: SimulationResult, destination: str | Path | TextIO) -> None:
    """Delegate simulation JSON export to the shared exporter."""
    _export_simulation_json(result, destination)


def export_stability_json(payload: StabilityResult | CGSweepResult, destination: str | Path | TextIO) -> None:
    """Delegate stability JSON export to the shared exporter."""
    _export_stability_json(payload, destination)


def export_stability_modes_csv(
    payload: StabilityResult | CGSweepResult,
    destination: str | Path | TextIO,
) -> None:
    """Delegate stability mode CSV export to the shared exporter."""
    _export_stability_modes_csv(payload, destination)


def export_stability_markdown(result: StabilityResult, destination: str | Path | TextIO) -> None:
    """Delegate stability Markdown export to the shared exporter."""
    _export_stability_markdown(result, destination)


def export_simulation_csv_text(result: SimulationResult) -> str:
    """Return simulation CSV export text for browser download controls."""
    stream = StringIO()
    export_simulation_csv(result, stream)
    return stream.getvalue()


def export_simulation_json_text(result: SimulationResult) -> str:
    """Return simulation JSON export text for browser download controls."""
    stream = StringIO()
    export_simulation_json(result, stream)
    return stream.getvalue()


def export_stability_json_text(payload: StabilityResult | CGSweepResult) -> str:
    """Return stability JSON export text for browser download controls."""
    stream = StringIO()
    export_stability_json(payload, stream)
    return stream.getvalue()


def export_stability_modes_csv_text(payload: StabilityResult | CGSweepResult) -> str:
    """Return stability mode CSV export text for browser download controls."""
    stream = StringIO()
    export_stability_modes_csv(payload, stream)
    return stream.getvalue()


def export_stability_markdown_text(result: StabilityResult) -> str:
    """Return stability Markdown export text for browser download controls."""
    stream = StringIO()
    export_stability_markdown(result, stream)
    return stream.getvalue()

