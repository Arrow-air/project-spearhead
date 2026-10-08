"""CSV export for structured simulation results."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import TextIO

from ..result import SimulationResult
from ..state import STATE_COLUMNS


def export_csv(result: SimulationResult, destination: str | Path | TextIO) -> None:
    """Write a simulation result time history to CSV."""
    fieldnames = ("time",) + STATE_COLUMNS + ("aircraft_id", "aero_database_version")
    should_close = False
    if hasattr(destination, "write"):
        stream = destination
    else:
        stream = Path(destination).open("w", newline="")
        should_close = True

    try:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(result.time.size):
            writer.writerow({
                "aircraft_id": result.config.params.aircraft_id,
                "aero_database_version": result.config.params.aero_database_version,
                **result.state_at_index(index).as_dict(),
            })
    finally:
        if should_close:
            stream.close()
