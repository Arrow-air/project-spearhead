"""Small shared state container for GUI pages."""

from __future__ import annotations

from dataclasses import dataclass

from ..config import SimulationConfig
from ..result import SimulationResult


@dataclass
class GuiState:
    """State shared by lightweight GUI pages in one local app process."""

    scenario_source: str | None = None
    config: SimulationConfig | None = None
    simulation_result: SimulationResult | None = None
    simulation_error: str | None = None

    def set_scenario(self, source: str, config: SimulationConfig) -> None:
        """Store the active scenario and clear dependent results."""
        self.scenario_source = source
        self.config = config
        self.simulation_result = None
        self.simulation_error = None

    def set_simulation_result(self, result: SimulationResult) -> None:
        """Store the latest simulation result."""
        self.simulation_result = result
        self.simulation_error = None

    def set_simulation_error(self, message: str) -> None:
        """Store the latest simulation error."""
        self.simulation_result = None
        self.simulation_error = message


APP_STATE = GuiState()

