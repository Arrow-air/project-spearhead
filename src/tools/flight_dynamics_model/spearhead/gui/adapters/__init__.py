"""Thin GUI adapters over the shared Spearhead backend APIs."""

from .exports import (
    export_simulation_csv,
    export_simulation_json,
    export_simulation_json_text,
    export_simulation_csv_text,
    export_stability_json,
    export_stability_json_text,
    export_stability_markdown,
    export_stability_markdown_text,
    export_stability_modes_csv,
    export_stability_modes_csv_text,
)
from .scenarios import load_scenario_config, scenario_summary
from .simulation import (
    build_simulation_config,
    run_simulation_for_gui,
    simulation_result_summary,
    simulation_time_history_rows,
)
from .stability import (
    cg_sweep_summary_rows,
    run_cg_sweep_for_gui,
    run_stability_for_gui,
    stability_mode_table_rows,
    static_derivative_table_rows,
)

__all__ = [
    "cg_sweep_summary_rows",
    "export_simulation_csv",
    "export_simulation_csv_text",
    "export_simulation_json",
    "export_simulation_json_text",
    "export_stability_json",
    "export_stability_json_text",
    "export_stability_markdown",
    "export_stability_markdown_text",
    "export_stability_modes_csv",
    "export_stability_modes_csv_text",
    "load_scenario_config",
    "build_simulation_config",
    "run_cg_sweep_for_gui",
    "run_simulation_for_gui",
    "run_stability_for_gui",
    "scenario_summary",
    "simulation_result_summary",
    "simulation_time_history_rows",
    "stability_mode_table_rows",
    "static_derivative_table_rows",
]
