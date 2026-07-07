import json

import numpy as np

from spearhead.config import SimulationConfig
from spearhead.gui.adapters import exports
from spearhead.gui.adapters.scenarios import load_scenario_config, scenario_summary
from spearhead.gui.adapters.simulation import (
    build_simulation_config,
    run_simulation_for_gui,
    simulation_result_summary,
    simulation_time_history_rows,
)
from spearhead.gui.adapters.stability import (
    cg_sweep_summary_rows,
    run_cg_sweep_for_gui,
    run_stability_for_gui,
    stability_mode_table_rows,
    static_derivative_table_rows,
)


def test_scenario_adapter_loads_valid_yaml_scenario():
    config = load_scenario_config("src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml")
    summary = scenario_summary(config)

    assert config.name == "pitch_doublet_nominal"
    assert summary["duration"] == 20.0
    assert summary["dt"] == 0.01
    assert summary["n_points"] is None
    assert summary["input_commands"][0]["surface"] == "elevator"


def test_simulation_adapter_returns_summary_and_table_rows():
    config = build_simulation_config(name="gui_sim", duration=0.1, n_points=3)
    result = run_simulation_for_gui(config)
    summary = simulation_result_summary(result)
    rows = simulation_time_history_rows(result)

    assert result.success
    assert summary["name"] == "gui_sim"
    assert summary["samples"] == 3
    assert summary["trim"]["success"] is True
    assert len(rows) == 3
    assert {"time", "pn", "theta", "r"}.issubset(rows[0])


def test_stability_adapter_returns_mode_and_static_derivative_rows():
    result = run_stability_for_gui(SimulationConfig(name="gui_stability", duration=0.01, n_points=2))
    mode_rows = stability_mode_table_rows(result)
    derivative_rows = static_derivative_table_rows(result)

    assert len(mode_rows) == 9
    assert {"subsystem", "mode", "eigenvalue_label", "eigenvalue_real", "classification"}.issubset(mode_rows[0])
    assert any(row["name"] == "Cm_alpha_per_rad" for row in derivative_rows)
    assert all("source" in row for row in derivative_rows)


def test_cg_sweep_adapter_returns_per_point_summary_rows():
    config = SimulationConfig(name="gui_cg_sweep", duration=0.01, n_points=2)
    sweep = run_cg_sweep_for_gui(
        config,
        (
            np.array([-0.15, 0.0, 0.15]),
            np.array([-0.05, 0.0, 0.15]),
        ),
    )
    rows = cg_sweep_summary_rows(sweep)

    assert len(rows) == 2
    assert rows[0]["case"] == "cg_0"
    assert rows[0]["success"] is True
    assert rows[0]["longitudinal_modes"] == 4
    assert rows[0]["lateral_directional_modes"] == 5


def test_export_adapter_text_helpers_use_existing_export_shapes():
    sim_result = run_simulation_for_gui(SimulationConfig(name="gui_export", duration=0.1, n_points=3))
    stability_result = run_stability_for_gui(SimulationConfig(name="gui_export_stability", duration=0.01, n_points=2))

    sim_csv = exports.export_simulation_csv_text(sim_result)
    sim_json = json.loads(exports.export_simulation_json_text(sim_result))
    stability_csv = exports.export_stability_modes_csv_text(stability_result)
    stability_json = json.loads(exports.export_stability_json_text(stability_result))
    stability_markdown = exports.export_stability_markdown_text(stability_result)

    assert sim_csv.splitlines()[0].startswith("time,pn,pe,pd")
    assert sim_json["name"] == "gui_export"
    assert "eigenvalue_real" in stability_csv
    assert stability_json["name"] == "gui_export_stability"
    assert stability_markdown.startswith("# Stability Analysis: gui_export_stability")


def test_export_adapter_delegates_to_existing_exporters(monkeypatch):
    calls = []

    def fake_export(result, destination):
        calls.append((result, destination))
        destination.write("delegated")

    monkeypatch.setattr(exports, "_export_simulation_csv", fake_export)

    result = run_simulation_for_gui(SimulationConfig(duration=0.1, n_points=3))
    text = exports.export_simulation_csv_text(result)

    assert text == "delegated"
    assert calls[0][0] is result
