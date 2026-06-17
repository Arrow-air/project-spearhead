import subprocess
import sys

import pytest

from spearhead.config import SimulationConfig
from spearhead.gui.adapters.simulation import run_simulation_for_gui
from spearhead.gui.adapters.stability import build_linearization_config, run_stability_for_gui


def test_gui_app_helpers_do_not_require_nicegui():
    from spearhead.gui.app import build_parser, get_example_scenarios

    scenarios = get_example_scenarios()
    args = build_parser().parse_args(["--host", "0.0.0.0", "--port", "9090", "--no-show"])

    assert any(path.name == "pitch_doublet.yaml" for path in scenarios)
    assert args.host == "0.0.0.0"
    assert args.port == 9090
    assert args.no_show is True


def test_gui_app_registers_with_nicegui_when_available():
    pytest.importorskip("nicegui")

    from spearhead.gui.app import create_app

    create_app()


def test_gui_module_help_smoke():
    completed = subprocess.run(
        [sys.executable, "-m", "spearhead.gui", "--help"],
        check=True,
        text=True,
        capture_output=True,
    )

    assert "python -m spearhead.gui" in completed.stdout


def test_simulation_form_helper_uses_adapter_constructor(monkeypatch):
    from spearhead.gui import app

    calls = []

    def fake_build_simulation_config(**overrides):
        calls.append(overrides)
        return SimulationConfig(**overrides)

    monkeypatch.setattr(app, "build_simulation_config", fake_build_simulation_config)
    base_config = SimulationConfig(name="base", duration=1.0, n_points=3)

    config = app.build_simulation_config_from_form(
        base_config,
        name="edited",
        duration=2.0,
        dt=0.1,
        n_points=None,
        start_from_trim=False,
        trim_target_speed=24.0,
        trim_altitude=90.0,
    )

    assert calls
    assert config.name == "edited"
    assert config.duration == 2.0
    assert config.dt == 0.1
    assert config.n_points is None
    assert config.params is base_config.params
    assert config.input_commands == base_config.input_commands


def test_stability_summary_helper_formats_backend_result():
    from spearhead.gui.app import stability_result_summary_rows

    result = run_stability_for_gui(
        SimulationConfig(name="stability_summary_test", duration=0.01, n_points=2),
        linearization_config=build_linearization_config(
            state_perturbation=1e-5,
            control_perturbation=1e-5,
        ),
    )
    rows = stability_result_summary_rows(result)
    values = {row["field"]: row["value"] for row in rows}

    assert values["success"] == "True"
    assert values["A_shape"] == "(12, 12)"
    assert values["B_shape"] == "(12, 4)"
    assert values["A_longitudinal_shape"] == "(4, 4)"
    assert "theta" in values["longitudinal_state_labels"]


def test_simulation_plot_helper_returns_serializable_plotly_figure():
    pytest.importorskip("plotly")

    from spearhead.gui.plotting import figure_to_json_dict, simulation_time_history_figure

    result = run_simulation_for_gui(SimulationConfig(name="plot_test", duration=0.1, n_points=3))
    figure = simulation_time_history_figure(result)
    payload = figure_to_json_dict(figure)

    assert "data" in payload
    assert len(payload["data"]) == 10
    assert payload["data"][0]["name"] == "u"


def test_stability_plot_helper_returns_serializable_plotly_figure():
    pytest.importorskip("plotly")

    from spearhead.gui.plotting import figure_to_json_dict, stability_eigenvalue_figure

    result = run_stability_for_gui(SimulationConfig(name="pole_plot_test", duration=0.01, n_points=2))
    figure = stability_eigenvalue_figure(result)
    payload = figure_to_json_dict(figure)

    assert "data" in payload
    assert len(payload["data"]) == 2
    assert {trace["name"] for trace in payload["data"]} == {"longitudinal", "lateral_directional"}
