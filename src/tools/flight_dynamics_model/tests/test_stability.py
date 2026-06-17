import json
import os
import subprocess
import sys

import numpy as np

from spearhead.config import SimulationConfig
from spearhead.stability import analyze_stability
from spearhead.stability.export import export_json, export_modes_csv
from spearhead.stability.linearize import linearize
from spearhead.stability.modes import modal_analysis
from spearhead.stability.sweep import run_cg_sweep


def test_linearization_output_dimensions():
    config = SimulationConfig(name="linearize_test", duration=0.01, n_points=2)

    result = linearize(config)

    assert result.A.shape == (12, 12)
    assert result.B.shape == (12, 4)
    assert result.A_longitudinal.shape == (4, 4)
    assert result.B_longitudinal.shape == (4, 2)
    assert result.A_lateral_directional.shape == (5, 5)
    assert result.B_lateral_directional.shape == (5, 2)


def test_modal_metrics_for_known_stable_system():
    A = np.array([[-1.0, 0.0], [0.0, -2.0]])

    result = modal_analysis(A, "synthetic", labels=("slow", "fast"))

    assert set(mode.classification for mode in result.modes) == {"stable"}
    assert all(mode.time_to_half_s is not None for mode in result.modes)
    assert all(mode.time_to_double_s is None for mode in result.modes)


def test_modal_metrics_for_known_unstable_system():
    A = np.array([[0.5, 0.0], [0.0, -1.0]])

    result = modal_analysis(A, "synthetic", labels=("unstable", "stable"))

    classifications = {mode.mode: mode.classification for mode in result.modes}
    assert classifications["unstable"] == "unstable"
    assert any(mode.time_to_double_s is not None for mode in result.modes)


def test_analyze_stability_returns_modes_and_derivatives():
    config = SimulationConfig(name="stability_test", duration=0.01, n_points=2)

    result = analyze_stability(config)

    assert len(result.longitudinal.modes) == 4
    assert len(result.lateral_directional.modes) == 5
    assert "Cm_alpha_per_rad" in result.static_derivatives.derivatives
    assert result.linearization.trim_state.shape == (12,)


def test_cg_sweep_executes_without_mutating_base_config():
    config = SimulationConfig(name="cg_sweep_test", duration=0.01, n_points=2)
    original_cg = config.params.cg_body_xyz.copy()
    cg_positions = [
        np.array([-0.15, 0.0, 0.15]),
        np.array([-0.05, 0.0, 0.15]),
    ]

    result = run_cg_sweep(config, cg_positions)

    assert len(result.points) == 2
    assert all(point.success for point in result.points)
    np.testing.assert_allclose(config.params.cg_body_xyz, original_cg)
    np.testing.assert_allclose(result.points[1].result.config.params.cg_body_xyz, cg_positions[1])


def test_stability_exports_json_and_csv(tmp_path):
    result = analyze_stability(SimulationConfig(name="export_test", duration=0.01, n_points=2))
    json_path = tmp_path / "stability.json"
    csv_path = tmp_path / "modes.csv"

    export_json(result, json_path)
    export_modes_csv(result, csv_path)

    payload = json.loads(json_path.read_text())
    assert payload["name"] == "export_test"
    assert "modes" in payload
    csv_text = csv_path.read_text()
    assert "eigenvalue_real" in csv_text
    assert "longitudinal" in csv_text


def test_stability_cli_smoke(tmp_path):
    json_path = tmp_path / "stability.json"
    env = os.environ.copy()
    env["PYTHONPATH"] = "src/tools/flight_dynamics_model"
    env["PYTHONDONTWRITEBYTECODE"] = "1"

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "spearhead.stability",
            "analyze",
            "src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml",
            "--json",
            str(json_path),
        ],
        check=True,
        env=env,
        text=True,
        capture_output=True,
    )

    assert "stability analysis complete" in completed.stdout
    assert json_path.exists()
