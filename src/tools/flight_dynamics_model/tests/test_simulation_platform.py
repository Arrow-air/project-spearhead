import csv
import json

import numpy as np

from spearhead.config import SimulationConfig
from spearhead.export import export_csv, export_json
from spearhead.realtime import RealtimeSimulation
from spearhead.scenarios import load_scenario
from spearhead.simulation import run_open_loop, run_simulation


def test_simulation_config_accepts_valid_initial_state():
    initial_state = np.zeros(12)
    config = SimulationConfig(duration=1.0, n_points=3, initial_state=initial_state)

    assert config.duration == 1.0
    np.testing.assert_allclose(config.initial_state, initial_state)


def test_load_scenario_converts_yaml_to_config(tmp_path):
    scenario_path = tmp_path / "scenario.yaml"
    scenario_path.write_text(
        "\n".join(
            [
                "name: test_pitch_doublet",
                "aircraft_case: nominal",
                "trim:",
                "  airspeed: 25",
                "  altitude: 120",
                "simulation:",
                "  duration: 1.0",
                "  dt: 0.1",
                "inputs:",
                "  elevator:",
                "    type: doublet",
                "    amplitude_deg: 5",
                "    start_time: 0.2",
                "    duration: 0.4",
            ]
        )
    )

    config = load_scenario(scenario_path)

    assert config.name == "test_pitch_doublet"
    assert config.duration == 1.0
    assert config.dt == 0.1
    assert config.n_points is None
    assert config.trim_target_speed == 25.0
    assert config.trim_altitude == 120.0
    assert len(config.input_commands) == 1
    np.testing.assert_allclose(config.input_commands[0].amplitude, np.deg2rad(5.0))


def test_run_simulation_returns_structured_result():
    config = SimulationConfig(duration=0.1, n_points=3)
    result = run_simulation(config)

    assert result.success
    assert result.time.shape == (3,)
    assert result.states.shape == (3, 12)
    assert result.final_state.time == result.time[-1]
    assert result.trim_result is not None


def test_run_open_loop_remains_compatible_with_existing_solution_shape():
    sol = run_open_loop(t_final=0.1, n_points=3)

    assert sol.success
    assert sol.y.shape == (12, 3)
    assert hasattr(sol, "initial_state")
    assert hasattr(sol, "control_override")
    assert hasattr(sol, "trim_result")


def test_csv_export_writes_time_history(tmp_path):
    result = run_simulation(SimulationConfig(duration=0.1, n_points=3))
    csv_path = tmp_path / "result.csv"

    export_csv(result, csv_path)

    with csv_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 3
    assert rows[0]["time"] == "0.0"
    assert "theta" in rows[0]


def test_json_export_writes_result_payload(tmp_path):
    result = run_simulation(SimulationConfig(name="json_case", duration=0.1, n_points=3))
    json_path = tmp_path / "result.json"

    export_json(result, json_path)

    payload = json.loads(json_path.read_text())
    assert payload["name"] == "json_case"
    assert payload["success"] is True
    assert len(payload["states"]) == 3
    assert payload["trim"]["success"] is True


def test_realtime_simulation_step_advances_state():
    sim = RealtimeSimulation(SimulationConfig(duration=0.2, dt=0.1, n_points=None))
    initial = sim.state

    next_state = sim.step()

    assert sim.running
    assert next_state.time == 0.1
    assert next_state.vector.shape == (12,)
    assert not np.array_equal(next_state.vector, initial.vector)
