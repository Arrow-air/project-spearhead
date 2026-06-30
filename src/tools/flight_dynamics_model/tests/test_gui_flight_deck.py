import numpy as np
import pytest

from spearhead.config import ControlInputCommand, SimulationConfig
from spearhead.gui.flight_deck import (
    COMPLETE,
    PAUSED,
    READY,
    RUNNING,
    FlightDeckSession,
    attitude_html,
    config_for_maneuver,
    maneuver_by_key,
    maneuver_options,
    predefined_maneuvers,
    state_telemetry,
)
from spearhead.state import SimulationState


def small_config(**overrides):
    values = {
        "name": "flight_deck_test",
        "duration": 0.06,
        "dt": 0.02,
        "n_points": None,
        "start_from_trim": False,
    }
    values.update(overrides)
    return SimulationConfig(**values)


def test_maneuver_options_include_scenario_inputs_when_present():
    config = small_config(
        input_commands=(
            ControlInputCommand(
                surface="elevator",
                kind="doublet",
                amplitude=np.deg2rad(2.0),
                start_time=1.0,
                duration=1.0,
            ),
        )
    )

    options = maneuver_options(config)

    assert options["scenario"] == "Scenario inputs"
    assert options["pitch_doublet"] == "Pitch doublet"
    assert options["elevator_step"] == "Elevator step"
    assert options["roll_doublet"] == "Roll doublet"
    assert options["aileron_step"] == "Aileron step"
    assert options["rudder_doublet"] == "Rudder doublet"
    assert options["rudder_step"] == "Rudder step"
    assert options["throttle_step"] == "Throttle step"


def test_config_for_maneuver_preserves_scenario_owned_fields():
    config = small_config(trim_target_speed=25.0, trim_altitude=120.0)
    maneuver = maneuver_by_key(config, "elevator_step")

    maneuver_config = config_for_maneuver(config, maneuver)

    assert maneuver_config.name == "flight_deck_test_elevator_step"
    assert maneuver_config.params is config.params
    assert maneuver_config.duration == config.duration
    assert maneuver_config.dt == config.dt
    assert maneuver_config.trim_target_speed == 25.0
    assert maneuver_config.input_commands == maneuver.commands


def test_flight_deck_session_transport_transitions():
    session = FlightDeckSession("unit", small_config())

    assert session.status == READY
    assert len(session.history) == 1

    session.play()
    assert session.status == RUNNING

    session.advance(0.02)
    assert session.current_state.time > 0.0

    session.pause()
    assert session.status == PAUSED

    session.step_once()
    assert session.status in {PAUSED, COMPLETE}
    assert len(session.history) >= 3

    session.reset()
    assert session.status == READY
    assert session.current_state.time == 0.0
    assert len(session.history) == 1


def test_flight_deck_session_records_completed_result():
    session = FlightDeckSession("unit", small_config(duration=0.04))

    session.play()
    session.advance(1.0)

    assert session.status == COMPLETE
    result = session.latest_completed_result
    assert result is not None
    assert result.success is True
    assert result.time[-1] == 0.04
    assert result.states.shape[1] == 12


def test_flight_deck_session_avoids_tiny_remainder_steps():
    session = FlightDeckSession("unit", small_config(duration=1.0, dt=0.01))

    session.play()
    session.advance(0.1)
    session.advance(0.1)

    assert session.status == RUNNING
    assert session.message == "Running"
    assert session.current_state.time == pytest.approx(0.2)


def test_flight_deck_session_bounds_rolling_history():
    session = FlightDeckSession("unit", small_config(duration=0.2), history_limit=3)

    session.play()
    session.advance(0.12)

    assert len(session.history) == 3
    assert len(session.run_history) > len(session.history)
    assert session.history[-1]["time"] > session.history[0]["time"]


def test_completed_result_uses_full_run_history_not_rolling_display_history():
    session = FlightDeckSession("unit", small_config(duration=0.08), history_limit=2)

    session.play()
    session.advance(1.0)

    assert session.status == COMPLETE
    assert len(session.history) == 2
    assert len(session.run_history) > len(session.history)
    assert session.latest_completed_result is not None
    assert len(session.latest_completed_result.time) == len(session.run_history)


def test_flight_deck_control_and_state_telemetry_are_display_ready():
    command = ControlInputCommand(
        surface="elevator",
        kind="step",
        amplitude=np.deg2rad(5.0),
        start_time=0.0,
        duration=1.0,
    )
    session = FlightDeckSession(
        "unit",
        small_config(input_commands=(command,)),
        maneuver_key="scenario",
    )

    telemetry = session.telemetry()
    controls = session.controls()

    assert telemetry["airspeed"] > 0.0
    assert telemetry["altitude"] == 100.0
    assert controls["elevator_deg"] == 5.0
    assert "throttle_pct" in controls
    assert {"elevator_deg", "aileron_deg", "rudder_deg", "throttle_pct"}.issubset(session.history[0])


def test_lateral_directional_maneuver_commands_use_existing_command_model():
    config = small_config()
    expected = {
        "roll_doublet": ("aileron", "doublet", 5.0),
        "aileron_step": ("aileron", "step", 3.0),
        "rudder_doublet": ("rudder", "doublet", 5.0),
        "rudder_step": ("rudder", "step", 3.0),
    }

    for key, (surface, kind, amplitude_deg) in expected.items():
        maneuver = maneuver_by_key(config, key)
        command = maneuver.commands[0]
        assert command.surface == surface
        assert command.kind == kind
        assert np.rad2deg(command.amplitude) == pytest.approx(amplitude_deg)


def test_maneuver_command_values_over_time():
    config = small_config()

    doublet = maneuver_by_key(config, "roll_doublet").commands[0]
    assert doublet.value_at(1.9) == 0.0
    assert np.rad2deg(doublet.value_at(2.1)) == pytest.approx(5.0)
    assert np.rad2deg(doublet.value_at(3.1)) == pytest.approx(-5.0)
    assert doublet.value_at(4.1) == 0.0

    throttle = maneuver_by_key(config, "throttle_step").commands[0]
    assert throttle.value_at(1.9) == 0.0
    assert throttle.value_at(2.0) == pytest.approx(0.08)


def test_state_telemetry_computes_alpha_beta_and_angles():
    state = SimulationState(
        1.5,
        np.array(
            [
                0.0,
                0.0,
                -50.0,
                10.0,
                1.0,
                1.0,
                np.deg2rad(10.0),
                np.deg2rad(5.0),
                np.deg2rad(30.0),
                0.0,
                0.0,
                0.0,
            ]
        ),
    )

    telemetry = state_telemetry(state)

    assert telemetry["time"] == 1.5
    assert telemetry["altitude"] == 50.0
    assert telemetry["roll_deg"] == pytest.approx(10.0)
    assert telemetry["pitch_deg"] == pytest.approx(5.0)
    assert telemetry["heading_deg"] == pytest.approx(30.0)
    assert telemetry["alpha_deg"] > 0.0
    assert telemetry["beta_deg"] > 0.0


def test_attitude_html_contains_pitch_and_roll_readout():
    html = attitude_html(pitch_deg=3.25, roll_deg=-4.5)

    assert "<span>PITCH</span><strong>+3.2</strong><small>deg</small>" in html
    assert "<span>ROLL</span><strong>-4.5</strong><small>deg</small>" in html
    assert "fd-pfd" in html
    assert "fd-pfd__aircraft" in html


def test_attitude_html_contains_pitch_ladder_and_roll_scale():
    html = attitude_html(pitch_deg=12.0, roll_deg=35.0)

    for label in ("+30", "+20", "+10", "0", "-10", "-20", "-30"):
        assert label in html
    assert "fd-pfd__roll-scale" in html
    assert "fd-pfd__roll-pointer" in html
    assert "rotate(-35.000deg)" in html


def test_predefined_maneuvers_work_without_scenario_commands():
    maneuvers = predefined_maneuvers(small_config())

    assert [maneuver.key for maneuver in maneuvers] == [
        "pitch_doublet",
        "elevator_step",
        "roll_doublet",
        "aileron_step",
        "rudder_doublet",
        "rudder_step",
        "throttle_step",
    ]
