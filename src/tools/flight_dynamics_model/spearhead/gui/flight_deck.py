"""Flight Deck session helpers for the NiceGUI realtime console."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import asin, atan2, degrees, sqrt
from typing import Any

import numpy as np

from ..config import ControlInputCommand, SimulationConfig
from ..controls import control_at
from ..realtime import RealtimeSimulation
from ..result import SimulationResult
from ..state import STATE_COLUMNS, SimulationState


READY = "READY"
RUNNING = "RUNNING"
PAUSED = "PAUSED"
COMPLETE = "COMPLETE"
ERROR = "ERROR"

DEFAULT_HISTORY_LIMIT = 600
DEFAULT_STEP_DT = 0.02
MIN_ADVANCE_DT = 1e-9
ANGLE_SURFACES = {"elevator", "de", "aileron", "da", "rudder", "ruddervator", "dr"}


@dataclass(frozen=True)
class ManeuverDefinition:
    """Small GUI-owned preset mapped to backend control commands."""

    key: str
    label: str
    commands: tuple[ControlInputCommand, ...]
    description: str


def predefined_maneuvers(config: SimulationConfig | None = None) -> list[ManeuverDefinition]:
    """Return the Flight Deck v1 maneuver menu."""
    maneuvers: list[ManeuverDefinition] = []
    if config is not None and config.input_commands:
        maneuvers.append(
            ManeuverDefinition(
                key="scenario",
                label="Scenario inputs",
                commands=config.input_commands,
                description="Use the input schedule from the selected scenario.",
            )
        )
    maneuvers.extend(
        [
            ManeuverDefinition(
                key="pitch_doublet",
                label="Pitch doublet",
                commands=(
                    ControlInputCommand(
                        surface="elevator",
                        kind="doublet",
                        amplitude=np.deg2rad(5.0),
                        start_time=2.0,
                        duration=2.0,
                    ),
                ),
                description="Elevator doublet, 5 deg, starting at 2 s.",
            ),
            ManeuverDefinition(
                key="elevator_step",
                label="Elevator step",
                commands=(
                    ControlInputCommand(
                        surface="elevator",
                        kind="step",
                        amplitude=np.deg2rad(3.0),
                        start_time=2.0,
                        duration=2.0,
                    ),
                ),
                description="Elevator step, 3 deg, starting at 2 s.",
            ),
            ManeuverDefinition(
                key="roll_doublet",
                label="Roll doublet",
                commands=(
                    ControlInputCommand(
                        surface="aileron",
                        kind="doublet",
                        amplitude=np.deg2rad(5.0),
                        start_time=2.0,
                        duration=2.0,
                    ),
                ),
                description="Aileron doublet, 5 deg, starting at 2 s.",
            ),
            ManeuverDefinition(
                key="aileron_step",
                label="Aileron step",
                commands=(
                    ControlInputCommand(
                        surface="aileron",
                        kind="step",
                        amplitude=np.deg2rad(3.0),
                        start_time=2.0,
                        duration=2.0,
                    ),
                ),
                description="Aileron step, 3 deg, starting at 2 s.",
            ),
            ManeuverDefinition(
                key="rudder_doublet",
                label="Rudder doublet",
                commands=(
                    ControlInputCommand(
                        surface="rudder",
                        kind="doublet",
                        amplitude=np.deg2rad(5.0),
                        start_time=2.0,
                        duration=2.0,
                    ),
                ),
                description="Rudder doublet, 5 deg, starting at 2 s.",
            ),
            ManeuverDefinition(
                key="rudder_step",
                label="Rudder step",
                commands=(
                    ControlInputCommand(
                        surface="rudder",
                        kind="step",
                        amplitude=np.deg2rad(3.0),
                        start_time=2.0,
                        duration=2.0,
                    ),
                ),
                description="Rudder step, 3 deg, starting at 2 s.",
            ),
            ManeuverDefinition(
                key="throttle_step",
                label="Throttle step",
                commands=(
                    ControlInputCommand(
                        surface="throttle",
                        kind="step",
                        amplitude=0.08,
                        start_time=2.0,
                        duration=2.0,
                    ),
                ),
                description="Throttle increment, 8 percent, starting at 2 s.",
            ),
        ]
    )
    return maneuvers


def maneuver_options(config: SimulationConfig | None = None) -> dict[str, str]:
    """Return NiceGUI select options keyed by maneuver id."""
    return {maneuver.key: maneuver.label for maneuver in predefined_maneuvers(config)}


def maneuver_by_key(config: SimulationConfig, key: str) -> ManeuverDefinition:
    """Return a maneuver preset by key, falling back to the first available option."""
    maneuvers = predefined_maneuvers(config)
    for maneuver in maneuvers:
        if maneuver.key == key:
            return maneuver
    return maneuvers[0]


def config_for_maneuver(config: SimulationConfig, maneuver: ManeuverDefinition) -> SimulationConfig:
    """Build a simulation config that keeps scenario-owned settings and swaps the input schedule."""
    return SimulationConfig(
        name=f"{config.name}_{maneuver.key}",
        duration=config.duration,
        dt=config.dt,
        n_points=config.n_points,
        params=config.params,
        initial_state=config.initial_state,
        control_override=config.control_override,
        input_commands=maneuver.commands,
        start_from_trim=config.start_from_trim,
        trim_target_speed=config.trim_target_speed,
        trim_altitude=config.trim_altitude,
        rtol=config.rtol,
        atol=config.atol,
    )


def nominal_step_dt(config: SimulationConfig) -> float:
    """Return a practical realtime step size for this config."""
    if config.dt is not None:
        return float(config.dt)
    if config.n_points is not None:
        return float(config.duration / (config.n_points - 1))
    return DEFAULT_STEP_DT


def state_telemetry(state: SimulationState) -> dict[str, float]:
    """Return display telemetry derived from the existing 12-state vector."""
    row = state.as_dict()
    u = row["u"]
    v = row["v"]
    w = row["w"]
    airspeed = sqrt(u * u + v * v + w * w)
    beta = degrees(asin(max(-1.0, min(1.0, v / airspeed)))) if airspeed > 1e-9 else 0.0
    return {
        "time": row["time"],
        "airspeed": airspeed,
        "altitude": -row["pd"],
        "pitch_deg": degrees(row["theta"]),
        "roll_deg": degrees(row["phi"]),
        "heading_deg": degrees(row["psi"]),
        "alpha_deg": degrees(atan2(w, u)) if abs(u) > 1e-9 or abs(w) > 1e-9 else 0.0,
        "beta_deg": beta,
        "p_deg_s": degrees(row["p"]),
        "q_deg_s": degrees(row["q"]),
        "r_deg_s": degrees(row["r"]),
    }


def control_telemetry(simulation: RealtimeSimulation | None, t: float = 0.0) -> dict[str, float]:
    """Return current control values from the backend control schedule."""
    control = control_at(t, None if simulation is None else simulation.control_override)
    return {
        "elevator_deg": degrees(control.de),
        "aileron_deg": degrees(control.da),
        "rudder_deg": degrees(control.dr),
        "throttle_pct": 100.0 * control.throttle,
    }


@dataclass
class FlightDeckSession:
    """Frontend state controller for realtime Flight Deck playback."""

    scenario_source: str
    base_config: SimulationConfig
    maneuver_key: str = ""
    speed_multiplier: float = 1.0
    history_limit: int = DEFAULT_HISTORY_LIMIT
    status: str = READY
    message: str = ""
    simulation: RealtimeSimulation | None = None
    history: list[dict[str, float]] = field(default_factory=list)
    run_history: list[dict[str, float]] = field(default_factory=list)
    latest_completed_result: SimulationResult | None = None
    _initial_state: np.ndarray | None = None

    def __post_init__(self) -> None:
        if not self.maneuver_key:
            self.maneuver_key = predefined_maneuvers(self.base_config)[0].key
        self.reset()

    @property
    def maneuver(self) -> ManeuverDefinition:
        """Return the currently selected maneuver definition."""
        return maneuver_by_key(self.base_config, self.maneuver_key)

    @property
    def config(self) -> SimulationConfig:
        """Return the active scenario config with the selected maneuver schedule."""
        return config_for_maneuver(self.base_config, self.maneuver)

    @property
    def current_state(self) -> SimulationState | None:
        """Return the latest realtime state."""
        return None if self.simulation is None else self.simulation.state

    def set_scenario(self, source: str, config: SimulationConfig) -> None:
        """Switch scenario and reset playback."""
        self.scenario_source = source
        self.base_config = config
        if self.maneuver_key not in maneuver_options(config):
            self.maneuver_key = predefined_maneuvers(config)[0].key
        self.latest_completed_result = None
        self.reset()

    def set_maneuver(self, key: str) -> None:
        """Switch maneuver and reset playback."""
        self.maneuver_key = key
        self.latest_completed_result = None
        self.reset()

    def reset(self) -> SimulationState:
        """Create a fresh realtime simulation at t=0."""
        self.simulation = RealtimeSimulation(self.config)
        self.status = READY
        self.message = "Ready"
        state = self.simulation.state
        self._initial_state = state.vector.copy()
        self.history = []
        self.run_history = []
        self._append_history(state)
        return state

    def play(self) -> None:
        """Start or resume realtime playback."""
        if self.status == COMPLETE:
            self.reset()
        if self.simulation is None:
            self.reset()
        self.status = RUNNING
        self.message = "Running"

    def pause(self) -> None:
        """Pause playback without discarding state."""
        if self.status == RUNNING:
            self.status = PAUSED
            self.message = "Paused"

    def toggle_playback(self) -> None:
        """Toggle between running and paused/ready states."""
        if self.status == RUNNING:
            self.pause()
        else:
            self.play()

    def step_once(self) -> SimulationState | None:
        """Advance one simulation step while staying paused."""
        if self.simulation is None or self.status == COMPLETE:
            self.reset()
        state = self._step(nominal_step_dt(self.config))
        if self.status != COMPLETE:
            self.status = PAUSED
            self.message = "Paused"
        return state

    def advance(self, wall_dt: float) -> SimulationState | None:
        """Advance simulation by wall-clock delta scaled by speed multiplier."""
        if self.status != RUNNING or self.simulation is None:
            return self.current_state
        target_dt = max(0.0, wall_dt) * self.speed_multiplier
        step_dt = nominal_step_dt(self.config)
        elapsed = 0.0
        state = self.current_state
        while target_dt - elapsed > MIN_ADVANCE_DT and self.status == RUNNING:
            actual_step = min(step_dt, target_dt - elapsed)
            if actual_step <= MIN_ADVANCE_DT:
                break
            state = self._step(actual_step)
            elapsed += actual_step
        return state

    def telemetry(self) -> dict[str, float]:
        """Return display telemetry for the current state."""
        state = self.current_state
        return {} if state is None else state_telemetry(state)

    def controls(self) -> dict[str, float]:
        """Return display control positions for the current time."""
        t = 0.0 if self.current_state is None else self.current_state.time
        return control_telemetry(self.simulation, t)

    def latest_result(self) -> SimulationResult | None:
        """Return a completed result, or a partial result if history exists."""
        if self.latest_completed_result is not None:
            return self.latest_completed_result
        if not self.run_history or self._initial_state is None:
            return None
        return self._result_from_history(success=self.status != ERROR, message=self.message)

    def _step(self, dt: float) -> SimulationState | None:
        if self.simulation is None:
            return None
        try:
            state = self.simulation.step(dt)
            self._append_history(state)
            if not self.simulation.running:
                self.status = COMPLETE
                self.message = "Complete"
                self.latest_completed_result = self._result_from_history(success=True, message="Realtime run complete")
            return state
        except Exception as exc:  # noqa: BLE001 - GUI should surface backend errors.
            self.status = ERROR
            self.message = str(exc)
            return self.current_state

    def _append_history(self, state: SimulationState) -> None:
        row = state.as_dict()
        row.update(self.controls())
        self.run_history.append(row)
        self.history.append(row)
        if len(self.history) > self.history_limit:
            del self.history[: len(self.history) - self.history_limit]

    def _result_from_history(self, *, success: bool, message: str) -> SimulationResult:
        time = np.asarray([row["time"] for row in self.run_history], dtype=float)
        states = np.asarray([[row[column] for column in STATE_COLUMNS] for row in self.run_history], dtype=float)
        initial_state = self._initial_state if self._initial_state is not None else states[0]
        return SimulationResult(
            config=self.config,
            time=time,
            states=states,
            success=success,
            message=message,
            initial_state=initial_state,
            control_override=None if self.simulation is None else self.simulation.control_override,
            trim_result=None if self.simulation is None else self.simulation.trim_result,
        )


def _pitch_ladder_html() -> str:
    marks = []
    pixels_per_degree = 7.0
    for pitch in (30, 20, 10, 0, -10, -20, -30):
        y = -pitch * pixels_per_degree
        width = 148 if pitch == 0 else 88
        label = f"{pitch:+d}" if pitch != 0 else "0"
        marks.append(
            (
                f'<div class="fd-pfd__pitch-mark" style="top: calc(50% + {y:.1f}px);">'
                f'<span class="fd-pfd__pitch-label fd-pfd__pitch-label--left">{label}</span>'
                f'<i style="width: {width}px;"></i>'
                f'<span class="fd-pfd__pitch-label fd-pfd__pitch-label--right">{label}</span>'
                "</div>"
            )
        )
    return "".join(marks)


def _roll_scale_html() -> str:
    marks = []
    for angle in (-60, -45, -30, -20, -10, 10, 20, 30, 45, 60):
        length_class = "fd-pfd__roll-tick--major" if abs(angle) in {30, 45, 60} else ""
        label = str(abs(angle)) if abs(angle) in {30, 60} else ""
        marks.append(
            (
                f'<div class="fd-pfd__roll-tick {length_class}" '
                f'style="transform: translateX(-50%) rotate({angle:.1f}deg);">'
                f"<span>{label}</span>"
                "</div>"
            )
        )
    return "".join(marks)


def attitude_html(*, pitch_deg: float, roll_deg: float) -> str:
    """Return an aviation-style artificial horizon / PFD attitude display."""
    pitch_offset = max(-220.0, min(220.0, pitch_deg * 7.0))
    roll = max(-75.0, min(75.0, roll_deg))
    ladder = _pitch_ladder_html()
    roll_scale = _roll_scale_html()
    return f"""
    <div class="fd-pfd" aria-label="Primary attitude display">
      <div class="fd-pfd__bezel"></div>
      <div class="fd-pfd__world" style="transform: rotate({-roll:.3f}deg);">
        <div class="fd-pfd__pitch" style="transform: translateY({pitch_offset:.3f}px);">
          <div class="fd-pfd__sky"></div>
          <div class="fd-pfd__ground"></div>
          <div class="fd-pfd__horizon"></div>
          <div class="fd-pfd__ladder">{ladder}</div>
        </div>
      </div>
      <div class="fd-pfd__roll-scale">
        {roll_scale}
        <div class="fd-pfd__roll-zero"></div>
      </div>
      <div class="fd-pfd__roll-pointer" style="transform: translateX(-50%) rotate({-roll:.3f}deg);"></div>
      <div class="fd-pfd__aircraft">
        <span class="fd-pfd__wing fd-pfd__wing--left"></span>
        <strong></strong>
        <span class="fd-pfd__wing fd-pfd__wing--right"></span>
        <em></em>
      </div>
      <div class="fd-pfd__readout">
        <div><span>ROLL</span><strong>{roll_deg:+.1f}</strong><small>deg</small></div>
        <div><span>PITCH</span><strong>{pitch_deg:+.1f}</strong><small>deg</small></div>
      </div>
    </div>
    """
