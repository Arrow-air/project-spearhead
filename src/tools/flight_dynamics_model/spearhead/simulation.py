"""Simulation entry points."""

from collections.abc import Iterable

import numpy as np
from scipy.integrate import solve_ivp

from .config import ControlInputCommand, SimulationConfig
from .controls import Control, ControlInput, control_at
from .dynamics import dynamics
from .params import AircraftParams
from .result import SimulationResult
from .trim import trim_fixed_wing_longitudinal


def default_initial_state() -> np.ndarray:
    """Return the default 12-state initial condition."""
    return np.array(
        [
            0.0,
            0.0,
            -100.0,
            22.0,
            0.0,
            0.0,
            np.deg2rad(0.0),
            np.deg2rad(2.0),
            np.deg2rad(0.0),
            0.0,
            0.0,
            0.0,
        ]
    )


def _time_eval(config: SimulationConfig) -> np.ndarray:
    if config.n_points is not None:
        return np.linspace(0.0, config.duration, config.n_points)

    if config.dt is None:
        raise ValueError("SimulationConfig requires either n_points or dt")

    t_eval = np.arange(0.0, config.duration + 0.5 * config.dt, config.dt)
    if t_eval[-1] < config.duration:
        t_eval = np.append(t_eval, config.duration)
    if t_eval[-1] > config.duration:
        t_eval[-1] = config.duration
    return t_eval


def _surface_to_channel(surface: str) -> str:
    normalized = surface.strip().lower()
    mapping = {
        "elevator": "de",
        "de": "de",
        "aileron": "da",
        "da": "da",
        "rudder": "dr",
        "ruddervator": "dr",
        "dr": "dr",
        "throttle": "throttle",
    }
    if normalized not in mapping:
        raise ValueError(f"Unsupported control surface: {surface}")
    return mapping[normalized]


def _control_with_commands(
    base_control: ControlInput | None,
    commands: Iterable[ControlInputCommand],
) -> ControlInput | None:
    commands = tuple(commands)
    if not commands:
        return base_control

    channels = tuple((_surface_to_channel(command.surface), command) for command in commands)

    def scheduled_control(t: float) -> Control:
        base = control_at(t, base_control)
        values = {
            "de": base.de,
            "da": base.da,
            "dr": base.dr,
            "throttle": base.throttle,
        }
        for channel, command in channels:
            values[channel] += command.value_at(t)
        return Control(**values)

    return scheduled_control


def _prepare_initial_state_and_control(
    config: SimulationConfig,
) -> tuple[np.ndarray, ControlInput | None, dict | None]:
    trim_result = None
    initial_state = config.initial_state
    control_override = config.control_override

    if initial_state is None and config.start_from_trim:
        trim_result = trim_fixed_wing_longitudinal(
            config.params,
            target_speed=config.trim_target_speed,
            altitude=config.trim_altitude,
        )
        initial_state = trim_result["x_trim"]
        if control_override is None:
            control_override = trim_result["control_trim"]

    y0 = default_initial_state() if initial_state is None else np.asarray(initial_state, dtype=float)
    control_override = _control_with_commands(control_override, config.input_commands)
    return y0, control_override, trim_result


def run_simulation(config: SimulationConfig) -> SimulationResult:
    """Run the shared nonlinear simulation backend."""
    y0, control_override, trim_result = _prepare_initial_state_and_control(config)
    t_eval = _time_eval(config)
    sol = solve_ivp(
        fun=lambda t, x: dynamics(t, x, config.params, control_override=control_override),
        t_span=(0.0, config.duration),
        y0=y0,
        t_eval=t_eval,
        rtol=config.rtol,
        atol=config.atol,
    )
    sol.initial_state = y0
    sol.control_override = control_override
    sol.trim_result = trim_result
    return SimulationResult.from_solve_ivp(
        config=config,
        solution=sol,
        initial_state=y0,
        control_override=control_override,
        trim_result=trim_result,
    )


def run_open_loop(
    t_final: float = 20.0,
    n_points: int = 1001,
    params: AircraftParams | None = None,
    initial_state: np.ndarray | None = None,
    control_override: Control | dict[str, float] | None = None,
    start_from_trim: bool = True,
    trim_target_speed: float = 22.0,
    trim_altitude: float = 100.0,
):
    """Run an open-loop nonlinear simulation.

    By default this starts from the preliminary longitudinal trim condition and
    holds the trimmed control input constant. Pass ``start_from_trim=False`` to
    use the original arbitrary default initial state and scheduled controls.
    """
    config = SimulationConfig(
        name="open_loop",
        duration=t_final,
        n_points=n_points,
        params=AircraftParams() if params is None else params,
        initial_state=initial_state,
        control_override=control_override,
        start_from_trim=start_from_trim,
        trim_target_speed=trim_target_speed,
        trim_altitude=trim_altitude,
    )
    return run_simulation(config).raw_solution
