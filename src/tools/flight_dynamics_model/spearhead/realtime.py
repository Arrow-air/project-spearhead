"""Reusable real-time stepping interface for the nonlinear FDM."""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .config import SimulationConfig
from .dynamics import dynamics
from .simulation import _prepare_initial_state_and_control
from .state import SimulationState


class RealtimeSimulation:
    """Incremental simulation wrapper suitable for future GUI event loops."""

    def __init__(self, config: SimulationConfig):
        self.config = config
        self._state, self.control_override, self.trim_result = _prepare_initial_state_and_control(config)
        self.time = 0.0
        self.running = True

    @property
    def state(self) -> SimulationState:
        """Return the current timestamped simulation state."""
        return SimulationState(self.time, self._state)

    def step(self, dt: float | None = None) -> SimulationState:
        """Advance the simulation by ``dt`` seconds using the shared dynamics."""
        if not self.running:
            return self.state

        step_dt = self.config.dt if dt is None else dt
        if step_dt is None:
            if self.config.n_points is None:
                raise ValueError("RealtimeSimulation.step requires dt when config.dt is not set")
            step_dt = self.config.duration / (self.config.n_points - 1)
        if step_dt <= 0.0:
            raise ValueError("RealtimeSimulation step dt must be positive")

        next_time = min(self.time + step_dt, self.config.duration)
        sol = solve_ivp(
            fun=lambda t, x: dynamics(t, x, self.config.params, control_override=self.control_override),
            t_span=(self.time, next_time),
            y0=self._state,
            t_eval=np.array([next_time]),
            rtol=self.config.rtol,
            atol=self.config.atol,
        )
        if not sol.success:
            self.running = False
            raise RuntimeError(sol.message)

        self.time = float(next_time)
        self._state = sol.y[:, -1]
        if self.time >= self.config.duration:
            self.running = False
        return self.state
