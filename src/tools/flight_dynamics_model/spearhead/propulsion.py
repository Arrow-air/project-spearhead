"""Shared pusher and lift propulsion models with static manufacturer curves."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from .controls import Control
from .params import AircraftParams

CURVE_FIELDS = ("thrust_n", "current_a", "power_w", "rpm", "torque_nm")


@dataclass(frozen=True)
class StaticPerformance:
    voltages: np.ndarray
    throttles: np.ndarray
    values: np.ndarray

    def lookup(self, voltage_v: float, throttle: float) -> dict[str, float]:
        """Interpolate voltage/throttle within the measured voltage envelope.

        Below the first measured throttle, use a provisional linear ramp from
        motor-off. This is not an ESC startup, idle, or transient model.
        """
        if not np.isfinite(voltage_v) or not self.voltages[0] <= voltage_v <= self.voltages[-1]:
            raise ValueError(f"Lift voltage {voltage_v} V is outside the curve range {self.voltages[0]}–{self.voltages[-1]} V")
        if not np.isfinite(throttle):
            raise ValueError("Lift throttle must be finite")
        throttle = float(np.clip(throttle, 0.0, 1.0))
        samples = np.array([
            [np.interp(throttle, self.throttles, values[:, column]) for column in range(len(CURVE_FIELDS))]
            for values in self.values
        ])
        return {
            field: float(np.interp(voltage_v, self.voltages, samples[:, column]))
            for column, field in enumerate(CURVE_FIELDS)
        }


@lru_cache(maxsize=8)
def _load_performance(path: str, modified_ns: int, size: int) -> StaticPerformance:
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    voltages = np.array(sorted({float(row["voltage_v"]) for row in rows}))
    throttles = np.array([0.0] + sorted({float(row["throttle"]) for row in rows}))
    if voltages.size < 2 or len(throttles) < 3 or throttles[-1] != 1.0 or np.any(np.diff(throttles) <= 0):
        raise ValueError("Lift performance needs an increasing voltage/throttle grid ending at full throttle")
    values = np.full((len(voltages), len(throttles), len(CURVE_FIELDS)), np.nan)
    values[:, 0, :] = 0.0
    for row in rows:
        i = int(np.searchsorted(voltages, float(row["voltage_v"])))
        j = int(np.searchsorted(throttles, float(row["throttle"])))
        if np.all(np.isfinite(values[i, j])):
            raise ValueError("Duplicate lift performance sample")
        values[i, j] = [float(row[field]) for field in CURVE_FIELDS]
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Lift performance grid must be complete, finite and nonnegative")
    if np.any(np.diff(values[:, :, 0], axis=1) <= 0):
        raise ValueError("Lift thrust must increase with throttle")
    return StaticPerformance(voltages, throttles, values)


def load_static_performance(path: str | Path) -> StaticPerformance:
    path = Path(path).resolve()
    info = path.stat()
    return _load_performance(str(path), info.st_mtime_ns, info.st_size)


def propulsion_model_placeholder(control: Control, params: AircraftParams) -> tuple[np.ndarray, np.ndarray]:
    """Keep the reference model's preliminary forward-thrust mapping."""
    throttle = np.clip(control.throttle, 0.0, 1.0)
    thrust = params.Tmax * throttle if params.pusher_installed else 0.0
    return np.array([thrust, 0.0, 0.0]), np.zeros(3)


def pusher_performance(control: Control, params: AircraftParams, axial_speed_m_s: float = 0.0):
    if not params.pusher_installed:
        return {"thrust_n": 0.0}
    if params.pusher_system is None:
        return {"thrust_n": float(propulsion_model_placeholder(control, params)[0][0])}
    from .pusher import estimate_pusher
    return estimate_pusher(params.pusher_system, params.battery.voltage_v, control.throttle,
                           axial_speed_m_s, params.rho)


def pusher_propulsion(control: Control, params: AircraftParams, axial_speed_m_s: float = 0.0):
    sample = pusher_performance(control, params, axial_speed_m_s)
    force = np.array([sample["thrust_n"], 0., 0.])
    moment = np.zeros(3)
    if params.pusher_installed and params.pusher_system is not None:
        system = params.pusher_system
        if system.position_body_m is not None:
            moment += np.cross(np.asarray(system.position_body_m)-params.cg_body_xyz, force)
        moment[0] += system.reaction_torque_sign*sample["torque_nm"]
    return force, moment


def lift_propulsion(control: Control, params: AircraftParams) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]]]:
    """Resolve rotor thrust and reaction torque about the aircraft CG."""
    lift = params.lift_system
    if lift is None:
        if control.lift_throttle != 0.0 or control.lift_motor_throttles is not None:
            raise ValueError("This aircraft has no configured lift system")
        return np.zeros(3), np.zeros(3), []
    commands = control.lift_motor_throttles
    if commands is None:
        commands = (control.lift_throttle,) * len(lift.positions_body_m)
    if len(commands) != len(lift.positions_body_m) or not np.all(np.isfinite(commands)):
        raise ValueError("Supply one finite lift throttle per configured rotor")
    if lift.performance_path is None:
        if any(command != 0.0 for command in commands):
            raise ValueError(f"No lift performance curve for {params.aircraft_label}; supply measured data before running lift motors")
        return np.zeros(3), np.zeros(3), []
    curve = load_static_performance(lift.performance_path)
    force = np.zeros(3)
    moment = np.zeros(3)
    samples = []
    for position, sign, command in zip(lift.positions_body_m, lift.reaction_torque_signs, commands):
        sample = curve.lookup(params.battery.voltage_v, command)
        rotor_force = np.array([0.0, 0.0, -sample["thrust_n"]])
        force += rotor_force
        moment += np.cross(np.asarray(position) - params.cg_body_xyz, rotor_force)
        moment[2] += sign * sample["torque_nm"]
        samples.append(sample)
    return force, moment, samples


def propulsion_model(control: Control, params: AircraftParams, axial_speed_m_s: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Combine independently commanded forward and lift propulsion."""
    pusher_force, pusher_moment = pusher_propulsion(control, params, axial_speed_m_s)
    lift_force, lift_moment, _ = lift_propulsion(control, params)
    return pusher_force + lift_force, pusher_moment + lift_moment


def hover_collective(params: AircraftParams) -> float:
    """Estimate collective for static weight support, without an altitude controller."""
    lift = params.lift_system
    if lift is None or lift.performance_path is None:
        raise ValueError("Hover collective requires a lift performance curve")
    curve = load_static_performance(lift.performance_path)
    thrust = np.array([curve.lookup(params.battery.voltage_v, t)["thrust_n"] for t in curve.throttles])
    per_rotor_weight = params.mass * params.g / len(lift.positions_body_m)
    if per_rotor_weight > thrust[-1]:
        raise ValueError("Aircraft weight exceeds the configured static lift capacity")
    return float(np.interp(per_rotor_weight, thrust, curve.throttles))


def hover_trim_control(params: AircraftParams) -> Control:
    """Solve physical per-rotor static hover loads for an explicit CG trial.

    This is an initial-condition diagnostic, not a mixer or feedback tune.
    Firmware still owns all actual motor outputs during a SITL run.
    """
    from scipy.optimize import least_squares
    collective = hover_collective(params)
    count = len(params.lift_system.positions_body_m)

    def control(commands):
        return Control(0, 0, 0, 0, lift_throttle=collective,
                       lift_motor_throttles=tuple(float(v) for v in commands))

    def residual(commands):
        force, moment, _ = lift_propulsion(control(commands), params)
        return np.r_[force[2]+params.mass*params.g, moment]

    fit = least_squares(residual, np.full(count, collective), bounds=(0, 1),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    if not fit.success or np.linalg.norm(residual(fit.x)) > 1e-6:
        raise ValueError("Configured CG has no balanced per-rotor hover trim")
    return control(fit.x)


def propulsion_telemetry(control: Control, params: AircraftParams, axial_speed_m_s: float = 0.0) -> dict[str, float]:
    """Report static bench estimates, without a battery or voltage-sag model."""
    _, _, samples = lift_propulsion(control, params)
    pusher = pusher_performance(control, params, axial_speed_m_s)
    return {
        "lift_thrust_n": sum(sample["thrust_n"] for sample in samples),
        "lift_current_a": sum(sample["current_a"] for sample in samples),
        "lift_power_w": sum(sample["power_w"] for sample in samples),
        **{f"pusher_{key}": value for key, value in pusher.items()},
    }
