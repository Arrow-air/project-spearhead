"""Air-start trim, shared ODE stepping, and body-FRD/NED sensor conversion."""

import math
import numpy as np
from scipy.optimize import least_squares

from ..controls import Control
from ..dynamics import dynamics
from ..force_moment import compute_force_moment_breakdown
from ..rotations import euler_body_to_ned
from . import SITL_HOME, ardupilot_isa_density


def common_control_alpha_bounds(params):
    """Intersect the actual base/control coverage instead of freezing old limits."""
    if params.deflection_database_path is None:
        return (-5., 0.)
    from ..aerodb import load_aerodynamic_database_cached
    from ..deflections import load_deflection_database_cached
    base = load_aerodynamic_database_cached(params.aero_database_path, params.aerodb_reference)
    controls = load_deflection_database_cached(params.deflection_database_path, params.aerodb_reference)
    lower = max(base.alpha_bounds_deg[0], *(s.alpha_deg[0] for s in controls.surfaces.values()))
    upper = min(base.alpha_bounds_deg[1], *(s.alpha_deg[-1] for s in controls.surfaces.values()))
    if lower >= upper:
        raise ValueError("No common measured alpha range for full trim")
    return float(lower), float(upper)


def full_trim(params, speed=28., altitude=100.):
    if not np.isfinite([speed, altitude]).all() or speed <= 0 or altitude <= 0:
        raise ValueError("Air-start requires positive finite speed and altitude")

    def unpack(z):
        alpha, theta, phi, de, da, dr, throttle = z
        x = np.array([0., 0., -altitude, speed*np.cos(alpha), 0., speed*np.sin(alpha),
                      phi, theta, 0., 0., 0., 0.])
        return x, Control(de, da, dr, throttle)

    def residual(z):
        x, control = unpack(z)
        dx = dynamics(0, x, params, control)
        return np.r_[dx[3:6], dx[9:12], dx[2]/speed]

    limits = np.array(params.control_surface_limits_deg) if params.control_deflection_policy == "linear" else np.full(3, 5.)
    alpha_min, alpha_max = common_control_alpha_bounds(params)
    bounds = (np.r_[np.deg2rad([alpha_min, -10, -10]), -np.deg2rad(limits), 0],
              np.r_[np.deg2rad([alpha_max, 10, 10]), np.deg2rad(limits), 1])
    initial_alpha = float(np.clip(-.15, alpha_min, alpha_max))
    initial_throttle = .85 if params.pusher_system is not None else .25
    # The elevator branch is strongly nonlinear outside measured coverage; a single
    # start can stall between solvable neighbours. Try a few elevator starts in order.
    starts = [-3.] + [value for value in (-12., -20., -28.) if value >= -limits[0]]
    for start_de in starts:
        fit = least_squares(residual, np.r_[np.deg2rad([initial_alpha, -.15, 0, start_de, 0, 0]), initial_throttle],
                            bounds=bounds, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        x, control = unpack(fit.x)
        loads = compute_force_moment_breakdown(0, x, params, control)
        if (fit.success and np.linalg.norm(loads["total"]["force_b"]) < .01
                and np.linalg.norm(loads["total"]["moment_b"]) < .01
                and abs(dynamics(0, x, params, control)[2]) < .001):
            return x, control
    raise ValueError("No full six-axis trim within configured surface travel and common measured alpha range")


def sensor_state(t, x, params, control, *, hold=False, ground=None, battery=None, wind_ned=None):
    rotation = euler_body_to_ned(*x[6:9])
    if hold:
        accel = -rotation.T @ np.array([0., 0., params.g])
    else:
        loads = compute_force_moment_breakdown(t, x, params, control, wind_ned=wind_ned)
        accel = (loads["total"]["force_b"] - loads["gravity"]["force_b"]) / params.mass
        if ground is not None:
            accel += ground.loads(x)[0] / params.mass
    sample = {"timestamp": t, "imu": {"gyro": x[9:12].tolist(), "accel_body": accel.tolist()},
              "position": x[:3].tolist(), "velocity": (rotation @ x[3:6]).tolist(),
              "attitude": x[6:9].tolist(), "airspeed": ardupilot_pitot_airspeed(x, params, wind_ned)}
    if battery is not None:
        sample["battery"] = battery
    return sample


def ardupilot_pitot_airspeed(x, params, wind_ned=None):
    """JSON airspeed that makes ArduPilot indicate the FDM's equivalent airspeed.

    Without this field, 4.7.1 JSON passes body TAS to its pitot model and converts it
    with its own standard-atmosphere density at the vehicle AMSL altitude, not the FDM
    density. 4.7.1 treats this field as true airspeed and divides it by that same ISA
    EAS2TAS, so scaling by sqrt(rho_FDM/rho_ISA) yields TAS*sqrt(rho_FDM/1.225). The
    20-degree pitot cone and cosine roll-off copy Aircraft::update_eas_airspeed.
    """
    v_air = np.asarray(x[3:6], dtype=float)
    if wind_ned is not None:
        v_air = v_air - euler_body_to_ned(*x[6:9]).T @ np.asarray(wind_ned, dtype=float)
    u, v, w = (float(value) for value in v_air)
    pitot = math.sqrt(u*u + v*v + w*w)
    pitot_aoa = math.atan2(math.hypot(v, w), u)
    max_pitot_aoa = math.radians(20)
    if pitot_aoa > math.radians(90):
        pitot = 0.
    elif pitot_aoa > max_pitot_aoa:
        pitot *= math.cos((pitot_aoa - max_pitot_aoa)*(math.pi/2)/(math.radians(90) - max_pitot_aoa))
    altitude_amsl = SITL_HOME[2] - float(x[2])
    return pitot*math.sqrt(params.rho/ardupilot_isa_density(altitude_amsl))


class Physics:
    def __init__(self, params, state, *, hold=True, max_step=.0025, ground=None, battery=None, wind=None):
        if not np.isfinite(state).all() or np.shape(state) != (12,):
            raise ValueError("State must contain 12 finite values")
        if not math.isfinite(max_step) or not 0 < max_step <= .0025:
            raise ValueError("Physics substep must be in (0, 0.0025] seconds")
        self.params, self.state = params, np.array(state, dtype=float, copy=True)
        self.time, self.hold, self.max_step = 0., hold, max_step
        self.ground = ground
        self.battery = battery
        # Wind starts at release so the stationary pitot calibration sees still air.
        self.wind = None if wind is None or wind.calm else wind
        self.release_time = None

    def wind_ned(self, t):
        if self.wind is None or self.release_time is None:
            return None
        return self.wind.ned(t - self.release_time)

    def derivative(self, t, x, control):
        result = dynamics(t, x, self.params, control, wind_ned=self.wind_ned(t))
        if self.ground is not None:
            force, moment, _ = self.ground.loads(x)
            result[3:6] += force/self.params.mass
            result[9:12] += np.linalg.solve(self.params.inertia, moment)
        return result

    def advance(self, dt, control):
        if not math.isfinite(dt) or not 0 < dt <= 1:
            raise ValueError("Invalid physics timestep")
        x = self.state.copy()
        if not self.hold and self.release_time is None:
            # Air starts begin at release. Ground starts begin once the lift motors
            # spin (after arming), i.e. after the firmware's pitot calibration.
            spinning = control.lift_throttle > 0 or any(v > 0 for v in (control.lift_motor_throttles or ()))
            if self.ground is None or spinning:
                self.release_time = self.time
        if self.hold:
            x[:3] += (euler_body_to_ned(*x[6:9]) @ x[3:6]) * dt
        else:
            n = math.ceil(dt / self.max_step)
            h = dt / n
            for i in range(n):
                t = self.time + i*h
                k1 = self.derivative(t, x, control)
                k2 = self.derivative(t+h/2, x+h*k1/2, control)
                k3 = self.derivative(t+h/2, x+h*k2/2, control)
                k4 = self.derivative(t+h, x+h*k3, control)
                x += h*(k1+2*k2+2*k3+k4)/6
                if not np.isfinite(x).all() or abs(x[7]) >= math.radians(85) or x[2] >= (0 if self.ground is None else .05):
                    raise ValueError("Physics diverged, reached ground, or approached Euler singularity")
        self.state, self.time = x, self.time + dt
        reading = None
        if self.battery is not None:
            # Lift commands are held over the frame, so the frame current is constant.
            current = self.battery.pack_current(control, self.params)
            self.battery.consume(current, dt)
            reading = self.battery.reading(current)
        return sensor_state(self.time, x, self.params, control, hold=self.hold, ground=self.ground,
                            battery=reading, wind_ned=self.wind_ned(self.time))
