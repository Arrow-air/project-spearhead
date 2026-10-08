"""Air-start trim, shared ODE stepping, and body-FRD/NED sensor conversion."""

import math
import numpy as np
from scipy.optimize import least_squares

from ..controls import Control
from ..dynamics import dynamics
from ..force_moment import compute_force_moment_breakdown
from ..rotations import euler_body_to_ned


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

    bounds = (np.r_[np.deg2rad([-5, -10, -10, -5, -5, -5]), 0],
              np.r_[np.deg2rad([0, 10, 10, 5, 5, 5]), 1])
    fit = least_squares(residual, np.r_[np.deg2rad([-.15, -.15, 0, -3, 0, 0]), .25],
                        bounds=bounds, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    x, control = unpack(fit.x)
    loads = compute_force_moment_breakdown(0, x, params, control)
    if (not fit.success or np.linalg.norm(loads["total"]["force_b"]) >= .01
            or np.linalg.norm(loads["total"]["moment_b"]) >= .01
            or abs(dynamics(0, x, params, control)[2]) >= .001):
        raise ValueError("No full six-axis trim within measured control coverage")
    return x, control


def sensor_state(t, x, params, control, *, hold=False):
    rotation = euler_body_to_ned(*x[6:9])
    if hold:
        accel = -rotation.T @ np.array([0., 0., params.g])
    else:
        loads = compute_force_moment_breakdown(t, x, params, control)
        accel = (loads["total"]["force_b"] - loads["gravity"]["force_b"]) / params.mass
    # Airspeed is omitted: 4.7.1 derives EAS/pitot from velocity and atmosphere.
    # No TAS value is injected into the pitot/EAS field at the wrong density.
    return {"timestamp": t, "imu": {"gyro": x[9:12].tolist(), "accel_body": accel.tolist()},
            "position": x[:3].tolist(), "velocity": (rotation @ x[3:6]).tolist(),
            "attitude": x[6:9].tolist()}


class Physics:
    def __init__(self, params, state, *, hold=True, max_step=.0025):
        if not np.isfinite(state).all() or np.shape(state) != (12,):
            raise ValueError("State must contain 12 finite values")
        if not math.isfinite(max_step) or not 0 < max_step <= .0025:
            raise ValueError("Physics substep must be in (0, 0.0025] seconds")
        self.params, self.state = params, np.array(state, dtype=float, copy=True)
        self.time, self.hold, self.max_step = 0., hold, max_step

    def advance(self, dt, control):
        if not math.isfinite(dt) or not 0 < dt <= 1:
            raise ValueError("Invalid physics timestep")
        x = self.state.copy()
        if self.hold:
            x[:3] += (euler_body_to_ned(*x[6:9]) @ x[3:6]) * dt
        else:
            n = math.ceil(dt / self.max_step)
            h = dt / n
            for i in range(n):
                t = self.time + i*h
                k1 = dynamics(t, x, self.params, control)
                k2 = dynamics(t+h/2, x+h*k1/2, self.params, control)
                k3 = dynamics(t+h/2, x+h*k2/2, self.params, control)
                k4 = dynamics(t+h, x+h*k3, self.params, control)
                x += h*(k1+2*k2+2*k3+k4)/6
                if not np.isfinite(x).all() or abs(x[7]) >= math.radians(85) or x[2] >= 0:
                    raise ValueError("Physics diverged, reached ground, or approached Euler singularity")
        self.state, self.time = x, self.time + dt
        return sensor_state(self.time, x, self.params, control, hold=self.hold)
