"""Plotting helpers for open-loop simulation results."""

import matplotlib.pyplot as plt
import numpy as np

from .aero import aero_model
from .controls import ControlInput, control_at
from .force_moment import COMPONENTS
from .params import AircraftParams


def _airspeed_alpha_beta(
    sol,
    params: AircraftParams,
    control_override: ControlInput | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    V = np.zeros_like(sol.t)
    alpha = np.zeros_like(sol.t)
    beta = np.zeros_like(sol.t)

    for i, t in enumerate(sol.t):
        control = control_at(float(t), control_override)
        _, _, info = aero_model(sol.y[3:6, i], sol.y[9:12, i], control, params)
        V[i] = info["V"]
        alpha[i] = info["alpha"]
        beta[i] = info["beta"]

    return V, alpha, beta


def _control_history(
    t: np.ndarray,
    control_override: ControlInput | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    de = np.zeros_like(t)
    da = np.zeros_like(t)
    dr = np.zeros_like(t)
    throttle = np.zeros_like(t)

    for i, ti in enumerate(t):
        control = control_at(float(ti), control_override)
        de[i] = control.de
        da[i] = control.da
        dr[i] = control.dr
        throttle[i] = control.throttle

    return de, da, dr, throttle


def plot_open_loop(
    sol,
    params: AircraftParams | None = None,
    control_override: ControlInput | None = None,
):
    """Create standard plots for an open-loop simulation."""
    if params is None:
        params = AircraftParams()
    if control_override is None:
        control_override = getattr(sol, "control_override", None)

    t = sol.t
    pn, pe, pd = sol.y[0], sol.y[1], sol.y[2]
    phi, theta, psi = sol.y[6], sol.y[7], sol.y[8]
    p, q, r = sol.y[9], sol.y[10], sol.y[11]
    V, alpha, beta = _airspeed_alpha_beta(sol, params, control_override=control_override)
    de, da, dr, throttle = _control_history(t, control_override=control_override)

    fig, axes = plt.subplots(4, 2, figsize=(12, 13), constrained_layout=True)

    axes[0, 0].plot(t, V)
    axes[0, 0].set_title("Airspeed")
    axes[0, 0].set_xlabel("Time [s]")
    axes[0, 0].set_ylabel("V [m/s]")
    axes[0, 0].grid(True)

    axes[0, 1].plot(t, -pd)
    axes[0, 1].set_title("Altitude")
    axes[0, 1].set_xlabel("Time [s]")
    axes[0, 1].set_ylabel("Altitude [m]")
    axes[0, 1].grid(True)

    axes[1, 0].plot(t, np.rad2deg(alpha), label="alpha")
    axes[1, 0].plot(t, np.rad2deg(beta), label="beta")
    axes[1, 0].set_title("Alpha / Beta")
    axes[1, 0].set_xlabel("Time [s]")
    axes[1, 0].set_ylabel("Angle [deg]")
    axes[1, 0].legend()
    axes[1, 0].grid(True)

    axes[1, 1].plot(t, np.rad2deg(phi), label="phi")
    axes[1, 1].plot(t, np.rad2deg(theta), label="theta")
    axes[1, 1].plot(t, np.rad2deg(psi), label="psi")
    axes[1, 1].set_title("Euler Angles")
    axes[1, 1].set_xlabel("Time [s]")
    axes[1, 1].set_ylabel("Angle [deg]")
    axes[1, 1].legend()
    axes[1, 1].grid(True)

    axes[2, 0].plot(t, np.rad2deg(p), label="p")
    axes[2, 0].plot(t, np.rad2deg(q), label="q")
    axes[2, 0].plot(t, np.rad2deg(r), label="r")
    axes[2, 0].set_title("Body Rates")
    axes[2, 0].set_xlabel("Time [s]")
    axes[2, 0].set_ylabel("Rate [deg/s]")
    axes[2, 0].legend()
    axes[2, 0].grid(True)

    ax_controls = axes[3, 0]
    ax_controls.plot(t, np.rad2deg(de), label="de")
    ax_controls.plot(t, np.rad2deg(da), label="da")
    ax_controls.plot(t, np.rad2deg(dr), label="dr")
    ax_controls.set_title("Control Inputs")
    ax_controls.set_xlabel("Time [s]")
    ax_controls.set_ylabel("Deflection [deg]")
    ax_controls.grid(True)

    ax_throttle = ax_controls.twinx()
    ax_throttle.plot(t, throttle, color="black", linestyle="--", label="throttle")
    ax_throttle.set_ylabel("Throttle [-]")

    control_lines, control_labels = ax_controls.get_legend_handles_labels()
    throttle_lines, throttle_labels = ax_throttle.get_legend_handles_labels()
    ax_controls.legend(control_lines + throttle_lines, control_labels + throttle_labels)

    fig.delaxes(axes[2, 1])
    fig.delaxes(axes[3, 1])
    ax3d = fig.add_subplot(4, 2, 8, projection="3d")
    ax3d.plot(pn, pe, -pd)
    ax3d.set_title("3D Trajectory")
    ax3d.set_xlabel("North [m]")
    ax3d.set_ylabel("East [m]")
    ax3d.set_zlabel("Altitude [m]")

    return fig


def plot_force_moment_breakdown(history):
    """Plot body-frame force and moment component histories."""
    t = history["time"]
    components = list(COMPONENTS)
    force_labels = ("Fx [N]", "Fy [N]", "Fz [N]")
    moment_labels = ("Mx [N*m]", "My [N*m]", "Mz [N*m]")
    titles = ("Body Force Fx", "Body Force Fy", "Body Force Fz", "Body Moment Mx", "Body Moment My", "Body Moment Mz")

    fig, axes = plt.subplots(3, 2, figsize=(13, 11), constrained_layout=True)
    axes_flat = axes.ravel()

    for axis_index in range(3):
        ax = axes_flat[axis_index]
        for component in components:
            ax.plot(t, history[component]["force_b"][:, axis_index], label=component)
        ax.set_title(titles[axis_index])
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(force_labels[axis_index])
        ax.grid(True)
        ax.legend()

    for axis_index in range(3):
        ax = axes_flat[axis_index + 3]
        for component in components:
            ax.plot(t, history[component]["moment_b"][:, axis_index], label=component)
        ax.set_title(titles[axis_index + 3])
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(moment_labels[axis_index])
        ax.grid(True)
        ax.legend()

    return fig
