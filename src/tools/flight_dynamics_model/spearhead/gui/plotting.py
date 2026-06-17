"""Plotting helpers for GUI pages."""

from __future__ import annotations

import numpy as np

from .adapters.simulation import simulation_time_history_rows
from ..result import SimulationResult
from ..stability.types import StabilityResult


def _require_plotly():
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as exc:
        raise RuntimeError(
            "Plotly is required for GUI plots. Install optional GUI dependencies with: "
            "python -m pip install -r src/tools/flight_dynamics_model/requirements-gui.txt"
        ) from exc
    return go, make_subplots


def _series(rows: list[dict[str, float]], field: str) -> list[float]:
    return [float(row[field]) for row in rows]


def simulation_time_history_figure(result: SimulationResult):
    """Return a Plotly figure for key simulation time histories."""
    go, make_subplots = _require_plotly()
    rows = simulation_time_history_rows(result)
    time = _series(rows, "time")

    fig = make_subplots(
        rows=4,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=(
            "Body Velocity",
            "Body Rates",
            "Euler Angles",
            "Altitude",
        ),
    )

    for field in ("u", "v", "w"):
        fig.add_trace(go.Scatter(x=time, y=_series(rows, field), mode="lines", name=field), row=1, col=1)

    for field in ("p", "q", "r"):
        values = np.rad2deg(_series(rows, field)).tolist()
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=f"{field} deg/s"), row=2, col=1)

    for field in ("phi", "theta", "psi"):
        values = np.rad2deg(_series(rows, field)).tolist()
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=f"{field} deg"), row=3, col=1)

    altitude = [-value for value in _series(rows, "pd")]
    fig.add_trace(go.Scatter(x=time, y=altitude, mode="lines", name="altitude"), row=4, col=1)

    fig.update_layout(
        height=820,
        margin={"l": 40, "r": 20, "t": 60, "b": 40},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "left", "x": 0},
    )
    fig.update_xaxes(title_text="Time [s]", row=4, col=1)
    fig.update_yaxes(title_text="m/s", row=1, col=1)
    fig.update_yaxes(title_text="deg/s", row=2, col=1)
    fig.update_yaxes(title_text="deg", row=3, col=1)
    fig.update_yaxes(title_text="m", row=4, col=1)
    return fig


def stability_eigenvalue_figure(result: StabilityResult):
    """Return a Plotly pole plot for longitudinal and lateral-directional modes."""
    go, _ = _require_plotly()
    fig = go.Figure()
    for analysis in (result.longitudinal, result.lateral_directional):
        eigenvalues = [mode.eigenvalue for mode in analysis.modes]
        fig.add_trace(
            go.Scatter(
                x=[float(eigenvalue.real) for eigenvalue in eigenvalues],
                y=[float(eigenvalue.imag) for eigenvalue in eigenvalues],
                mode="markers+text",
                name=analysis.subsystem,
                text=[mode.mode for mode in analysis.modes],
                textposition="top center",
            )
        )

    fig.add_hline(y=0.0, line_width=1, line_dash="dash", line_color="#94a3b8")
    fig.add_vline(x=0.0, line_width=1, line_dash="dash", line_color="#ef4444")
    fig.update_layout(
        height=520,
        margin={"l": 50, "r": 20, "t": 40, "b": 45},
        xaxis_title="Real [1/s]",
        yaxis_title="Imaginary [rad/s]",
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "left", "x": 0},
    )
    return fig


def figure_to_json_dict(figure) -> dict:
    """Return a serializable Plotly figure dictionary."""
    return figure.to_plotly_json()
