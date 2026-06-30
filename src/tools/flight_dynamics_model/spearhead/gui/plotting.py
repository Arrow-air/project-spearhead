"""Plotting helpers for GUI pages."""

from __future__ import annotations

import numpy as np

from .adapters.simulation import simulation_time_history_rows
from .adapters.stability import cg_sweep_summary_rows
from ..result import SimulationResult
from ..stability.types import CGSweepResult, StabilityResult


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


def _optional_series(rows: list[dict[str, float]], field: str, default: float = 0.0) -> list[float]:
    return [float(row.get(field, default)) for row in rows]


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


def _flight_deck_base_figure(title: str):
    go, _ = _require_plotly()
    fig = go.Figure()
    fig.update_layout(
        title={"text": title, "font": {"size": 13}},
        height=255,
        margin={"l": 42, "r": 14, "t": 38, "b": 34},
        paper_bgcolor="#0f172a",
        plot_bgcolor="#111827",
        font={"color": "#dbeafe", "size": 11},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "left", "x": 0},
        hovermode="x unified",
    )
    fig.update_xaxes(title_text="Time [s]", gridcolor="#243247", zerolinecolor="#334155")
    fig.update_yaxes(gridcolor="#243247", zerolinecolor="#334155")
    return go, fig


def flight_deck_attitude_figure(history: list[dict[str, float]]):
    """Return the rolling attitude plot for the Flight Deck."""
    go, fig = _flight_deck_base_figure("Attitude")
    time = _series(history, "time") if history else [0.0]
    for field in ("phi", "theta", "psi"):
        values = np.rad2deg(_series(history, field)).tolist() if history else [0.0]
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=f"{field} deg"))
    fig.update_yaxes(title_text="deg")
    return fig


def flight_deck_rates_figure(history: list[dict[str, float]]):
    """Return the rolling body-rates plot for the Flight Deck."""
    go, fig = _flight_deck_base_figure("Rates")
    time = _series(history, "time") if history else [0.0]
    for field in ("p", "q", "r"):
        values = np.rad2deg(_series(history, field)).tolist() if history else [0.0]
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=f"{field} deg/s"))
    fig.update_yaxes(title_text="deg/s")
    return fig


def flight_deck_controls_figure(history: list[dict[str, float]]):
    """Return the rolling control-input plot for the Flight Deck."""
    go, fig = _flight_deck_base_figure("Control Inputs")
    time = _series(history, "time") if history else [0.0]
    fields = (
        ("elevator_deg", "elevator deg"),
        ("aileron_deg", "aileron deg"),
        ("rudder_deg", "rudder deg"),
        ("throttle_pct", "throttle %"),
    )
    for field, name in fields:
        values = _optional_series(history, field) if history else [0.0]
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=name))
    fig.update_yaxes(title_text="deg / %")
    return fig


def flight_deck_energy_figure(history: list[dict[str, float]]):
    """Return the rolling airspeed/altitude plot for the Flight Deck."""
    go, fig = _flight_deck_base_figure("Airspeed / Altitude")
    time = _series(history, "time") if history else [0.0]
    if history:
        airspeed = [
            float(np.sqrt(row["u"] * row["u"] + row["v"] * row["v"] + row["w"] * row["w"]))
            for row in history
        ]
        altitude = [-float(row["pd"]) for row in history]
    else:
        airspeed = [0.0]
        altitude = [0.0]
    fig.add_trace(go.Scatter(x=time, y=airspeed, mode="lines", name="airspeed m/s"))
    fig.add_trace(go.Scatter(x=time, y=altitude, mode="lines", name="altitude m"))
    fig.update_yaxes(title_text="m/s / m")
    return fig


def flight_deck_live_figures(history: list[dict[str, float]]) -> dict[str, object]:
    """Return Flight Deck live plot cards keyed by card name."""
    return {
        "attitude": flight_deck_attitude_figure(history),
        "rates": flight_deck_rates_figure(history),
        "controls": flight_deck_controls_figure(history),
        "airspeed_altitude": flight_deck_energy_figure(history),
    }


def flight_deck_live_figure(history: list[dict[str, float]]):
    """Return a legacy combined Flight Deck plot."""
    go, make_subplots = _require_plotly()
    time = _series(history, "time") if history else [0.0]

    fig = make_subplots(
        rows=4,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.07,
        subplot_titles=("Attitude", "Rates", "Control Inputs", "Airspeed / Altitude"),
    )

    for field in ("phi", "theta", "psi"):
        values = np.rad2deg(_series(history, field)).tolist() if history else [0.0]
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=f"{field} deg"), row=1, col=1)
    for field in ("p", "q", "r"):
        values = np.rad2deg(_series(history, field)).tolist() if history else [0.0]
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=f"{field} deg/s"), row=2, col=1)
    for field, name in (
        ("elevator_deg", "elevator deg"),
        ("aileron_deg", "aileron deg"),
        ("rudder_deg", "rudder deg"),
        ("throttle_pct", "throttle %"),
    ):
        values = _optional_series(history, field) if history else [0.0]
        fig.add_trace(go.Scatter(x=time, y=values, mode="lines", name=name), row=3, col=1)
    if history:
        airspeed = [
            float(np.sqrt(row["u"] * row["u"] + row["v"] * row["v"] + row["w"] * row["w"]))
            for row in history
        ]
        altitude = [-float(row["pd"]) for row in history]
    else:
        airspeed = [0.0]
        altitude = [0.0]
    fig.add_trace(go.Scatter(x=time, y=airspeed, mode="lines", name="airspeed m/s"), row=4, col=1)
    fig.add_trace(go.Scatter(x=time, y=altitude, mode="lines", name="altitude m"), row=4, col=1)

    fig.update_layout(
        height=620,
        margin={"l": 42, "r": 14, "t": 44, "b": 34},
        paper_bgcolor="#0f172a",
        plot_bgcolor="#111827",
        font={"color": "#dbeafe", "size": 11},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "left", "x": 0},
    )
    fig.update_xaxes(title_text="Time [s]", gridcolor="#243247", zerolinecolor="#334155", row=4, col=1)
    fig.update_yaxes(gridcolor="#243247", zerolinecolor="#334155")
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


def cg_sweep_trend_figure(result: CGSweepResult):
    """Return a Plotly trend of CG x against the most unstable pole real part."""
    go, _ = _require_plotly()
    rows = [
        row
        for row in cg_sweep_summary_rows(result)
        if row["success"] and row["most_unstable_real"] is not None
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[row["cg_x"] for row in rows],
            y=[row["most_unstable_real"] for row in rows],
            mode="lines+markers",
            name="max real eigenvalue",
        )
    )
    fig.add_hline(y=0.0, line_width=1, line_dash="dash", line_color="#ef4444")
    fig.update_layout(
        height=420,
        margin={"l": 50, "r": 20, "t": 40, "b": 45},
        xaxis_title="CG x body [m]",
        yaxis_title="Max real eigenvalue [1/s]",
    )
    return fig


def figure_to_json_dict(figure) -> dict:
    """Return a serializable Plotly figure dictionary."""
    return figure.to_plotly_json()
