"""Minimal NiceGUI shell for the Spearhead Flight Dynamics Model."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from .adapters.exports import (
    export_simulation_csv_text,
    export_simulation_json_text,
    export_stability_json_text,
    export_stability_markdown_text,
    export_stability_modes_csv_text,
)
from .adapters.scenarios import load_scenario_config, scenario_summary
from .adapters.simulation import (
    build_simulation_config,
    run_simulation_for_gui,
    simulation_result_summary,
)
from .adapters.stability import (
    build_linearization_config,
    cg_sweep_summary_rows,
    run_cg_sweep_for_gui,
    run_stability_for_gui,
    stability_mode_table_rows,
    static_derivative_table_rows,
)
from .flight_deck import (
    COMPLETE,
    ERROR,
    PAUSED,
    READY,
    RUNNING,
    FlightDeckSession,
    attitude_html,
    maneuver_options,
)
from .plotting import (
    cg_sweep_trend_figure,
    flight_deck_live_figures,
    simulation_time_history_figure,
    stability_eigenvalue_figure,
)
from .state import APP_STATE
from ..config import ControlInputCommand, SimulationConfig


FDM_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_SCENARIOS_DIR = FDM_ROOT / "scenarios"
DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8080

_pages_registered = False
_flight_deck_session: FlightDeckSession | None = None


def get_example_scenarios() -> list[Path]:
    """Return bundled YAML scenarios that can be loaded from the GUI."""
    return sorted(EXAMPLE_SCENARIOS_DIR.glob("*.yaml"))


def _require_nicegui():
    try:
        from nicegui import ui
    except ImportError as exc:
        raise RuntimeError(
            "NiceGUI is required to run the Spearhead GUI. "
            "Install optional GUI dependencies with: "
            "python -m pip install -r src/tools/flight_dynamics_model/requirements-gui.txt"
        ) from exc
    return ui


def _summary_rows(summary: dict[str, Any]) -> list[dict[str, str]]:
    rows = []
    for key, value in summary.items():
        if isinstance(value, (dict, list)):
            display_value = json.dumps(value, indent=2)
        else:
            display_value = "" if value is None else str(value)
        rows.append({"field": key, "value": display_value})
    return rows


def _essential_scenario_rows(summary: dict[str, Any]) -> list[dict[str, str]]:
    """Return the first-pass scenario fields most users need."""
    timing = f"{summary['duration']} s"
    if summary["dt"] is not None:
        timing = f"{timing}, dt={summary['dt']} s"
    elif summary["n_points"] is not None:
        timing = f"{timing}, samples={summary['n_points']}"
    return [
        {"field": "name", "value": str(summary["name"])},
        {"field": "timing", "value": timing},
        {"field": "trim", "value": f"{summary['trim_target_speed']} m/s at {summary['trim_altitude']} m"},
        {"field": "start_from_trim", "value": str(summary["start_from_trim"])},
        {"field": "input_commands", "value": str(len(summary["input_commands"]))},
    ]


def _scenario_options() -> dict[str, Path]:
    return {path.name: path for path in get_example_scenarios()}


def _load_uploaded_config(upload_event) -> tuple[str, dict[str, Any]]:
    suffix = Path(getattr(upload_event, "name", "scenario.yaml")).suffix or ".yaml"
    content = upload_event.content.read()
    if isinstance(content, str):
        content = content.encode()

    temp_path = None
    try:
        with tempfile.NamedTemporaryFile("wb", suffix=suffix, delete=False) as stream:
            stream.write(content)
            temp_path = Path(stream.name)
        config = load_scenario_config(temp_path)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)

    source = getattr(upload_event, "name", "uploaded scenario")
    APP_STATE.set_scenario(source, config)
    return source, scenario_summary(config)


def _load_example_config(name: str) -> SimulationConfig:
    scenario_options = _scenario_options()
    if name not in scenario_options:
        raise ValueError(f"Unknown example scenario: {name}")
    config = load_scenario_config(scenario_options[name])
    APP_STATE.set_scenario(name, config)
    return config


def _active_or_example_config(example_name: str | None = None) -> SimulationConfig:
    if APP_STATE.config is not None:
        return APP_STATE.config
    scenario_options = _scenario_options()
    selected_name = example_name or next(iter(scenario_options), None)
    if selected_name is None:
        raise ValueError(f"No example scenarios found in {EXAMPLE_SCENARIOS_DIR}")
    return _load_example_config(selected_name)


def _active_flight_deck_session(example_name: str | None = None) -> FlightDeckSession:
    """Return the shared Flight Deck session, creating it from the active scenario."""
    global _flight_deck_session
    config = _active_or_example_config(example_name)
    source = APP_STATE.scenario_source or config.name
    if _flight_deck_session is None:
        _flight_deck_session = FlightDeckSession(source, config)
    elif _flight_deck_session.base_config is not config:
        _flight_deck_session.set_scenario(source, config)
    return _flight_deck_session


def build_simulation_config_from_form(
    base_config: SimulationConfig,
    *,
    name: str,
    duration: float,
    dt: float | None,
    n_points: int | None,
    start_from_trim: bool,
    trim_target_speed: float,
    trim_altitude: float,
) -> SimulationConfig:
    """Build a simulation config from GUI form values while preserving backend-owned fields."""
    return build_simulation_config(
        name=name,
        duration=duration,
        dt=dt,
        n_points=n_points,
        params=base_config.params,
        initial_state=base_config.initial_state,
        control_override=base_config.control_override,
        input_commands=base_config.input_commands,
        start_from_trim=start_from_trim,
        trim_target_speed=trim_target_speed,
        trim_altitude=trim_altitude,
        rtol=base_config.rtol,
        atol=base_config.atol,
    )


def _command_rows(commands: tuple[ControlInputCommand, ...]) -> list[dict[str, Any]]:
    return [
        {
            "surface": command.surface,
            "kind": command.kind,
            "amplitude": command.amplitude,
            "start_time": command.start_time,
            "duration": command.duration,
        }
        for command in commands
    ]


def _simulation_summary_rows(summary: dict[str, Any]) -> list[dict[str, str]]:
    trim = summary.get("trim") or {}
    rows = [
        {"field": "success", "value": str(summary["success"])},
        {"field": "message", "value": str(summary["message"])},
        {"field": "samples", "value": str(summary["samples"])},
        {"field": "final_time", "value": str(summary["end_time"])},
    ]
    if trim:
        rows.extend(
            [
                {"field": "trim_success", "value": str(trim["success"])},
                {"field": "trim_alpha_deg", "value": f"{trim['alpha_deg']:.6g}"},
                {"field": "trim_theta_deg", "value": f"{trim['theta_deg']:.6g}"},
                {"field": "trim_elevator_deg", "value": f"{trim['de_deg']:.6g}"},
                {"field": "trim_throttle", "value": f"{trim['throttle']:.6g}"},
            ]
        )
    return rows


def stability_result_summary_rows(result) -> list[dict[str, str]]:
    """Return stability result summary rows for GUI tables."""
    trim = result.trim_result
    linearization = result.linearization
    return [
        {"field": "success", "value": "True"},
        {"field": "trim_success", "value": str(bool(trim["success"]))},
        {"field": "trim_message", "value": str(trim["message"])},
        {"field": "trim_alpha_deg", "value": f"{float(trim['alpha_deg']):.6g}"},
        {"field": "trim_theta_deg", "value": f"{float(trim['theta_deg']):.6g}"},
        {"field": "trim_elevator_deg", "value": f"{float(trim['de_deg']):.6g}"},
        {"field": "trim_throttle", "value": f"{float(trim['throttle']):.6g}"},
        {"field": "A_shape", "value": str(tuple(linearization.A.shape))},
        {"field": "B_shape", "value": str(tuple(linearization.B.shape))},
        {"field": "A_longitudinal_shape", "value": str(tuple(linearization.A_longitudinal.shape))},
        {"field": "B_longitudinal_shape", "value": str(tuple(linearization.B_longitudinal.shape))},
        {"field": "A_lateral_directional_shape", "value": str(tuple(linearization.A_lateral_directional.shape))},
        {"field": "B_lateral_directional_shape", "value": str(tuple(linearization.B_lateral_directional.shape))},
        {"field": "longitudinal_state_labels", "value": ", ".join(linearization.longitudinal_state_labels)},
        {"field": "lateral_state_labels", "value": ", ".join(linearization.lateral_state_labels)},
    ]


def parse_cg_x_values(text: str) -> list[float]:
    """Parse comma-separated CG x values from a GUI text input."""
    values: list[float] = []
    for item in text.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        values.append(float(stripped))
    if not values:
        raise ValueError("Enter at least one CG x value.")
    return values


def cg_positions_from_x_values(config: SimulationConfig, x_values: list[float]) -> list[np.ndarray]:
    """Build body-frame CG vectors by varying x and preserving the scenario y/z CG."""
    base_cg = np.asarray(config.params.cg_body_xyz, dtype=float)
    return [np.array([x_value, base_cg[1], base_cg[2]], dtype=float) for x_value in x_values]


def _render_page_shell(ui) -> None:
    ui.colors(primary="#0f766e", secondary="#475569", accent="#38bdf8")
    ui.add_head_html(
        """
        <style>
        body { background: #08111f; }
        .fd-console {
            background: radial-gradient(circle at 50% 20%, #172554 0, #0f172a 38%, #08111f 100%);
            color: #dbeafe;
            min-height: calc(100vh - 58px);
        }
        .fd-panel {
            background: rgba(15, 23, 42, 0.86);
            border: 1px solid rgba(56, 189, 248, 0.18);
            border-radius: 8px;
            box-shadow: 0 16px 40px rgba(0, 0, 0, 0.28);
        }
        .fd-card {
            background: rgba(15, 23, 42, 0.92);
            border: 1px solid rgba(148, 163, 184, 0.20);
            border-radius: 8px;
        }
        .fd-muted { color: #94a3b8; }
        .fd-maneuver-badge {
            display: inline-flex;
            align-items: center;
            gap: 10px;
            min-width: 220px;
            padding: 8px 14px;
            border-radius: 8px;
            border: 1px solid rgba(56, 189, 248, 0.34);
            background: linear-gradient(90deg, rgba(8, 47, 73, 0.92), rgba(20, 83, 45, 0.84));
            box-shadow: 0 0 24px rgba(56, 189, 248, 0.14);
        }
        .fd-maneuver-badge small {
            color: #93c5fd;
            font: 700 10px ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .fd-maneuver-badge strong {
            color: #f8fafc;
            font: 800 20px ui-monospace, SFMono-Regular, Menlo, monospace;
            letter-spacing: 0;
        }
        .fd-pfd {
            position: relative;
            height: 520px;
            min-height: 440px;
            overflow: hidden;
            border-radius: 8px;
            border: 1px solid rgba(125, 211, 252, 0.38);
            background: #020617;
            box-shadow:
                inset 0 0 0 1px rgba(15, 23, 42, 0.9),
                inset 0 0 70px rgba(2, 6, 23, 0.8),
                0 22px 60px rgba(0, 0, 0, 0.36);
        }
        .fd-pfd__bezel {
            position: absolute;
            inset: 12px;
            z-index: 8;
            pointer-events: none;
            border: 1px solid rgba(226, 232, 240, 0.16);
            border-radius: 8px;
            box-shadow: inset 0 0 30px rgba(2, 6, 23, 0.9);
        }
        .fd-pfd__world {
            position: absolute;
            inset: -32%;
            transform-origin: center center;
            transition: transform 90ms linear;
        }
        .fd-pfd__pitch {
            position: absolute;
            inset: -18%;
            transform-origin: center center;
            transition: transform 90ms linear;
        }
        .fd-pfd__sky {
            position: absolute;
            inset: 0 0 50% 0;
            background:
                radial-gradient(circle at 50% 18%, rgba(186, 230, 253, 0.42), transparent 24%),
                linear-gradient(#0284c7 0%, #075985 64%, #0c4a6e 100%);
        }
        .fd-pfd__ground {
            position: absolute;
            inset: 50% 0 0 0;
            background:
                linear-gradient(rgba(251, 191, 36, 0.14), transparent 22%),
                linear-gradient(#92400e 0%, #78350f 54%, #451a03 100%);
        }
        .fd-pfd__horizon {
            position: absolute;
            top: 50%;
            left: 0;
            right: 0;
            height: 3px;
            background: #f8fafc;
            box-shadow: 0 0 14px rgba(248, 250, 252, 0.7);
        }
        .fd-pfd__ladder {
            position: absolute;
            inset: 0;
            z-index: 3;
        }
        .fd-pfd__pitch-mark {
            position: absolute;
            left: 50%;
            display: flex;
            align-items: center;
            gap: 12px;
            transform: translate(-50%, -50%);
            color: rgba(241, 245, 249, 0.9);
            font: 700 13px ui-monospace, SFMono-Regular, Menlo, monospace;
            letter-spacing: 0;
        }
        .fd-pfd__pitch-mark i {
            display: block;
            height: 2px;
            background: rgba(248, 250, 252, 0.92);
            box-shadow: 0 0 8px rgba(248, 250, 252, 0.32);
        }
        .fd-pfd__pitch-label {
            width: 32px;
        }
        .fd-pfd__pitch-label--left {
            text-align: right;
        }
        .fd-pfd__pitch-label--right {
            text-align: left;
        }
        .fd-pfd__roll-scale {
            position: absolute;
            top: 26px;
            left: 50%;
            z-index: 12;
            width: 360px;
            height: 180px;
            transform: translateX(-50%);
            pointer-events: none;
        }
        .fd-pfd__roll-tick {
            position: absolute;
            top: 0;
            left: 50%;
            width: 2px;
            height: 18px;
            transform-origin: 50% 160px;
            background: rgba(226, 232, 240, 0.82);
        }
        .fd-pfd__roll-tick--major {
            height: 26px;
            background: #f8fafc;
        }
        .fd-pfd__roll-tick span {
            position: absolute;
            top: 28px;
            left: 50%;
            transform: translateX(-50%);
            color: #e0f2fe;
            font: 700 11px ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .fd-pfd__roll-zero {
            position: absolute;
            top: -1px;
            left: 50%;
            width: 0;
            height: 0;
            transform: translateX(-50%);
            border-left: 8px solid transparent;
            border-right: 8px solid transparent;
            border-bottom: 13px solid #f8fafc;
        }
        .fd-pfd__roll-pointer {
            position: absolute;
            top: 35px;
            left: 50%;
            z-index: 13;
            width: 0;
            height: 0;
            transform-origin: 50% 151px;
            border-left: 9px solid transparent;
            border-right: 9px solid transparent;
            border-top: 16px solid #facc15;
            filter: drop-shadow(0 0 8px rgba(250, 204, 21, 0.55));
        }
        .fd-pfd__aircraft {
            position: absolute;
            top: 50%;
            left: 50%;
            z-index: 15;
            display: flex;
            align-items: center;
            gap: 12px;
            transform: translate(-50%, -50%);
        }
        .fd-pfd__wing {
            width: 124px;
            height: 7px;
            background: #facc15;
            box-shadow: 0 0 14px rgba(250, 204, 21, 0.66);
        }
        .fd-pfd__wing--left {
            border-radius: 999px 0 0 999px;
        }
        .fd-pfd__wing--right {
            border-radius: 0 999px 999px 0;
        }
        .fd-pfd__aircraft strong {
            width: 18px;
            height: 18px;
            border-radius: 999px;
            border: 4px solid #facc15;
            background: rgba(2, 6, 23, 0.72);
        }
        .fd-pfd__aircraft em {
            position: absolute;
            top: 23px;
            left: 50%;
            width: 5px;
            height: 42px;
            transform: translateX(-50%);
            background: #facc15;
            border-radius: 999px;
            box-shadow: 0 0 12px rgba(250, 204, 21, 0.58);
        }
        .fd-pfd__readout {
            position: absolute;
            left: 18px;
            right: 18px;
            bottom: 16px;
            z-index: 16;
            display: flex;
            justify-content: space-between;
            color: #dbeafe;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            letter-spacing: 0;
        }
        .fd-pfd__readout div {
            display: flex;
            align-items: baseline;
            gap: 7px;
            padding: 7px 10px;
            border: 1px solid rgba(125, 211, 252, 0.22);
            border-radius: 8px;
            background: rgba(2, 6, 23, 0.62);
        }
        .fd-pfd__readout span,
        .fd-pfd__readout small {
            color: #93c5fd;
            font-size: 11px;
            font-weight: 700;
        }
        .fd-pfd__readout strong {
            color: #f8fafc;
            font-size: 20px;
        }
        .fd-telemetry-card {
            min-height: 78px;
            padding: 10px 12px;
            border: 1px solid rgba(148, 163, 184, 0.20);
            border-radius: 8px;
            background: linear-gradient(180deg, rgba(15, 23, 42, 0.96), rgba(2, 6, 23, 0.88));
        }
        .fd-telemetry-card__value {
            color: #f8fafc;
            font: 800 26px ui-monospace, SFMono-Regular, Menlo, monospace;
            line-height: 1;
            letter-spacing: 0;
        }
        .fd-telemetry-card__unit {
            margin-top: 3px;
            color: #7dd3fc;
            font: 700 11px ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .fd-telemetry-card__label {
            margin-top: 8px;
            color: #94a3b8;
            font: 800 11px ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .fd-transport .q-btn { min-width: 88px; }
        .fd-control-meter {
            position: relative;
            height: 12px;
            overflow: hidden;
            border-radius: 999px;
            background: rgba(15, 23, 42, 0.95);
            border: 1px solid rgba(148, 163, 184, 0.25);
        }
        .fd-control-meter__center {
            position: absolute;
            top: 0;
            bottom: 0;
            left: 50%;
            width: 1px;
            background: rgba(226, 232, 240, 0.72);
        }
        .fd-control-meter__fill {
            position: absolute;
            top: 2px;
            bottom: 2px;
            border-radius: 999px;
            background: linear-gradient(90deg, #38bdf8, #5eead4);
            box-shadow: 0 0 10px rgba(94, 234, 212, 0.35);
        }
        .fd-timeline {
            position: relative;
            width: min(360px, 34vw);
            height: 10px;
            overflow: hidden;
            border-radius: 999px;
            border: 1px solid rgba(125, 211, 252, 0.28);
            background: rgba(2, 6, 23, 0.86);
        }
        .fd-timeline__fill {
            height: 100%;
            border-radius: 999px;
            background: linear-gradient(90deg, #38bdf8, #5eead4);
            box-shadow: 0 0 16px rgba(56, 189, 248, 0.42);
        }
        .fd-preset {
            border-left: 3px solid transparent;
            color: #94a3b8;
        }
        .fd-preset--active {
            border-left-color: #5eead4;
            color: #e0f2fe;
            background: rgba(20, 83, 45, 0.24);
        }
        .fd-plot-grid {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 12px;
        }
        .fd-plot-card {
            min-height: 292px;
            padding: 8px;
            border: 1px solid rgba(56, 189, 248, 0.18);
            border-radius: 8px;
            background: rgba(15, 23, 42, 0.72);
        }
        </style>
        """
    )
    with ui.header().classes("items-center justify-between bg-slate-950 text-white"):
        ui.label("Spearhead Flight Deck").classes("text-lg font-semibold")
        with ui.row().classes("items-center gap-2"):
            ui.link("Flight Deck", "/").classes("text-white no-underline")
            ui.link("Analyze", "/analyze").classes("text-white no-underline")
            ui.link("Scenario Library", "/scenario-library").classes("text-white no-underline")
    ui.separator()


def _page_intro(ui, title: str, description: str) -> None:
    ui.label(title).classes("text-2xl font-semibold")
    ui.label(description).classes("text-slate-600")


def _section_header(ui, title: str, description: str | None = None) -> None:
    ui.label(title).classes("text-lg font-semibold")
    if description:
        ui.label(description).classes("text-slate-600")


def _empty_state(ui, text: str):
    return ui.label(text).classes("text-slate-500")


def _status_classes(status: str) -> str:
    colors = {
        READY: "bg-sky-950 text-sky-200 border-sky-700",
        RUNNING: "bg-emerald-950 text-emerald-200 border-emerald-700",
        PAUSED: "bg-amber-950 text-amber-200 border-amber-700",
        COMPLETE: "bg-indigo-950 text-indigo-200 border-indigo-700",
        ERROR: "bg-red-950 text-red-200 border-red-700",
    }
    return colors.get(status, "bg-slate-900 text-slate-200 border-slate-700")


def _control_bar_value(value: float, *, centered: bool) -> float:
    if centered:
        return max(0.0, min(1.0, 0.5 + value / 30.0))
    return max(0.0, min(1.0, value / 100.0))


def _control_bar_html(value: float, *, centered: bool) -> str:
    meter_value = _control_bar_value(value, centered=centered)
    if centered:
        left = min(50.0, meter_value * 100.0)
        width = abs(meter_value - 0.5) * 100.0
    else:
        left = 0.0
        width = meter_value * 100.0
    center = '<div class="fd-control-meter__center"></div>' if centered else ""
    return (
        '<div class="fd-control-meter">'
        f"{center}"
        f'<div class="fd-control-meter__fill" style="left: {left:.2f}%; width: {width:.2f}%;"></div>'
        "</div>"
    )


def _maneuver_badge_html(label: str) -> str:
    return (
        '<div class="fd-maneuver-badge">'
        "<small>MANEUVER</small>"
        f"<strong>{label.upper()}</strong>"
        "</div>"
    )


def _telemetry_card_html(*, label: str, value: float, unit: str) -> str:
    return (
        '<div class="fd-telemetry-card">'
        f'<div class="fd-telemetry-card__value">{value:.1f}</div>'
        f'<div class="fd-telemetry-card__unit">{unit}</div>'
        f'<div class="fd-telemetry-card__label">{label}</div>'
        "</div>"
    )


def _maneuver_metadata(config: SimulationConfig, maneuver_key: str) -> list[dict[str, str]]:
    from .flight_deck import maneuver_by_key

    maneuver = maneuver_by_key(config, maneuver_key)
    rows = [{"field": "mode", "value": maneuver.label}]
    for index, command in enumerate(maneuver.commands, start=1):
        amplitude = command.amplitude
        unit = ""
        if command.surface.strip().lower() in {"elevator", "de", "aileron", "da", "rudder", "ruddervator", "dr"}:
            amplitude = float(np.rad2deg(amplitude))
            unit = " deg"
        elif command.surface.strip().lower() == "throttle":
            amplitude = 100.0 * amplitude
            unit = "%"
        rows.append(
            {
                "field": f"cmd {index}",
                "value": (
                    f"{command.surface} {command.kind} {amplitude:+.2f}{unit}, "
                    f"t={command.start_time:.1f}s, dur={command.duration:.1f}s"
                ),
            }
        )
    return rows


def _render_maneuver_metadata(ui, container, rows: list[dict[str, str]]) -> None:
    container.clear()
    with container:
        for row in rows:
            with ui.column().classes("w-full gap-0 fd-card px-3 py-2"):
                ui.label(row["field"].upper()).classes("text-xs fd-muted")
                ui.label(row["value"]).classes("text-sm text-sky-100")


def _render_preset_list(ui, container, options: dict[str, str], active_key: str) -> None:
    container.clear()
    with container:
        ui.label("PRESETS").classes("text-xs font-semibold text-sky-200")
        for key, label in options.items():
            active_class = " fd-preset--active" if key == active_key else ""
            ui.label(label.upper()).classes(f"w-full fd-preset{active_class} px-2 py-1 text-xs font-mono")


def _timeline_html(elapsed: float, total: float) -> str:
    progress = 0.0 if total <= 0.0 else max(0.0, min(1.0, elapsed / total))
    return f'<div class="fd-timeline"><div class="fd-timeline__fill" style="width: {100.0 * progress:.1f}%;"></div></div>'


def _render_live_plot_cards(ui, container, history: list[dict[str, float]]) -> None:
    container.clear()
    with container:
        with ui.element("div").classes("fd-plot-grid w-full"):
            for figure in flight_deck_live_figures(history).values():
                with ui.element("div").classes("fd-plot-card"):
                    ui.plotly(figure).classes("w-full")


def _register_flight_deck_page(ui) -> None:
    @ui.page("/")
    def flight_deck_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = APP_STATE.scenario_source
        if selected_name not in scenario_options:
            selected_name = next(iter(scenario_options), None)

        try:
            session = _active_flight_deck_session(selected_name)
        except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
            with ui.column().classes("w-full max-w-5xl mx-auto gap-4 p-4"):
                _page_intro(ui, "Flight Deck", "Realtime aircraft simulation console.")
                ui.label(f"Unable to initialize Flight Deck: {exc}").classes("text-red-600")
            return

        telemetry_cards: dict[str, Any] = {}
        control_labels: dict[str, Any] = {}
        control_bars: dict[str, Any] = {}
        control_centered: dict[str, bool] = {}
        plot_container = None
        attitude = None
        last_plot_time = -1.0

        with ui.column().classes("fd-console w-full gap-3 p-3"):
            with ui.row().classes("w-full items-center justify-between fd-panel px-4 py-3"):
                with ui.row().classes("items-center gap-3"):
                    ui.label("FLIGHT DECK").classes("text-xl font-semibold text-sky-100")
                    maneuver_badge = ui.html(_maneuver_badge_html(session.maneuver.label))
                    scenario_badge = ui.label(session.scenario_source).classes("text-sm fd-muted")
                with ui.row().classes("items-center gap-3"):
                    flight_condition_label = ui.label("IAS 0.0 m/s | ALT 0.0 m").classes(
                        "text-sm font-mono text-sky-100"
                    )
                    status_badge = ui.label(session.status).classes(
                        f"text-sm font-semibold border rounded px-3 py-1 {_status_classes(session.status)}"
                    )
                    time_label = ui.label("T+0.00 s").classes("text-sm font-mono text-sky-100")
                    speed_label = ui.label("1.0x").classes("text-sm font-mono text-sky-100")

            with ui.grid(columns="220px minmax(620px, 1.4fr) 250px").classes("w-full gap-3"):
                with ui.column().classes("fd-panel gap-2 p-3"):
                    ui.label("SETUP").classes("text-xs font-semibold text-sky-200")
                    scenario_select = ui.select(
                        options=list(scenario_options),
                        value=selected_name,
                        label="Scenario",
                    ).props("dark filled dense").classes("w-full")
                    maneuver_select = ui.select(
                        options=maneuver_options(session.base_config),
                        value=session.maneuver_key,
                        label="Maneuver",
                    ).props("dark filled dense").classes("w-full")
                    preset_container = ui.column().classes("w-full gap-1")
                    _render_preset_list(
                        ui,
                        preset_container,
                        maneuver_options(session.base_config),
                        session.maneuver_key,
                    )
                    metadata_container = ui.column().classes("w-full gap-2")
                    _render_maneuver_metadata(
                        ui,
                        metadata_container,
                        _maneuver_metadata(session.base_config, session.maneuver_key),
                    )
                    ui.label("Trim").classes("text-xs font-semibold text-sky-200")
                    trim_label = ui.label(
                        f"{session.base_config.trim_target_speed:g} m/s at "
                        f"{session.base_config.trim_altitude:g} m"
                    ).classes("text-sm text-slate-300")

                with ui.column().classes("fd-panel gap-3 p-3"):
                    with ui.row().classes("w-full items-center justify-between"):
                        ui.label("PRIMARY FLIGHT DISPLAY").classes("text-xs font-semibold text-sky-200")
                        status_message = ui.label(session.message).classes("text-sm text-slate-300")
                    values = session.telemetry()
                    attitude = ui.html(
                        attitude_html(
                            pitch_deg=values.get("pitch_deg", 0.0),
                            roll_deg=values.get("roll_deg", 0.0),
                        )
                    ).classes("w-full")

                with ui.column().classes("fd-panel gap-2 p-3"):
                    ui.label("TELEMETRY").classes("text-xs font-semibold text-sky-200")
                    telemetry_fields = [
                        ("airspeed", "IAS", "m/s"),
                        ("altitude", "ALT", "m"),
                        ("pitch_deg", "PITCH", "deg"),
                        ("roll_deg", "ROLL", "deg"),
                        ("heading_deg", "HDG", "deg"),
                        ("alpha_deg", "ALPHA", "deg"),
                        ("beta_deg", "BETA", "deg"),
                    ]
                    for key, label, unit in telemetry_fields:
                        telemetry_cards[key] = ui.html(
                            _telemetry_card_html(label=label, value=0.0, unit=unit)
                        ).classes("w-full")

                    ui.label("CONTROLS").classes("text-xs font-semibold text-sky-200 mt-2")
                    for key, label, centered in [
                        ("elevator_deg", "ELEV", True),
                        ("aileron_deg", "AIL", True),
                        ("rudder_deg", "RUD", True),
                        ("throttle_pct", "THR", False),
                    ]:
                        with ui.column().classes("w-full gap-1"):
                            with ui.row().classes("w-full justify-between"):
                                ui.label(label).classes("text-xs fd-muted")
                                control_labels[key] = ui.label("0.0").classes("text-xs font-mono text-sky-100")
                            control_centered[key] = centered
                            control_bars[key] = ui.html(_control_bar_html(0.0, centered=centered)).classes("w-full")

            with ui.row().classes("w-full items-center justify-between fd-panel fd-transport px-4 py-3"):
                with ui.row().classes("items-center gap-2"):
                    reset_button = ui.button("RESET")
                    step_button = ui.button("STEP")
                    play_button = ui.button("PLAY").props("color=accent")
                with ui.row().classes("items-center gap-3"):
                    ui.label("SPEED").classes("text-xs fd-muted")
                    speed_select = ui.select(
                        options={0.5: "0.5x", 1.0: "1x", 2.0: "2x", 5.0: "5x"},
                        value=session.speed_multiplier,
                    ).props("dark filled dense").classes("w-28")
                    timeline_meter = ui.html(_timeline_html(0.0, session.config.duration))
                timeline_label = ui.label("0.00 / 0.00 s").classes("text-sm font-mono text-sky-100")

            with ui.column().classes("w-full fd-panel p-2"):
                ui.label("LIVE RESPONSE STRIPS").classes("text-xs font-semibold text-sky-200")
                plot_container = ui.column().classes("w-full")
                _render_live_plot_cards(ui, plot_container, session.history)

        def sync_latest_result() -> None:
            result = session.latest_result()
            if result is not None:
                APP_STATE.set_simulation_result(result)

        def render(*, force_plot: bool = False) -> None:
            nonlocal last_plot_time
            values = session.telemetry()
            controls = session.controls()
            scenario_badge.text = session.scenario_source
            maneuver_badge.content = _maneuver_badge_html(session.maneuver.label)
            maneuver_badge.update()
            status_badge.text = session.status
            status_badge.classes(replace=f"text-sm font-semibold border rounded px-3 py-1 {_status_classes(session.status)}")
            time_value = values.get("time", 0.0)
            time_label.text = f"T+{time_value:.2f} s"
            speed_label.text = f"{session.speed_multiplier:g}x"
            timeline_label.text = f"{time_value:.2f} / {session.config.duration:.2f} s"
            timeline_meter.content = _timeline_html(time_value, session.config.duration)
            timeline_meter.update()
            flight_condition_label.text = (
                f"IAS {values.get('airspeed', 0.0):.1f} m/s | ALT {values.get('altitude', 0.0):.0f} m"
            )
            play_button.text = "PAUSE" if session.status == RUNNING else "PLAY"
            status_message.text = session.message
            trim_label.text = (
                f"{session.base_config.trim_target_speed:g} m/s at "
                f"{session.base_config.trim_altitude:g} m"
            )

            units = {
                "airspeed": "m/s",
                "altitude": "m",
                "pitch_deg": "deg",
                "roll_deg": "deg",
                "heading_deg": "deg",
                "alpha_deg": "deg",
                "beta_deg": "deg",
            }
            labels = {
                "airspeed": "IAS",
                "altitude": "ALT",
                "pitch_deg": "PITCH",
                "roll_deg": "ROLL",
                "heading_deg": "HDG",
                "alpha_deg": "ALPHA",
                "beta_deg": "BETA",
            }
            for key, card in telemetry_cards.items():
                card.content = _telemetry_card_html(
                    label=labels[key],
                    value=values.get(key, 0.0),
                    unit=units[key],
                )
                card.update()
            for key, label in control_labels.items():
                label.text = f"{controls.get(key, 0.0):+.1f}"
            for key, bar in control_bars.items():
                bar.content = _control_bar_html(
                    controls.get(key, 0.0),
                    centered=control_centered[key],
                )
                bar.update()
            if attitude is not None:
                attitude.content = attitude_html(
                    pitch_deg=values.get("pitch_deg", 0.0),
                    roll_deg=values.get("roll_deg", 0.0),
                )
                attitude.update()

            if force_plot or time_value - last_plot_time >= 0.25:
                _render_live_plot_cards(ui, plot_container, session.history)
                last_plot_time = time_value
            sync_latest_result()

        def load_scenario() -> None:
            try:
                if not scenario_select.value:
                    raise ValueError("No example scenario is selected.")
                config = _load_example_config(str(scenario_select.value))
                session.set_scenario(str(scenario_select.value), config)
                maneuver_select.options = maneuver_options(config)
                maneuver_select.value = session.maneuver_key
                maneuver_select.update()
                _render_preset_list(ui, preset_container, maneuver_options(config), session.maneuver_key)
                _render_maneuver_metadata(ui, metadata_container, _maneuver_metadata(session.base_config, session.maneuver_key))
                render(force_plot=True)
            except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                session.status = ERROR
                session.message = f"Failed to load scenario: {exc}"
                ui.notify(session.message, type="negative")
                render(force_plot=True)

        def select_maneuver() -> None:
            session.set_maneuver(str(maneuver_select.value))
            _render_preset_list(ui, preset_container, maneuver_options(session.base_config), session.maneuver_key)
            _render_maneuver_metadata(ui, metadata_container, _maneuver_metadata(session.base_config, session.maneuver_key))
            render(force_plot=True)

        def reset_clicked() -> None:
            session.reset()
            render(force_plot=True)

        def step_clicked() -> None:
            session.step_once()
            render(force_plot=True)

        def play_clicked() -> None:
            session.toggle_playback()
            render(force_plot=True)

        def speed_changed() -> None:
            session.speed_multiplier = float(speed_select.value)
            render()

        def tick() -> None:
            if session.status == RUNNING:
                session.advance(0.1)
                render()

        scenario_select.on("update:model-value", lambda _: load_scenario())
        maneuver_select.on("update:model-value", lambda _: select_maneuver())
        speed_select.on("update:model-value", lambda _: speed_changed())
        reset_button.on_click(reset_clicked)
        step_button.on_click(step_clicked)
        play_button.on_click(play_clicked)
        ui.timer(0.1, tick)
        render(force_plot=True)


def _register_analyze_page(ui) -> None:
    @ui.page("/analyze")
    def analyze_page() -> None:
        _render_page_shell(ui)
        with ui.column().classes("w-full max-w-6xl mx-auto gap-4 p-4"):
            _page_intro(ui, "Analyze", "Review the latest run and open the advanced analysis tools.")

            with ui.grid(columns=3).classes("w-full gap-3"):
                with ui.card().classes("gap-2 p-4 shadow-sm"):
                    ui.label("Simulation Review").classes("text-lg font-semibold")
                    ui.label("Batch run, plots, and exports.").classes("text-slate-600")
                    ui.link("Open Simulation", "/simulation").classes("text-blue-700 no-underline")
                with ui.card().classes("gap-2 p-4 shadow-sm"):
                    ui.label("Stability").classes("text-lg font-semibold")
                    ui.label("Trim, linearization, modes, derivatives.").classes("text-slate-600")
                    ui.link("Open Stability", "/stability").classes("text-blue-700 no-underline")
                with ui.card().classes("gap-2 p-4 shadow-sm"):
                    ui.label("Advanced CG Sweep").classes("text-lg font-semibold")
                    ui.label("Runtime CG stability sweep.").classes("text-slate-600")
                    ui.link("Open CG Sweep", "/cg-sweep").classes("text-blue-700 no-underline")

            result = APP_STATE.simulation_result
            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Latest Flight Deck Run")
                if result is None:
                    _empty_state(ui, "Press Play on the Flight Deck to capture a run.")
                    return
                summary = simulation_result_summary(result)
                ui.label(
                    f"{summary['name']} | samples={summary['samples']} | "
                    f"t={summary['start_time']:.2f}-{summary['end_time']:.2f}s | "
                    f"success={summary['success']}"
                ).classes("text-slate-600")
                ui.plotly(simulation_time_history_figure(result)).classes("w-full")


def _register_scenario_page(ui) -> None:
    @ui.page("/scenario-library")
    def scenario_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = next(iter(scenario_options), None)

        with ui.column().classes("w-full max-w-5xl mx-auto gap-4 p-4"):
            _page_intro(ui, "Scenario", "Choose a YAML scenario for the simulation and stability pages.")

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Load Scenario", "Start with the bundled example or upload a YAML file.")
                with ui.row().classes("items-end gap-3"):
                    selected = ui.select(
                        options=list(scenario_options),
                        value=selected_name,
                        label="Example scenario",
                    ).classes("min-w-80")
                    ui.button("Load Scenario", on_click=lambda: load_selected()).props("color=primary")

                ui.upload(
                    label="Upload YAML scenario",
                    auto_upload=True,
                    on_upload=lambda event: load_uploaded(event),
                ).props("accept=.yaml,.yml").classes("w-full")

            summary_columns = [
                {"name": "field", "label": "Field", "field": "field", "align": "left"},
                {"name": "value", "label": "Value", "field": "value", "align": "left"},
            ]
            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Scenario Summary", "Essential fields from the active scenario.")
                status = ui.label("Load a scenario first.").classes("text-slate-600")
                empty = _empty_state(ui, "Load a scenario first.")
                essential_table = ui.table(columns=summary_columns, rows=[], row_key="field").classes("w-full")
                essential_table.visible = False

                with ui.expansion("Raw / config details").classes("w-full"):
                    raw_table = ui.table(columns=summary_columns, rows=[], row_key="field").classes("w-full")

            def show_summary(source: str, summary: dict[str, Any]) -> None:
                essential_table.rows = _essential_scenario_rows(summary)
                essential_table.visible = True
                essential_table.update()
                raw_table.rows = _summary_rows(summary)
                raw_table.update()
                empty.visible = False
                status.text = f"Loaded {source}"

            def show_error(message: str) -> None:
                essential_table.rows = []
                essential_table.visible = False
                essential_table.update()
                raw_table.rows = []
                raw_table.update()
                empty.visible = True
                status.text = message
                ui.notify(message, type="negative")

            def load_selected() -> None:
                if not selected.value:
                    show_error("No example scenario is selected.")
                    return
                try:
                    config = _load_example_config(selected.value)
                    show_summary(selected.value, scenario_summary(config))
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    show_error(f"Failed to load scenario: {exc}")

            def load_uploaded(upload_event) -> None:
                try:
                    source, summary = _load_uploaded_config(upload_event)
                    show_summary(source, summary)
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    show_error(f"Failed to load uploaded scenario: {exc}")

            if APP_STATE.config is not None:
                show_summary(APP_STATE.scenario_source or APP_STATE.config.name, scenario_summary(APP_STATE.config))
            elif selected_name is None:
                status.text = f"No example scenarios found in {EXAMPLE_SCENARIOS_DIR}"
            else:
                raw_table.rows = []
                raw_table.update()


def _register_simulation_page(ui) -> None:
    @ui.page("/simulation")
    def simulation_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = APP_STATE.scenario_source
        if selected_name not in scenario_options:
            selected_name = next(iter(scenario_options), None)

        active_config = APP_STATE.config
        if active_config is None and selected_name is not None:
            active_config = load_scenario_config(scenario_options[selected_name])

        with ui.column().classes("w-full max-w-6xl mx-auto gap-4 p-4"):
            _page_intro(ui, "Simulation", "Run the loaded scenario through the shared nonlinear simulation backend.")

            if active_config is None:
                _empty_state(ui, "Load a scenario first.")
                return

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Run Simulation", "Use the active scenario, or load one of the bundled examples.")
                status = ui.label(
                    f"Loaded scenario: {APP_STATE.scenario_source or active_config.name}"
                    if APP_STATE.config is not None
                    else "Load a scenario first, or use the bundled example below."
                ).classes("text-slate-600")
                with ui.row().classes("items-end gap-3"):
                    selected = ui.select(
                        options=list(scenario_options),
                        value=selected_name,
                        label="Example scenario",
                    ).classes("min-w-80")
                    ui.button("Load Scenario", on_click=lambda: load_selected_for_form())
                    ui.button("Run Simulation", on_click=lambda: run_clicked()).props("color=primary")

                with ui.expansion("Advanced settings").classes("w-full"):
                    with ui.grid(columns=3).classes("w-full gap-3"):
                        name_input = ui.input("Scenario name", value=active_config.name)
                        duration_input = ui.number("Duration [s]", value=active_config.duration, min=0.001, step=0.5)
                        dt_input = ui.input(
                            "Time step dt [s]",
                            value="" if active_config.dt is None else str(active_config.dt),
                        )
                        n_points_input = ui.input(
                            "Sample count",
                            value="" if active_config.n_points is None else str(active_config.n_points),
                        )
                        trim_toggle = ui.switch("Start from trim", value=active_config.start_from_trim)
                        trim_speed_input = ui.number(
                            "Trim airspeed [m/s]",
                            value=active_config.trim_target_speed,
                            min=0.001,
                        )
                        trim_altitude_input = ui.number("Trim altitude [m]", value=active_config.trim_altitude)

                    command_columns = [
                        {"name": "surface", "label": "Surface", "field": "surface", "align": "left"},
                        {"name": "kind", "label": "Type", "field": "kind", "align": "left"},
                        {"name": "amplitude", "label": "Amplitude", "field": "amplitude", "align": "right"},
                        {"name": "start_time", "label": "Start [s]", "field": "start_time", "align": "right"},
                        {"name": "duration", "label": "Duration [s]", "field": "duration", "align": "right"},
                    ]
                    commands_table = ui.table(
                        columns=command_columns,
                        rows=_command_rows(active_config.input_commands),
                        row_key="surface",
                    ).classes("w-full")

            summary_columns = [
                {"name": "field", "label": "Field", "field": "field", "align": "left"},
                {"name": "value", "label": "Value", "field": "value", "align": "left"},
            ]
            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Result Summary", "Simulation status, sample count, final time, and trim result.")
                summary_empty = _empty_state(ui, "Run simulation to see results.")
                summary_table = ui.table(columns=summary_columns, rows=[], row_key="field").classes("w-full")
                summary_table.visible = False

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Plots and Exports", "Inspect time histories and download simulation outputs.")
                plot_empty = _empty_state(ui, "Run simulation to see plots.")
                plot_container = ui.column().classes("w-full")
                export_row = ui.row().classes("gap-3")
                export_row.visible = False

            def load_selected_for_form() -> None:
                try:
                    if not selected.value:
                        raise ValueError("No example scenario is selected.")
                    config = _load_example_config(selected.value)
                    name_input.value = config.name
                    duration_input.value = config.duration
                    dt_input.value = "" if config.dt is None else config.dt
                    n_points_input.value = "" if config.n_points is None else config.n_points
                    trim_toggle.value = config.start_from_trim
                    trim_speed_input.value = config.trim_target_speed
                    trim_altitude_input.value = config.trim_altitude
                    commands_table.rows = _command_rows(config.input_commands)
                    commands_table.update()
                    summary_table.rows = []
                    summary_table.visible = False
                    summary_table.update()
                    summary_empty.visible = True
                    plot_container.clear()
                    plot_empty.visible = True
                    render_exports()
                    status.text = f"Loaded scenario: {selected.value}"
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    status.text = f"Failed to load scenario: {exc}"
                    ui.notify(status.text, type="negative")

            def config_for_form() -> SimulationConfig:
                if APP_STATE.config is not None:
                    return APP_STATE.config
                if not selected.value:
                    raise ValueError("Load a scenario first.")
                return _load_example_config(selected.value)

            def form_config() -> SimulationConfig:
                base_config = config_for_form()
                raw_dt = str(dt_input.value).strip()
                raw_n_points = str(n_points_input.value).strip()
                dt = float(raw_dt) if raw_dt else None
                n_points = int(raw_n_points) if raw_n_points else None
                return build_simulation_config_from_form(
                    base_config,
                    name=str(name_input.value).strip() or base_config.name,
                    duration=float(duration_input.value),
                    dt=dt,
                    n_points=n_points,
                    start_from_trim=bool(trim_toggle.value),
                    trim_target_speed=float(trim_speed_input.value),
                    trim_altitude=float(trim_altitude_input.value),
                )

            def render_exports() -> None:
                export_row.clear()
                result = APP_STATE.simulation_result
                if result is None:
                    export_row.visible = False
                    return
                export_row.visible = True
                with export_row:
                    ui.button(
                        "Export CSV",
                        on_click=lambda: ui.download(
                            export_simulation_csv_text(result),
                            filename=f"{result.config.name}_simulation.csv",
                        ),
                    )
                    ui.button(
                        "Export JSON",
                        on_click=lambda: ui.download(
                            export_simulation_json_text(result),
                            filename=f"{result.config.name}_simulation.json",
                        ),
                    )

            def run_clicked() -> None:
                try:
                    config = form_config()
                    APP_STATE.set_scenario(APP_STATE.scenario_source or config.name, config)
                    status.text = "Running simulation..."
                    result = run_simulation_for_gui(config)
                    APP_STATE.set_simulation_result(result)
                    summary = simulation_result_summary(result)
                    summary_table.rows = _simulation_summary_rows(summary)
                    summary_table.visible = True
                    summary_table.update()
                    summary_empty.visible = False
                    plot_container.clear()
                    with plot_container:
                        ui.plotly(simulation_time_history_figure(result)).classes("w-full")
                    plot_empty.visible = False
                    render_exports()
                    status.text = f"Simulation complete: success={result.success}"
                    if not result.success:
                        ui.notify(result.message, type="warning")
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    message = f"Simulation failed: {exc}"
                    APP_STATE.set_simulation_error(message)
                    summary_table.rows = []
                    summary_table.visible = False
                    summary_table.update()
                    summary_empty.visible = True
                    plot_container.clear()
                    plot_empty.visible = True
                    render_exports()
                    status.text = message
                    ui.notify(message, type="negative")


def _register_stability_page(ui) -> None:
    @ui.page("/stability")
    def stability_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = APP_STATE.scenario_source
        if selected_name not in scenario_options:
            selected_name = next(iter(scenario_options), None)

        active_config = APP_STATE.config
        if active_config is None and selected_name is not None:
            active_config = load_scenario_config(scenario_options[selected_name])

        with ui.column().classes("w-full max-w-6xl mx-auto gap-4 p-4"):
            _page_intro(ui, "Stability", "Analyze the active scenario with the shared trim, linearization, and modal tools.")

            if active_config is None:
                _empty_state(ui, "Load a scenario first.")
                return

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Run Stability Analysis", "Use the active scenario, or load one of the bundled examples.")
                status = ui.label(
                    f"Loaded scenario: {APP_STATE.scenario_source or active_config.name}"
                    if APP_STATE.config is not None
                    else "Load a scenario first, or use the bundled example below."
                ).classes("text-slate-600")
                with ui.row().classes("items-end gap-3"):
                    selected = ui.select(
                        options=list(scenario_options),
                        value=selected_name,
                        label="Example scenario",
                    ).classes("min-w-80")
                    ui.button("Load Scenario", on_click=lambda: load_selected_for_analysis())
                    ui.button("Run Stability Analysis", on_click=lambda: run_clicked()).props("color=primary")

                with ui.expansion("Advanced settings").classes("w-full"):
                    with ui.grid(columns=2).classes("w-full gap-3"):
                        state_perturbation_input = ui.number(
                            "State perturbation",
                            value=1e-5,
                            min=1e-12,
                        )
                        control_perturbation_input = ui.number(
                            "Control perturbation",
                            value=1e-5,
                            min=1e-12,
                        )

            summary_columns = [
                {"name": "field", "label": "Field", "field": "field", "align": "left"},
                {"name": "value", "label": "Value", "field": "value", "align": "left"},
            ]
            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Analysis Summary", "Trim condition and full/reduced matrix dimensions.")
                summary_empty = _empty_state(ui, "Run stability analysis to see summary.")
                summary_table = ui.table(columns=summary_columns, rows=[], row_key="field").classes("w-full")
                summary_table.visible = False

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Modes, Derivatives, and Exports", "Inspect poles, modal metrics, static derivatives, and reports.")
                analysis_empty = _empty_state(ui, "Run stability analysis to see modes.")

                mode_columns = [
                    {"name": "subsystem", "label": "Subsystem", "field": "subsystem", "align": "left"},
                    {"name": "mode", "label": "Mode", "field": "mode", "align": "left"},
                    {"name": "eigenvalue_real", "label": "Real", "field": "eigenvalue_real", "align": "right"},
                    {"name": "eigenvalue_imag", "label": "Imaginary", "field": "eigenvalue_imag", "align": "right"},
                    {
                        "name": "natural_frequency_rad_s",
                        "label": "wn [rad/s]",
                        "field": "natural_frequency_rad_s",
                        "align": "right",
                    },
                    {"name": "damping_ratio", "label": "Damping", "field": "damping_ratio", "align": "right"},
                    {
                        "name": "oscillation_period_s",
                        "label": "Period [s]",
                        "field": "oscillation_period_s",
                        "align": "right",
                    },
                    {"name": "time_to_half_s", "label": "Half [s]", "field": "time_to_half_s", "align": "right"},
                    {"name": "time_to_double_s", "label": "Double [s]", "field": "time_to_double_s", "align": "right"},
                    {"name": "classification", "label": "Class", "field": "classification", "align": "left"},
                ]
                modes_table = ui.table(columns=mode_columns, rows=[], row_key="eigenvalue_label").classes("w-full")
                modes_table.visible = False

                derivative_columns = [
                    {"name": "name", "label": "Derivative", "field": "name", "align": "left"},
                    {"name": "value", "label": "Value", "field": "value", "align": "right"},
                    {"name": "alpha_deg", "label": "Alpha [deg]", "field": "alpha_deg", "align": "right"},
                    {"name": "beta_deg", "label": "Beta [deg]", "field": "beta_deg", "align": "right"},
                    {"name": "source", "label": "Source", "field": "source", "align": "left"},
                ]
                derivatives_table = ui.table(columns=derivative_columns, rows=[], row_key="name").classes("w-full")
                derivatives_table.visible = False
                plot_container = ui.column().classes("w-full")
                export_row = ui.row().classes("gap-3")
                export_row.visible = False

            def load_selected_for_analysis() -> None:
                try:
                    if not selected.value:
                        raise ValueError("No example scenario is selected.")
                    config = _load_example_config(selected.value)
                    summary_table.rows = []
                    summary_table.visible = False
                    summary_table.update()
                    summary_empty.visible = True
                    modes_table.rows = []
                    modes_table.visible = False
                    modes_table.update()
                    derivatives_table.rows = []
                    derivatives_table.visible = False
                    derivatives_table.update()
                    plot_container.clear()
                    analysis_empty.visible = True
                    render_exports()
                    status.text = f"Loaded scenario: {selected.value}"
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    status.text = f"Failed to load scenario: {exc}"
                    ui.notify(status.text, type="negative")

            def render_exports() -> None:
                export_row.clear()
                result = APP_STATE.stability_result
                if result is None:
                    export_row.visible = False
                    return
                export_row.visible = True
                with export_row:
                    ui.button(
                        "Export JSON",
                        on_click=lambda: ui.download(
                            export_stability_json_text(result),
                            filename=f"{result.config.name}_stability.json",
                        ),
                    )
                    ui.button(
                        "Export Modes CSV",
                        on_click=lambda: ui.download(
                            export_stability_modes_csv_text(result),
                            filename=f"{result.config.name}_stability_modes.csv",
                        ),
                    )
                    ui.button(
                        "Export Markdown",
                        on_click=lambda: ui.download(
                            export_stability_markdown_text(result),
                            filename=f"{result.config.name}_stability.md",
                        ),
                    )

            def run_clicked() -> None:
                try:
                    config = _active_or_example_config(selected.value if selected.value else None)
                    linearization_config = build_linearization_config(
                        state_perturbation=float(state_perturbation_input.value),
                        control_perturbation=float(control_perturbation_input.value),
                    )
                    status.text = "Running stability analysis..."
                    result = run_stability_for_gui(config, linearization_config=linearization_config)
                    APP_STATE.set_stability_result(result)
                    summary_table.rows = stability_result_summary_rows(result)
                    summary_table.visible = True
                    summary_table.update()
                    summary_empty.visible = False
                    modes_table.rows = stability_mode_table_rows(result)
                    modes_table.visible = True
                    modes_table.update()
                    derivatives_table.rows = static_derivative_table_rows(result)
                    derivatives_table.visible = True
                    derivatives_table.update()
                    plot_container.clear()
                    with plot_container:
                        ui.plotly(stability_eigenvalue_figure(result)).classes("w-full")
                    analysis_empty.visible = False
                    render_exports()
                    status.text = "Stability analysis complete: success=True"
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    message = f"Stability analysis failed: {exc}"
                    APP_STATE.set_stability_error(message)
                    summary_table.rows = [{"field": "success", "value": "False"}, {"field": "error", "value": str(exc)}]
                    summary_table.visible = True
                    summary_table.update()
                    summary_empty.visible = False
                    modes_table.rows = []
                    modes_table.visible = False
                    modes_table.update()
                    derivatives_table.rows = []
                    derivatives_table.visible = False
                    derivatives_table.update()
                    plot_container.clear()
                    analysis_empty.visible = True
                    render_exports()
                    status.text = message
                    ui.notify(message, type="negative")


def _register_cg_sweep_page(ui) -> None:
    @ui.page("/cg-sweep")
    def cg_sweep_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = APP_STATE.scenario_source
        if selected_name not in scenario_options:
            selected_name = next(iter(scenario_options), None)

        active_config = APP_STATE.config
        if active_config is None and selected_name is not None:
            active_config = load_scenario_config(scenario_options[selected_name])

        with ui.column().classes("w-full max-w-6xl mx-auto gap-4 p-4"):
            _page_intro(ui, "CG Sweep", "Run runtime center-of-gravity stability sweeps through the shared backend.")

            if active_config is None:
                _empty_state(ui, "Load a scenario first.")
                return

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Run CG Sweep", "Enter CG x positions and run the existing stability sweep.")
                ui.label(
                    "This sweep varies runtime `cg_body_xyz` with the nominal ADB. "
                    "It does not use separate per-CG aerodynamic databases."
                ).classes("text-slate-600")
                status = ui.label(
                    f"Loaded scenario: {APP_STATE.scenario_source or active_config.name}"
                    if APP_STATE.config is not None
                    else "Load a scenario first, or use the bundled example below."
                ).classes("text-slate-600")
                with ui.row().classes("items-end gap-3"):
                    selected = ui.select(
                        options=list(scenario_options),
                        value=selected_name,
                        label="Example scenario",
                    ).classes("min-w-80")
                    ui.button("Load Scenario", on_click=lambda: load_selected_for_sweep())
                    ui.button("Run CG Sweep", on_click=lambda: run_clicked()).props("color=primary")

                cg_x_input = ui.input(
                    "CG x values [m]",
                    value="-0.15, -0.10, -0.05, 0.0",
                    placeholder="-0.15, -0.10, -0.05, 0.0",
                ).classes("w-full")

                with ui.expansion("Advanced settings").classes("w-full"):
                    with ui.grid(columns=2).classes("w-full gap-3"):
                        state_perturbation_input = ui.number(
                            "State perturbation",
                            value=1e-5,
                            min=1e-12,
                        )
                        control_perturbation_input = ui.number(
                            "Control perturbation",
                            value=1e-5,
                            min=1e-12,
                        )

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Sweep Results", "Per-point success, mode counts, and most unstable pole.")
                results_empty = _empty_state(ui, "Run CG sweep to see results.")
                result_columns = [
                    {"name": "case", "label": "Case", "field": "case", "align": "left"},
                    {"name": "cg_x", "label": "CG x [m]", "field": "cg_x", "align": "right"},
                    {"name": "success", "label": "Success", "field": "success", "align": "left"},
                    {"name": "total_modes", "label": "Modes", "field": "total_modes", "align": "right"},
                    {
                        "name": "most_unstable_real",
                        "label": "Max Real",
                        "field": "most_unstable_real",
                        "align": "right",
                    },
                    {"name": "error", "label": "Error", "field": "error", "align": "left"},
                ]
                results_table = ui.table(columns=result_columns, rows=[], row_key="case").classes("w-full")
                results_table.visible = False

            with ui.card().classes("w-full gap-3 p-4 shadow-sm"):
                _section_header(ui, "Trend and Exports", "Inspect CG sensitivity and download sweep outputs.")
                trend_empty = _empty_state(ui, "Run CG sweep to see trend.")
                trend_container = ui.column().classes("w-full")
                export_row = ui.row().classes("gap-3")
                export_row.visible = False

            def load_selected_for_sweep() -> None:
                try:
                    if not selected.value:
                        raise ValueError("No example scenario is selected.")
                    _load_example_config(selected.value)
                    results_table.rows = []
                    results_table.visible = False
                    results_table.update()
                    results_empty.visible = True
                    trend_container.clear()
                    trend_empty.visible = True
                    render_exports()
                    status.text = f"Loaded scenario: {selected.value}"
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    status.text = f"Failed to load scenario: {exc}"
                    ui.notify(status.text, type="negative")

            def render_exports() -> None:
                export_row.clear()
                result = APP_STATE.cg_sweep_result
                if result is None:
                    export_row.visible = False
                    return
                export_row.visible = True
                with export_row:
                    ui.button(
                        "Export JSON",
                        on_click=lambda: ui.download(
                            export_stability_json_text(result),
                            filename=f"{result.base_config.name}_cg_sweep.json",
                        ),
                    )
                    ui.button(
                        "Export CSV",
                        on_click=lambda: ui.download(
                            export_stability_modes_csv_text(result),
                            filename=f"{result.base_config.name}_cg_sweep_modes.csv",
                        ),
                    )

            def run_clicked() -> None:
                try:
                    config = _active_or_example_config(selected.value if selected.value else None)
                    x_values = parse_cg_x_values(str(cg_x_input.value))
                    cg_positions = cg_positions_from_x_values(config, x_values)
                    linearization_config = build_linearization_config(
                        state_perturbation=float(state_perturbation_input.value),
                        control_perturbation=float(control_perturbation_input.value),
                    )
                    status.text = "Running CG sweep..."
                    result = run_cg_sweep_for_gui(
                        config,
                        cg_positions,
                        linearization_config=linearization_config,
                    )
                    APP_STATE.set_cg_sweep_result(result)
                    results_table.rows = cg_sweep_summary_rows(result)
                    results_table.visible = True
                    results_table.update()
                    results_empty.visible = False
                    trend_container.clear()
                    with trend_container:
                        ui.plotly(cg_sweep_trend_figure(result)).classes("w-full")
                    trend_empty.visible = False
                    render_exports()
                    status.text = (
                        f"CG sweep complete: successful={len(result.successful_points)} "
                        f"total={len(result.points)}"
                    )
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    message = f"CG sweep failed: {exc}"
                    APP_STATE.set_cg_sweep_error(message)
                    results_table.rows = []
                    results_table.visible = False
                    results_table.update()
                    results_empty.visible = True
                    trend_container.clear()
                    trend_empty.visible = True
                    render_exports()
                    status.text = message
                    ui.notify(message, type="negative")


def create_app() -> None:
    """Register GUI pages with NiceGUI."""
    global _pages_registered
    if _pages_registered:
        return
    ui = _require_nicegui()
    _register_flight_deck_page(ui)
    _register_analyze_page(ui)
    _register_scenario_page(ui)
    _register_simulation_page(ui)
    _register_stability_page(ui)
    _register_cg_sweep_page(ui)
    _pages_registered = True


def run_gui(*, host: str = DEFAULT_HOST, port: int = DEFAULT_PORT, reload: bool = False, show: bool = True) -> None:
    """Run the NiceGUI app."""
    ui = _require_nicegui()
    create_app()
    ui.run(title="Spearhead Flight Dynamics Model", host=host, port=port, reload=reload, show=show)


def build_parser() -> argparse.ArgumentParser:
    """Build the GUI command-line parser."""
    parser = argparse.ArgumentParser(prog="python -m spearhead.gui")
    parser.add_argument("--host", default=DEFAULT_HOST, help="host interface for the local GUI server")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help="port for the local GUI server")
    parser.add_argument("--reload", action="store_true", help="reload the GUI when source files change")
    parser.add_argument("--no-show", action="store_true", help="do not open a browser window automatically")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the Spearhead GUI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run_gui(host=args.host, port=args.port, reload=args.reload, show=not args.no_show)
    except RuntimeError as exc:
        parser.exit(1, f"{exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
