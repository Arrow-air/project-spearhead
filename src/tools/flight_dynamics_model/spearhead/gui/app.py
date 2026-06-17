"""Minimal NiceGUI shell for the Spearhead Flight Dynamics Model."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Any

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
    run_stability_for_gui,
    stability_mode_table_rows,
    static_derivative_table_rows,
)
from .plotting import simulation_time_history_figure, stability_eigenvalue_figure
from .state import APP_STATE
from ..config import ControlInputCommand, SimulationConfig


FDM_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_SCENARIOS_DIR = FDM_ROOT / "scenarios"
DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8080

_pages_registered = False


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


def _render_page_shell(ui) -> None:
    ui.colors(primary="#2563eb", secondary="#475569", accent="#0f766e")
    with ui.header().classes("items-center justify-between"):
        ui.label("Spearhead FDM").classes("text-lg font-semibold")
        with ui.row().classes("items-center gap-2"):
            ui.link("Scenario", "/").classes("text-white no-underline")
            ui.link("Simulation", "/simulation").classes("text-white no-underline")
            ui.link("Stability", "/stability").classes("text-white no-underline")
    ui.separator()


def _register_scenario_page(ui) -> None:
    @ui.page("/")
    def scenario_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = next(iter(scenario_options), None)

        with ui.column().classes("w-full max-w-5xl mx-auto gap-4 p-4"):
            ui.label("Scenario").classes("text-2xl font-semibold")
            status = ui.label("Load an example YAML scenario or upload one.").classes("text-slate-600")
            columns = [
                {"name": "field", "label": "Field", "field": "field", "align": "left"},
                {"name": "value", "label": "Value", "field": "value", "align": "left"},
            ]
            table = ui.table(columns=columns, rows=[], row_key="field").classes("w-full")

            def show_summary(source: str, summary: dict[str, Any]) -> None:
                table.rows = _summary_rows(summary)
                table.update()
                status.text = f"Loaded {source}"

            def show_error(message: str) -> None:
                table.rows = []
                table.update()
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

            with ui.row().classes("items-end gap-3"):
                selected = ui.select(
                    options=list(scenario_options),
                    value=selected_name,
                    label="Example scenario",
                ).classes("min-w-80")
                ui.button("Load", on_click=load_selected)

            ui.upload(
                label="Upload YAML scenario",
                auto_upload=True,
                on_upload=load_uploaded,
            ).props("accept=.yaml,.yml").classes("w-full")

            if selected_name is not None:
                load_selected()
            else:
                status.text = f"No example scenarios found in {EXAMPLE_SCENARIOS_DIR}"


def _register_simulation_page(ui) -> None:
    @ui.page("/simulation")
    def simulation_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = APP_STATE.scenario_source
        if selected_name not in scenario_options:
            selected_name = next(iter(scenario_options), None)

        with ui.column().classes("w-full max-w-6xl mx-auto gap-4 p-4"):
            ui.label("Simulation").classes("text-2xl font-semibold")
            status = ui.label("Load a scenario, edit safe options, then run the shared backend.").classes("text-slate-600")

            try:
                active_config = _active_or_example_config(selected_name)
            except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                status.text = str(exc)
                active_config = None

            with ui.row().classes("items-end gap-3"):
                selected = ui.select(
                    options=list(scenario_options),
                    value=selected_name,
                    label="Example scenario",
                ).classes("min-w-80")

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
                        status.text = f"Loaded {selected.value}"
                    except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                        status.text = f"Failed to load scenario: {exc}"
                        ui.notify(status.text, type="negative")

                ui.button("Load", on_click=load_selected_for_form)

            if active_config is None:
                return

            def config_for_form() -> SimulationConfig:
                return _active_or_example_config(selected.value if selected.value else None)

            with ui.grid(columns=3).classes("w-full gap-3"):
                name_input = ui.input("Scenario name", value=active_config.name)
                duration_input = ui.number("Duration [s]", value=active_config.duration, min=0.001, step=0.5)
                dt_input = ui.input("Time step dt [s]", value="" if active_config.dt is None else str(active_config.dt))
                n_points_input = ui.input(
                    "Sample count",
                    value="" if active_config.n_points is None else str(active_config.n_points),
                )
                trim_toggle = ui.switch("Start from trim", value=active_config.start_from_trim)
                trim_speed_input = ui.number("Trim airspeed [m/s]", value=active_config.trim_target_speed, min=0.001)
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
            summary_table = ui.table(columns=summary_columns, rows=[], row_key="field").classes("w-full")
            plot_container = ui.column().classes("w-full")
            export_row = ui.row().classes("gap-3")
            export_row.visible = False

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
                        "Download CSV",
                        on_click=lambda: ui.download(
                            export_simulation_csv_text(result),
                            filename=f"{result.config.name}_simulation.csv",
                        ),
                    )
                    ui.button(
                        "Download JSON",
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
                    summary_table.update()
                    plot_container.clear()
                    with plot_container:
                        ui.plotly(simulation_time_history_figure(result)).classes("w-full")
                    render_exports()
                    status.text = f"Simulation complete: success={result.success}"
                    if not result.success:
                        ui.notify(result.message, type="warning")
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    message = f"Simulation failed: {exc}"
                    APP_STATE.set_simulation_error(message)
                    summary_table.rows = []
                    summary_table.update()
                    plot_container.clear()
                    render_exports()
                    status.text = message
                    ui.notify(message, type="negative")

            ui.button("Run Simulation", on_click=run_clicked).classes("w-fit")


def _register_stability_page(ui) -> None:
    @ui.page("/stability")
    def stability_page() -> None:
        _render_page_shell(ui)
        scenario_options = _scenario_options()
        selected_name = APP_STATE.scenario_source
        if selected_name not in scenario_options:
            selected_name = next(iter(scenario_options), None)

        with ui.column().classes("w-full max-w-6xl mx-auto gap-4 p-4"):
            ui.label("Stability").classes("text-2xl font-semibold")
            status = ui.label("Load a scenario, set finite-difference options, then analyze stability.").classes(
                "text-slate-600"
            )

            try:
                active_config = _active_or_example_config(selected_name)
            except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                status.text = str(exc)
                active_config = None

            with ui.row().classes("items-end gap-3"):
                selected = ui.select(
                    options=list(scenario_options),
                    value=selected_name,
                    label="Example scenario",
                ).classes("min-w-80")

                def load_selected_for_analysis() -> None:
                    try:
                        if not selected.value:
                            raise ValueError("No example scenario is selected.")
                        config = _load_example_config(selected.value)
                        scenario_label.text = f"Active scenario: {config.name}"
                        status.text = f"Loaded {selected.value}"
                    except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                        status.text = f"Failed to load scenario: {exc}"
                        ui.notify(status.text, type="negative")

                ui.button("Load", on_click=load_selected_for_analysis)

            if active_config is None:
                return

            scenario_label = ui.label(f"Active scenario: {active_config.name}").classes("text-slate-700")

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
            summary_table = ui.table(columns=summary_columns, rows=[], row_key="field").classes("w-full")

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
                {"name": "oscillation_period_s", "label": "Period [s]", "field": "oscillation_period_s", "align": "right"},
                {"name": "time_to_half_s", "label": "Half [s]", "field": "time_to_half_s", "align": "right"},
                {"name": "time_to_double_s", "label": "Double [s]", "field": "time_to_double_s", "align": "right"},
                {"name": "classification", "label": "Class", "field": "classification", "align": "left"},
            ]
            modes_table = ui.table(columns=mode_columns, rows=[], row_key="eigenvalue_label").classes("w-full")

            derivative_columns = [
                {"name": "name", "label": "Derivative", "field": "name", "align": "left"},
                {"name": "value", "label": "Value", "field": "value", "align": "right"},
                {"name": "alpha_deg", "label": "Alpha [deg]", "field": "alpha_deg", "align": "right"},
                {"name": "beta_deg", "label": "Beta [deg]", "field": "beta_deg", "align": "right"},
                {"name": "source", "label": "Source", "field": "source", "align": "left"},
            ]
            derivatives_table = ui.table(columns=derivative_columns, rows=[], row_key="name").classes("w-full")
            plot_container = ui.column().classes("w-full")
            export_row = ui.row().classes("gap-3")
            export_row.visible = False

            def render_exports() -> None:
                export_row.clear()
                result = APP_STATE.stability_result
                if result is None:
                    export_row.visible = False
                    return
                export_row.visible = True
                with export_row:
                    ui.button(
                        "Download JSON",
                        on_click=lambda: ui.download(
                            export_stability_json_text(result),
                            filename=f"{result.config.name}_stability.json",
                        ),
                    )
                    ui.button(
                        "Download Modes CSV",
                        on_click=lambda: ui.download(
                            export_stability_modes_csv_text(result),
                            filename=f"{result.config.name}_stability_modes.csv",
                        ),
                    )
                    ui.button(
                        "Download Markdown",
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
                    summary_table.update()
                    modes_table.rows = stability_mode_table_rows(result)
                    modes_table.update()
                    derivatives_table.rows = static_derivative_table_rows(result)
                    derivatives_table.update()
                    plot_container.clear()
                    with plot_container:
                        ui.plotly(stability_eigenvalue_figure(result)).classes("w-full")
                    render_exports()
                    status.text = "Stability analysis complete: success=True"
                except Exception as exc:  # noqa: BLE001 - GUI should surface validation errors.
                    message = f"Stability analysis failed: {exc}"
                    APP_STATE.set_stability_error(message)
                    summary_table.rows = [{"field": "success", "value": "False"}, {"field": "error", "value": str(exc)}]
                    summary_table.update()
                    modes_table.rows = []
                    modes_table.update()
                    derivatives_table.rows = []
                    derivatives_table.update()
                    plot_container.clear()
                    render_exports()
                    status.text = message
                    ui.notify(message, type="negative")

            ui.button("Run Stability Analysis", on_click=run_clicked).classes("w-fit")


def create_app() -> None:
    """Register GUI pages with NiceGUI."""
    global _pages_registered
    if _pages_registered:
        return
    ui = _require_nicegui()
    _register_scenario_page(ui)
    _register_simulation_page(ui)
    _register_stability_page(ui)
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
