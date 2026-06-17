"""Minimal NiceGUI shell for the Spearhead Flight Dynamics Model."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Any

from .adapters.scenarios import load_scenario_config, scenario_summary


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

    return getattr(upload_event, "name", "uploaded scenario"), scenario_summary(config)


def _render_page_shell(ui) -> None:
    ui.colors(primary="#2563eb", secondary="#475569", accent="#0f766e")
    with ui.header().classes("items-center justify-between"):
        ui.label("Spearhead FDM").classes("text-lg font-semibold")
        with ui.row().classes("items-center gap-2"):
            ui.link("Scenario", "/").classes("text-white no-underline")
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
                    path = scenario_options[selected.value]
                    config = load_scenario_config(path)
                    show_summary(path.name, scenario_summary(config))
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


def create_app() -> None:
    """Register GUI pages with NiceGUI."""
    global _pages_registered
    if _pages_registered:
        return
    ui = _require_nicegui()
    _register_scenario_page(ui)
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
