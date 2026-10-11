"""Web Flight Deck for the external ArduPlane runner, without a second FDM."""

import asyncio
import math

from .flight_deck import attitude_html
from .external_view import ExternalView
from .plotting import _flight_deck_base_figure
from .flight_setup import AIRCRAFT, POINTS, MANEUVERS, MIN_DURATION, FlightSetup
from ..sitl.ground_control import QUADPLANE_MODES, PLANE_MODES, MODE_HELP, MODE_NAMES, session_modes


def sitl_figures(snapshot):
    truth = snapshot["history"]
    figures = {}
    for axis, index in [("Roll", 6), ("Pitch", 7)]:
        go, figure = _flight_deck_base_figure(f"{axis} · FDM")
        figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth],
                                   y=[math.degrees(r["state"][index]) for r in truth],
                                   name="FDM", line={"color": "#38bdf8"}))
        figure.update_yaxes(title_text="deg")
        figures[axis] = figure
    go, figure = _flight_deck_base_figure("Firmware outputs")
    for index, label in enumerate(["Aileron", "Elevator", "Pusher", "Rudder"]):
        figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth],
                                   y=[r["pwm"][index] for r in truth], name=label))
    figure.update_yaxes(title_text="PWM · µs")
    figures["PWM"] = figure
    go, figure = _flight_deck_base_figure("Physical elevator")
    figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth],
                               y=[math.degrees(r["control"]["de"]) for r in truth],
                               name="Elevator · deg", line={"color": "#fbbf24"}))
    manifest = snapshot.get("manifest")
    if manifest is not None:
        figure.add_hline(y=math.degrees(manifest["trim_control"]["de"]),
                         line_dash="dot", line_color="#94a3b8", annotation_text="Trim")
    figure.update_yaxes(title_text="Surface angle · deg")
    figures["Elevator"] = figure
    go, figure = _flight_deck_base_figure("FDM speed / altitude")
    for label, values in [("TAS · m/s", [math.sqrt(sum(v*v for v in r["state"][3:6])) for r in truth]),
                          ("Altitude · m", [-r["state"][2] for r in truth])]:
        figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth], y=values, name=label))
    figures["Energy"] = figure
    if (snapshot.get("manifest") or {}).get("config", {}).get("flight_regime") == "hover":
        go, figure = _flight_deck_base_figure("Quad lift motors · firmware outputs")
        for index, label in [(4, "FR · M1"), (7, "RR · M4"), (5, "RL · M2"), (6, "FL · M3")]:
            figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth], y=[r["pwm"][index] for r in truth], name=label))
        figure.update_yaxes(title_text="Lift PWM · µs")
        figures["PWM"] = figure
    for figure in figures.values():
        figure.update_layout(margin={"l": 42, "r": 14, "t": 72, "b": 34},
                             title={"y": .98, "x": .03})
    return figures


def register_sitl_page(ui, manager, render_shell):
    @ui.page("/embed")
    @ui.page("/")
    @ui.page("/sitl")
    def sitl_page():
        render_shell(ui)
        preset = FlightSetup(**(manager.snapshot(include_history=False).get("flight_setup") or {}))
        with ui.dialog() as history_dialog, ui.card().classes("fd-panel text-slate-200 w-full max-w-xl p-5 gap-4"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label("Recorded flights").classes("text-lg font-semibold")
                ui.button(icon="close", on_click=history_dialog.close).props('flat round dense aria-label="Close recordings"')
            ui.label("Review the last physical state, messages and plots. Recording controls are read-only.").classes("fd-muted text-sm")
            recording = ui.select(manager.completed_runs(), label="Recorded flights").props("dark filled dense").classes("w-full")
            open_recording = ui.button("Open recording", icon="history")
        embedded = ui.context.client.request.url.path == "/embed"
        if embedded:
            ui.add_body_html('''<script>
            (() => {
              let previous = 0;
              const origins = ['https://www.arrowair.com', 'https://arrowair.com', location.origin];
              if (window.parent === window) return;
              let parentOrigin;
              try { parentOrigin = new URL(document.referrer || location.origin).origin; } catch { return; }
              if (!origins.includes(parentOrigin)) return;
              const connect = () => {
                const root = document.querySelector('.fd-embed');
                if (!root) { requestAnimationFrame(connect); return; }
                new ResizeObserver(() => {
                  const height = Math.ceil(root.getBoundingClientRect().height);
                  if (height === previous) return;
                  previous = height;
                  window.parent.postMessage({type:'arrow-flight-deck:size',height}, parentOrigin);
                }).observe(root);
              };
              connect();
            })();
            </script>''')
        with ui.column().classes("fd-console fd-workspace w-full gap-3 p-4" + (" fd-embed" if embedded else "")):
            with ui.row().classes("w-full items-center justify-between"):
                with ui.column().classes("gap-1"):
                    ui.label("Spearhead flight simulator").classes("text-2xl fd-title text-sky-100")
                    aircraft_readout = ui.label("").classes("fd-muted text-xs")
                with ui.column().classes("items-end gap-1"):
                    status = ui.label("READY").classes("text-xs font-semibold rounded-full bg-slate-800 px-3 py-1 text-sky-200")
                    progress = ui.label("").classes("fd-muted text-xs")
            with ui.column().classes("w-full fd-panel p-4 gap-3"):
                with ui.element("div").classes("fd-setup-grid"):
                    aircraft = ui.select(AIRCRAFT, value=preset.aircraft, label="Aircraft").props("dark outlined dense").classes("w-full")
                    initial_point = ui.select(POINTS[preset.aircraft], value=preset.initial_point, label="Initial point").props("dark outlined dense").classes("w-full")
                    maneuver = ui.select({key: value for key, value in MANEUVERS[preset.aircraft].items()
                                          if (key == "takeoff_hover_land") == (preset.initial_point == "ground")},
                                         value=preset.maneuver, label="Maneuver").props("dark outlined dense").classes("w-full")
                    duration = ui.number("Simulation time (s)", value=preset.duration_s, min=15, max=600, step=5).props("dark outlined dense").classes("w-full")
                with ui.row().classes("w-full items-center gap-3"):
                    with ui.element("div").classes("fd-condition-grid"):
                        altitude = ui.number("Initial altitude (m)", value=preset.altitude_m, min=6, max=110, step=1).props("dark outlined dense").classes("w-full")
                        speed = ui.number("Initial airspeed (m/s)", value=preset.speed_m_s, min=0, max=30, step=1).props("dark outlined dense").classes("w-full")
                    ui.space()
                    start = ui.button("Start flight", icon="play_arrow")
                    stop = ui.button("Stop run", icon="stop").props("outline color=red-4")
                    history_button = ui.button("Recordings", icon="history", on_click=history_dialog.open).props("flat color=blue-grey-3")
                setup_hint = ui.label("").classes("fd-muted text-xs")
                setup_error = ui.label("").classes("text-amber-200 text-sm").props('role="status"')
            with ui.column().classes("w-full fd-panel p-3 gap-2"):
                with ui.row().classes("w-full items-center gap-3"):
                    with ui.column().classes("gap-0 mr-4"):
                        current_mode = ui.label("Mode: —").classes("text-lg font-semibold text-sky-100")
                        arm_state = ui.label("Arming: unknown").classes("fd-muted text-xs")
                    mode_select = ui.select(QUADPLANE_MODES, value=18, label="Requested mode").props("dark filled dense").classes("w-64")
                    apply_mode = ui.button("Change mode", icon="swap_horiz")
                    command_status = ui.label("").classes("text-sm text-sky-200")
                mode_help = ui.label(MODE_HELP[18]).classes("fd-muted text-xs")
                control_hint = ui.label("").classes("text-xs text-slate-300")
            with ui.row().classes("w-full fd-health gap-x-5 gap-y-1 items-center px-1"):
                link_status = ui.label("Link: offline").classes("text-slate-300")
                gps_status = ui.label("GPS: unknown").classes("text-slate-300")
                ekf_status = ui.label("EKF: unknown").classes("text-slate-300")
                system_status = ui.label("System: unknown").classes("text-slate-300")
                battery_status = ui.label("Simulated battery: unknown").classes("text-slate-300")
            error = ui.label("").classes("w-full fd-alert fd-error text-red-200 text-sm whitespace-pre-wrap").props('role="alert"')
            warning = ui.label("").classes("w-full fd-alert text-amber-200 text-xs").props('role="status"')
            with ui.element("div").classes("fd-instruments w-full gap-3"):
                with ui.column().classes("fd-panel p-3 w-full min-w-0 gap-2"):
                    ui.label("Attitude").classes("text-sm font-semibold text-sky-100")
                    fdm_display = ui.html(attitude_html(pitch_deg=0, roll_deg=0)).classes("w-full")
                    with ui.column().classes("w-full items-center justify-center gap-2 rounded-lg bg-slate-950").style("height:330px") as attitude_empty:
                        ui.icon("flight", size="40px").classes("text-slate-600")
                        ui.label("Ready for your next flight").classes("text-slate-300")
                        ui.label("Start a session or open a recording").classes("fd-muted text-xs")
                    truth_readout = ui.label("Waiting for flight data").classes("font-mono text-sm text-slate-300")
                    ui.label("Physical model state").classes("fd-muted text-xs")
                with ui.column().classes("fd-panel p-3 w-full min-w-0 gap-2"):
                    ui.label("Aircraft view").classes("text-sm font-semibold text-sky-100")
                    external_view = ExternalView(ui)
                    external_view.scene.classes("fd-scene")
            with ui.column().classes("w-full fd-panel gap-0"):
                with ui.tabs().props("no-caps inline-label align=left active-color=light-blue-3").classes("w-full") as detail_tabs:
                    messages_tab = ui.tab("Messages", icon="chat_bubble_outline")
                    plots_tab = ui.tab("Flight plots", icon="show_chart")
                    model_tab = ui.tab("Model & run details", icon="tune")
                with ui.tab_panels(detail_tabs, value=messages_tab).props("animated=false").classes("w-full bg-transparent text-slate-200"):
                    with ui.tab_panel(messages_tab).classes("p-3 gap-3"):
                        with ui.row().classes("w-full items-center justify-between"):
                            message_count = ui.label("No firmware messages yet.").classes("fd-muted text-xs")
                            message_filter = ui.select(["All messages", "Warnings and errors"], value="Warnings and errors", label="Show").props("dark filled dense").classes("w-56")
                        message_table = ui.table(columns=[
                            {"name": "time", "label": "Time", "field": "time", "align": "left"},
                            {"name": "level", "label": "Severity", "field": "level", "align": "left"},
                            {"name": "source", "label": "Source", "field": "source", "align": "left"},
                            {"name": "text", "label": "Message", "field": "text", "align": "left", "style": "white-space: normal"},
                        ], rows=[], row_key="id", pagination=5).props("dark flat dense wrap-cells").classes("w-full")
                        message_table.add_slot('body-cell-level', '<q-td :props="props"><q-badge :color="props.row.severity <= 3 ? \'red\' : props.row.severity === 4 ? \'amber-9\' : \'blue-grey\'">{{ props.value }}</q-badge></q-td>')
                    with ui.tab_panel(plots_tab).classes("p-3"):
                        with ui.element("div").classes("grid grid-cols-1 lg:grid-cols-2 w-full gap-3"):
                            plots = {name: ui.plotly(figure).classes("w-full")
                                     for name, figure in sitl_figures(manager.snapshot()).items()}
                    with ui.tab_panel(model_tab).classes("p-4 gap-3"):
                        ui.label("Aircraft & model assumptions").classes("font-semibold")
                        mode_scope = ui.label("").classes("text-sm text-slate-300")
                        mode_limits = ui.label("").classes("text-sm text-amber-200")
                        ui.label("Diagnostic CG and provisional aircraft model. Control data cover β=0 only; flight outside that plane uses endpoint assumptions.").classes("text-sm text-amber-200")
                        travel_readout = ui.label("New flights: assumed ±20° surface travel · measured ±5° · linear extrapolation beyond measured deflection").classes("text-sm text-amber-200")
                        aero_readout = ui.label("").classes("text-sm font-mono text-slate-300")
                        ui.separator()
                        ui.label("Run information").classes("font-semibold")
                        connection = ui.label("Firmware not running").classes("font-mono text-xs text-slate-300")
                        report = ui.label("").classes("text-sky-100 text-sm")
                        run_label = ui.label("").classes("fd-muted text-xs break-all")
                        ui.label("Closing this tab keeps the flight running. Stop run ends the simulation session.").classes("fd-muted text-xs")

        async def start_clicked():
            start.disable()
            try:
                await asyncio.to_thread(manager.start_configured, selected_setup())
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            refresh()

        def selected_setup():
            return FlightSetup(aircraft.value, initial_point.value, maneuver.value,
                               duration.value, altitude.value, speed.value).validate()

        syncing_setup = [False]
        setup_run = [None]

        def aircraft_changed():
            if syncing_setup[0]:
                return
            hover = aircraft.value == "pt1_as_built"
            initial_point.set_options(POINTS[aircraft.value], value="hover20" if hover else "trim28")
            maneuver.set_options({key: value for key, value in MANEUVERS[aircraft.value].items() if key != "takeoff_hover_land"}, value="interactive")
            altitude.value, speed.value = (20, 0) if hover else (100, 28)
            refresh()

        def point_changed():
            if syncing_setup[0]:
                return
            point = initial_point.value
            if point == "ground":
                maneuver.set_options({"takeoff_hover_land": MANEUVERS[aircraft.value]["takeoff_hover_land"]}, value="takeoff_hover_land")
                altitude.value = min(30, max(3, altitude.value or 20))
            else:
                options = {key: value for key, value in MANEUVERS[aircraft.value].items() if key != "takeoff_hover_land"}
                maneuver.set_options(options, value=maneuver.value if maneuver.value in options else "interactive")
            if point.startswith("hover"):
                altitude.value = int(point[5:])
            elif point.startswith("trim"):
                speed.value = int(point[4:])
            refresh()

        def conditions_changed():
            if syncing_setup[0]:
                return
            point = initial_point.value
            if point.startswith("hover") and altitude.value != int(point[5:]):
                initial_point.value = "custom_hover"
            elif point.startswith("trim") and speed.value != int(point[4:]):
                initial_point.value = "custom_trim"
            refresh()

        def maneuver_changed():
            if syncing_setup[0]:
                return
            minimum = MIN_DURATION.get(maneuver.value, 15)
            if duration.value is not None and duration.value < minimum:
                duration.value = minimum
            refresh()

        async def stop_clicked():
            stop.disable()
            await asyncio.to_thread(manager.stop)
            refresh()

        start.on_click(start_clicked)
        stop.on_click(stop_clicked)
        async def mode_clicked():
            try:
                manager.set_mode(mode_select.value)
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            refresh()
        apply_mode.on_click(mode_clicked)
        mode_select.on_value_change(lambda e: mode_help.set_text(MODE_HELP.get(e.value, "Choose a flight mode.")))
        async def recording_clicked():
            if recording.value is None:
                return
            try:
                await asyncio.to_thread(manager.load_completed, recording.value)
                history_dialog.close()
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            refresh()
        open_recording.on_click(recording_clicked)
        last_plot = [None, None]
        last_messages = [None]

        def refresh():
            view = manager.snapshot(include_history=False)
            active = view["active"]
            # Another browser or an embedded view can start the shared session.
            # While locked, show its actual recipe instead of this client's draft.
            if active and view.get("flight_setup") and setup_run[0] != view.get("run_id"):
                recipe = FlightSetup(**view["flight_setup"]).validate()
                syncing_setup[0] = True
                try:
                    aircraft.value = recipe.aircraft
                    initial_point.set_options(POINTS[recipe.aircraft], value=recipe.initial_point)
                    maneuver.set_options({key: value for key, value in MANEUVERS[recipe.aircraft].items()
                                          if (key == "takeoff_hover_land") == (recipe.initial_point == "ground")}, value=recipe.maneuver)
                    altitude.max = 110 if recipe.aircraft == "pt1_as_built" else 1000
                    speed.min = 0 if recipe.aircraft == "pt1_as_built" else 22
                    duration.min = MIN_DURATION.get(recipe.maneuver, 15)
                    duration.value = recipe.duration_s
                    altitude.value = recipe.altitude_m
                    speed.value = recipe.speed_m_s
                    setup_run[0] = view.get("run_id")
                finally:
                    syncing_setup[0] = False
            hover_setup = aircraft.value == "pt1_as_built"
            ground_setup = initial_point.value == "ground"
            try:
                selected_setup()
                setup_error.text = ""
            except ValueError as exc:
                setup_error.text = str(exc)
            setup_error.set_visibility(bool(setup_error.text))
            start.set_enabled(not active and not setup_error.text)
            stop.set_enabled(active)
            for field in (aircraft, initial_point, maneuver, duration, altitude):
                field.set_enabled(not active)
            speed.set_enabled(not active and not hover_setup)
            # Keep numeric limits numeric; props strings also reach server-side sanitization.
            altitude.min = 3 if ground_setup else 6
            altitude.max = 30 if ground_setup else 110 if hover_setup else 1000
            altitude.label = "Takeoff altitude (m)" if ground_setup else "Initial altitude (m)"
            duration.label = "Hover time (s)" if ground_setup else "Simulation time (s)"
            speed.min = 0 if hover_setup else 22
            duration.min = MIN_DURATION.get(maneuver.value, 15)
            setup_hint.text = ("Air-start hover · KDE lift · MAD 150KV / APC 20×12 pusher estimate. PT1 forward starts await trim validation." if hover_setup else
                               "Air-start at six-axis trim · provisional pusher · Spearhead hover awaits MAD performance data.")
            setup_hint.text += " Simulation time is free-flight time; calibration and arming run first."
            if ground_setup:
                setup_hint.text = "Start on the ground → VTOL takeoff → timed hover → landing. Climb and landing are additional to hover time. Flat ground and provisional landing gear."
            if maneuver.value in MIN_DURATION:
                setup_hint.text += f" This maneuver sequence needs at least {MIN_DURATION[maneuver.value]} s."
            recording.set_enabled(not active)
            open_recording.set_enabled(not active and recording.value is not None)
            history_button.set_enabled(not active)
            options = manager.completed_runs()
            if options != recording.options:
                recording.set_options(options)
            status.text = view["status"]
            manifest = view.get("manifest") or {}
            aircraft_readout.text = "ArduPlane 4.7.1 · October AeroDB · Provisional simulation"
            if manifest:
                aircraft_readout.text += " · Displaying " + manifest.get("aircraft", {}).get("label", "recorded aircraft")
            context = manifest if manifest else {"parameter_overlay": {"Q_ENABLE": int(hover_setup)}}
            mode_options = session_modes(context)
            if mode_select.options != mode_options:
                mode_select.set_options(mode_options, value=18 if 18 in mode_options else 5)
            if manifest.get("ground_contact") and view.get("mode") in mode_options:
                mode_select.value = view["mode"]
            mode_help.text = MODE_HELP.get(mode_select.value, "Choose a flight mode.")
            forward_mode = 18 in mode_options and (mode_select.value in PLANE_MODES or mode_select.value == 25)
            if forward_mode:
                mode_help.text += " PT1 electric pusher is estimated; transition and sustained forward flight remain unvalidated."
            mode_help.classes(replace="text-xs text-amber-200" if forward_mode else "fd-muted text-xs")
            mass = manifest.get("aircraft", {}).get("mass_kg", 21)
            mode_scope.text = (f"PT1 Quad-X · {mass:g} kg · four KDE lift motors · neutral pilot sticks" if 18 in mode_options else
                               "Fixed-wing reference aircraft · neutral pilot sticks")
            mode_limits.text = ("PT1 pusher: provisional motor/propeller estimate; sustainable motor-current guideline is 24 A. Transition remains unvalidated. "
                                "No ground contact or mission upload. AUTOLAND is unavailable for QuadPlane." if 18 in mode_options else
                                "Air-start simulation: no runway/contact model, mission upload or guided-target entry.")
            if manifest.get("ground_contact"):
                mode_limits.text = "AUTO flies the built-in takeoff, hover and landing mission. Flat-ground gear is provisional; forward transition remains unvalidated."
                if mode_select.value in (10, 13, 20, 21):
                    mode_help.text = "This ground session includes physical gear contact and an uploaded VTOL mission. Mode changes are disabled while the automatic sequence runs."
            error.text = view.get("error") or ""
            error.set_visibility(bool(error.text))
            historical = not active and view.get("run_id") is not None
            current_mode.text = "Mode: " + (view.get("mode_name") or "—")
            armed = view.get("armed")
            arm_state.text = ("Last sample: " if historical else "") + ("Armed" if armed else "Disarmed" if armed is not None else "Arming: unknown")
            command = view.get("command") or {}
            pending = command.get("status") in ("queued", "sent", "accepted")
            enabled = (active and view.get("interactive") and view["status"] == "FREE FLIGHT"
                       and not view.get("stale_truth") and not view.get("stale_mavlink") and not view.get("error")
                       and view.get("mode_controls_ready", True))
            mode_select.set_enabled(enabled and not pending)
            apply_mode.set_enabled(enabled and not pending)
            command_text = command.get("text", "")
            if command.get("status") in ("rejected", "timeout"):
                target = MODE_NAMES.get(command.get("mode"), "Mode")
                command_text = f"{target} request {command['status']}: {command_text}"
            command_status.text = command_text
            command_status.classes(replace="text-sm text-amber-200" if command.get("status") in ("rejected", "timeout") else "text-sm text-sky-200")
            control_hint.text = ("Mode changes are confirmed by the firmware heartbeat." if enabled else
                                 "Recording: controls disabled." if historical else
                                 "Scripted test: mode is controlled by the test." if active and manifest and not view.get("interactive") else
                                 "Calibrating and arming — mode controls unlock when ready." if active else
                                 "Start an interactive session to change flight modes.")
            stale = view.get("stale_mavlink") or view.get("stale_truth")
            link_status.text = "Link: recorded" if historical else "Link: stale" if active and stale else "Link: live" if active else "Link: offline"
            health = view.get("health") or {}
            fix, satellites = health.get("gps_fix"), health.get("satellites")
            gps_status.text = "GPS: " + ("unknown" if fix is None else f"{'3D fix' if fix >= 3 else 'no 3D fix'} · {satellites if satellites is not None else '—'} satellites")
            ekf = health.get("ekf_ok")
            ekf_status.text = "EKF: " + ("unknown" if ekf is None else "ready" if ekf else "not ready")
            ekf_status.classes(replace="text-amber-200" if ekf is False else "text-slate-300")
            faults = health.get("sensor_faults", [])
            system_status.text = "System: " + (", ".join(faults) if faults else "checks healthy" if health.get("system_known") else "unknown")
            system_status.classes(replace="text-amber-200" if faults else "text-slate-300")
            link_status.classes(replace="text-red-300" if active and stale else "text-slate-300")
            voltage, remaining = health.get("battery_v"), health.get("battery_pct")
            battery_status.text = "Battery (sim): " + ("unknown" if voltage is None else f"{voltage:.1f} V") + ("" if remaining is None else f" · {remaining}%")
            status.text = ("RECORDING · " if historical else "") + view["status"]
            messages = view.get("notifications", [])
            visible = [m for m in messages if message_filter.value == "All messages" or m["severity"] <= 4]
            stamp_messages = (view.get("run_id"), message_filter.value, repr(visible))
            if last_messages[0] != stamp_messages:
                message_table.rows = [{**m, "id": i, "time": f"T+{m['time_s']:.1f} s"} for i, m in reversed(list(enumerate(visible)))]
                message_table.update()
                last_messages[0] = stamp_messages
            warning_count = sum(m['severity'] <= 4 for m in messages)
            messages_tab.props(f'label="Messages · {warning_count} warning events"' if warning_count else 'label="Messages"')
            message_count.text = f"{len(visible)} displayed · {len(messages)} recent messages · latest 200 retained"
            row = view["truth"]
            model = (view.get("manifest") or {}).get("control_model")
            if model:
                travel = model["assumed_travel_deg"]
                travel_readout.text = (f"Assumed travel: elevator ±{travel['elevator']:g}° · aileron ±{travel['aileron']:g}° · rudder ±{travel['rudder']:g}° · "
                                       f"Measured ±5° · {'linear extrapolation' if model['deflection_policy'] == 'linear' else 'endpoint hold'} · uncalibrated")
            elif row is not None:
                travel_readout.text = "Recorded run uses an earlier provisional travel setup; see its manifest. New flights use ±20°."
            external_view.update(view)
            attitude_empty.set_visibility(row is None)
            if row is None:
                progress.text = "Starting firmware and FDM…" if active else "Select a flight test to begin."
                connection.text = "Waiting for firmware" if active else "Firmware not running"
                warning.text = ""
                warning.set_visibility(False)
                truth_readout.text = "Waiting for flight data"
                aero_readout.text = ""
                fdm_display.set_visibility(False)
                report.text = ""
                run_label.text = "" if view["run_id"] is None else f"Run: {view['run_id']}"
                return
            t = row["time_s"]
            free = view["free_time_s"]
            progress.text = f"T+{t:.1f} s · " + ("calibrating / arming on ground" if free is None and manifest.get("ground_contact") else
                                              "calibrating / arming" if free is None else
                                              f"{free:.1f} s since takeoff command" if manifest.get("ground_contact") else f"{free:.1f} s free flight")
            x = row["state"]
            fdm_display.set_visibility(True)
            roll, pitch = math.degrees(x[6]), math.degrees(x[7])
            fdm_display.content = attitude_html(pitch_deg=pitch, roll_deg=roll)
            fdm_display.update()
            tas = math.sqrt(sum(v*v for v in x[3:6]))
            truth_readout.text = f"Airspeed {tas:.2f} m/s · Altitude {-x[2]:.2f} m"
            aero_readout.text = f"α {row['alpha_deg']:.3f}° · β {row['beta_deg']:.3f}° · Elevator {math.degrees(row['control']['de']):+.2f}°"
            connection.text = (("LIVE · " if active else "RUN ENDED · Last sample: ") +
                               f"{'ARMED' if view['armed'] else 'DISARMED'} · "
                               f"Mode {view.get('mode_name') or view['mode']} · "
                               f"EKF flags {view['ekf_flags']} · Physics {row['physics_rate_hz']} Hz · "
                               f"{'BOOTSTRAP HOLD' if row['hold'] else 'FIRMWARE CONTROL'}")
            notices = []
            if active and view["stale_truth"]:
                notices.append("FDM connection stale")
            if active and view["stale_mavlink"]:
                notices.append("MAVLink connection stale")
            if row["control_coverage_exceeded"]:
                notices.append("Outside measured control-data plane")
            extrapolated = (row.get("control_coverage") or {}).get("extrapolated_surfaces", [])
            if extrapolated:
                notices.append("Deflection extrapolated: " + ", ".join(extrapolated))
            warning.text = " · ".join(notices)
            if historical and notices:
                warning.text = "Recorded model warning · " + warning.text
            warning.set_visibility(bool(warning.text))
            run_label.text = f"Run: {view['run_id']} · Logs: {view['directory']}"
            result = view["report"]
            if result is not None:
                rms = result.get("attitude_rms_error_deg")
                if result.get("test_kind") == "quad_hover" and rms:
                    report.text = (f"{view['status']} · Nominal hover · Roll / pitch RMS {rms[0]:.3f}° / {rms[1]:.3f}° · "
                                   f"Altitude error ≤{result['settled_max_altitude_error_m']:.2f} m · "
                                   f"Horizontal drift {result['horizontal_displacement_m']:.1f} m (QHOVER)")
                elif result.get("test_kind") == "manual_elevator_demo":
                    report.text = (f"{view['status']} · MANUAL elevator input · "
                                   f"Pitch excursion {result['pitch_excursion_deg']:.2f}° · "
                                   f"Final pitch error {result['final_5s_pitch_error_deg']:.2f}°")
                else:
                    report.text = (f"{view['status']} · " + (f"FDM target tracking · Roll / pitch RMS {rms[0]:.4f}° / {rms[1]:.4f}°" if rms else result.get("error", "Run ended")))
            else:
                report.text = ""
            # Display updates only; no physics or MAVLink calls from a browser timer.
            stamp = (view["run_id"], int(t/2))
            if tuple(last_plot) != stamp:
                for name, figure in sitl_figures(manager.snapshot()).items():
                    plots[name].figure = figure
                    plots[name].update()
                last_plot[:] = stamp

        aircraft.on_value_change(lambda: aircraft_changed())
        initial_point.on_value_change(lambda: point_changed())
        maneuver.on_value_change(lambda: maneuver_changed())
        altitude.on_value_change(lambda: conditions_changed())
        speed.on_value_change(lambda: conditions_changed())
        duration.on_value_change(lambda: refresh())
        refresh()
        ui.timer(.25, refresh)
