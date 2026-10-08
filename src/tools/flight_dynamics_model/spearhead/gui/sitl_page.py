"""Web Flight Deck for the external ArduPlane runner, without a second FDM."""

import asyncio
import math

from .flight_deck import attitude_html
from .plotting import _flight_deck_base_figure
from .sitl_session import SCENARIOS


def sitl_figures(snapshot):
    truth, estimates = snapshot["history"], snapshot["estimates"]
    if truth:
        estimates = [row for row in estimates if row["time_s"] >= truth[0]["time_s"]]
    figures = {}
    for axis, index in [("Roll", 6), ("Pitch", 7)]:
        go, figure = _flight_deck_base_figure(f"{axis} · FDM / ArduPilot")
        figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth],
                                   y=[math.degrees(r["state"][index]) for r in truth],
                                   name="FDM", line={"color": "#38bdf8"}))
        figure.add_trace(go.Scatter(x=[r["time_s"] for r in estimates],
                                   y=[r[axis.lower()+"_deg"] for r in estimates],
                                   name="ArduPilot", line={"color": "#fb923c"}))
        figure.update_yaxes(title_text="deg")
        figures[axis] = figure
    go, figure = _flight_deck_base_figure("Firmware outputs")
    for index, label in enumerate(["Aileron", "Elevator", "Pusher", "Rudder"]):
        figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth],
                                   y=[r["pwm"][index] for r in truth], name=label))
    figure.update_yaxes(title_text="PWM · µs")
    figures["PWM"] = figure
    go, figure = _flight_deck_base_figure("FDM speed / altitude")
    for label, values in [("TAS · m/s", [math.sqrt(sum(v*v for v in r["state"][3:6])) for r in truth]),
                          ("Altitude · m", [-r["state"][2] for r in truth])]:
        figure.add_trace(go.Scatter(x=[r["time_s"] for r in truth], y=values, name=label))
    figures["Energy"] = figure
    for figure in figures.values():
        figure.update_layout(margin={"l": 42, "r": 14, "t": 72, "b": 34},
                             title={"y": .98, "x": .03})
    return figures


def register_sitl_page(ui, manager, render_shell):
    @ui.page("/sitl")
    def sitl_page():
        render_shell(ui)
        with ui.column().classes("fd-console w-full gap-4 p-4"):
            with ui.row().classes("w-full fd-panel p-4 items-center justify-between"):
                with ui.column().classes("gap-1"):
                    ui.label("ARDUPILOT FLIGHT DECK").classes("text-xl font-semibold text-sky-100")
                    ui.label("ArduPlane 4.7.1 · October AeroDB · reference design").classes("fd-muted")
                status = ui.label("READY").classes("font-mono text-sky-200 text-lg")
            with ui.row().classes("w-full fd-panel p-4 items-center gap-4"):
                scenario = ui.select(SCENARIOS, value="fbwa", label="Flight test").props("dark filled").classes("w-80")
                start = ui.button("Start flight", icon="play_arrow")
                stop = ui.button("Stop run", icon="stop").props("color=secondary")
                progress = ui.label("Select a flight test to begin.").classes("text-slate-300")
            with ui.row().classes("w-full items-center gap-3"):
                recording = ui.select(manager.completed_runs(), label="Recorded flights").props("dark filled dense").classes("w-80")
                open_recording = ui.button("Open recording", icon="history").props("color=secondary")
            ui.label("Diagnostic CG and provisional aircraft model. Control data cover β=0 only; "
                     "flight outside that plane uses endpoint assumptions.").classes("text-sm text-amber-200")
            with ui.grid(columns=2).classes("w-full gap-4"):
                with ui.column().classes("fd-panel p-4 w-full"):
                    ui.label("FDM · PHYSICAL STATE").classes("text-sm font-semibold text-sky-200")
                    fdm_display = ui.html(attitude_html(pitch_deg=0, roll_deg=0)).classes("w-full")
                    truth_readout = ui.label("Waiting for FDM").classes("font-mono text-sm text-slate-300")
                with ui.column().classes("fd-panel p-4 w-full"):
                    ui.label("ARDUPILOT · AHRS / EKF3").classes("text-sm font-semibold text-orange-200")
                    ekf_display = ui.html(attitude_html(pitch_deg=0, roll_deg=0)).classes("w-full")
                    estimate_readout = ui.label("Waiting for EKF3").classes("font-mono text-sm text-slate-300")
            connection = ui.label("Firmware not running").classes("font-mono text-sm text-slate-300")
            warning = ui.label("").classes("text-amber-200 text-sm")
            error = ui.label("").classes("text-red-300 text-sm whitespace-pre-wrap")
            report = ui.label("").classes("text-sky-100 text-sm")
            with ui.grid(columns=2).classes("w-full gap-4"):
                plots = {name: ui.plotly(figure).classes("w-full")
                         for name, figure in sitl_figures(manager.snapshot()).items()}
            run_label = ui.label("").classes("fd-muted text-xs")
            ui.label("The backend owns the run. Closing this tab keeps the flight running; "
                     "Stop run ends this simulation session.").classes("fd-muted text-sm")

        async def start_clicked():
            start.disable()
            try:
                await asyncio.to_thread(manager.start, scenario.value)
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            refresh()

        async def stop_clicked():
            stop.disable()
            await asyncio.to_thread(manager.stop)
            refresh()

        start.on_click(start_clicked)
        stop.on_click(stop_clicked)
        async def recording_clicked():
            if recording.value is None:
                return
            try:
                await asyncio.to_thread(manager.load_completed, recording.value)
            except Exception as exc:
                ui.notify(str(exc), type="negative")
            refresh()
        open_recording.on_click(recording_clicked)
        last_plot = [None, None]

        def refresh():
            view = manager.snapshot()
            active = view["active"]
            start.set_enabled(not active)
            stop.set_enabled(active)
            scenario.set_enabled(not active)
            recording.set_enabled(not active)
            open_recording.set_enabled(not active)
            options = manager.completed_runs()
            if options != recording.options:
                recording.set_options(options)
            status.text = view["status"]
            error.text = view.get("error") or ""
            row = view["truth"]
            if row is None:
                progress.text = "Starting firmware and FDM…" if active else "Select a flight test to begin."
                connection.text = "Waiting for firmware" if active else "Firmware not running"
                truth_readout.text = "Waiting for FDM"
                estimate_readout.text = "Waiting for ArduPilot"
                fdm_display.set_visibility(False)
                ekf_display.set_visibility(False)
                report.text = ""
                run_label.text = "" if view["run_id"] is None else f"Run: {view['run_id']}"
                return
            t = row["time_s"]
            free = view["free_time_s"]
            progress.text = f"T+{t:.1f} s · " + ("calibrating / arming" if free is None else f"{free:.1f} s free flight")
            x = row["state"]
            fdm_display.set_visibility(True)
            roll, pitch = math.degrees(x[6]), math.degrees(x[7])
            fdm_display.content = attitude_html(pitch_deg=pitch, roll_deg=roll)
            fdm_display.update()
            tas = math.sqrt(sum(v*v for v in x[3:6]))
            truth_readout.text = f"TAS {tas:.2f} m/s · ALT {-x[2]:.2f} m · α {row['alpha_deg']:.3f}° · β {row['beta_deg']:.3f}°"
            estimate = view["estimate"]
            ekf_display.set_visibility(estimate is not None)
            if estimate is not None:
                ekf_display.content = attitude_html(pitch_deg=estimate["pitch_deg"], roll_deg=estimate["roll_deg"])
                ekf_display.update()
                estimate_readout.text = (f"Roll {estimate['roll_deg']:+.3f}° · Pitch {estimate['pitch_deg']:+.3f}° "
                                         f"· sampled at T+{estimate['time_s']:.2f} s")
            else:
                estimate_readout.text = "Estimate unavailable / stale"
            connection.text = (("LIVE · " if active else "RUN ENDED · Last sample: ") +
                               f"{'ARMED' if view['armed'] else 'DISARMED'} · "
                               f"Mode { {0: 'MANUAL', 5: 'FBWA'}.get(view['mode'], view['mode'])} · "
                               f"EKF flags {view['ekf_flags']} · Physics {row['physics_rate_hz']} Hz · "
                               f"{'BOOTSTRAP HOLD' if row['hold'] else 'FIRMWARE CONTROL'}")
            notices = []
            if view["stale_truth"]:
                notices.append("FDM connection stale")
            if view["stale_mavlink"]:
                notices.append("MAVLink connection stale")
            if row["control_coverage_exceeded"]:
                notices.append("Outside measured control-data plane")
            warning.text = " · ".join(notices)
            run_label.text = f"Run: {view['run_id']} · Logs: {view['directory']}"
            result = view["report"]
            if result is not None:
                rms = result.get("attitude_rms_error_deg")
                report.text = (f"{view['status']} · " + (f"Roll / pitch RMS {rms[0]:.4f}° / {rms[1]:.4f}°" if rms else result.get("error", "Run ended")))
            else:
                report.text = ""
            # Display updates only; no physics or MAVLink calls from a browser timer.
            stamp = (view["run_id"], int(t))
            if tuple(last_plot) != stamp:
                for name, figure in sitl_figures(view).items():
                    plots[name].figure = figure
                    plots[name].update()
                last_plot[:] = stamp

        refresh()
        ui.timer(.25, refresh)
