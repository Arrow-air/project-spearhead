"""Passive third-person FDM display, using a local north/west/up scene."""

import math

import numpy as np

from ..rotations import euler_body_to_ned


NED_TO_SCENE = np.diag([1., -1., -1.])
GRID_SPACING_M = 25.
MAX_TRAIL_POINTS = 121


def scene_pose(state):
    """Rebase horizontal position around the aircraft; keep true altitude.

    Both FRD-to-NED and NED-to-NWU are proper rotations. Sending the matrix
    directly avoids mixing the scene's Euler convention with aerospace Euler.
    """
    x = np.asarray(state, dtype=float)
    if x.shape != (12,) or not np.isfinite(x).all():
        raise ValueError("External view requires a finite 12-state FDM sample")
    return {
        "position": [0., 0., float(-x[2])],
        "rotation": (NED_TO_SCENE @ euler_body_to_ned(*x[6:9])).tolist(),
        "ground_offset": [float(-x[0] % GRID_SPACING_M),
                          float(x[1] % GRID_SPACING_M), 0.],
    }


def trail_points(history, current):
    """At most 121 recent free-flight points in the same rebased frame."""
    rows = [r for r in history if not r["hold"] and r["time_s"] <= current["time_s"]]
    rows = rows[-600::5]
    if not current["hold"] and (not rows or rows[-1]["time_s"] != current["time_s"]):
        rows.append(current)
    origin = np.array([current["state"][0], current["state"][1], 0.])
    return [(NED_TO_SCENE @ (np.asarray(r["state"][:3])-origin)).tolist()
            for r in rows[-MAX_TRAIL_POINTS:]]


def camera_pose(state, mode):
    """Heading-following camera with world-up, preserving visible roll/pitch."""
    offset = np.array([-7., 6., 3.2] if mode != "side" else [0., 10., 2.5])
    heading = -state[8]  # positive NED yaw turns toward negative scene Y
    c, s = math.cos(heading), math.sin(heading)
    offset[:2] = [c*offset[0]-s*offset[1], s*offset[0]+c*offset[1]]
    return [float(offset[0]), float(offset[1]), float(-state[2]+offset[2])]


class ExternalView:
    def __init__(self, ui):
        # NiceGUI stays optional for the pure coordinate helpers above.
        from nicegui.elements.scene.scene_object3d import Object3D

        class Ground(Object3D, component="external_ground.js"):
            pass

        class FlightPath(Object3D, component="external_path.js"):
            def __init__(self):
                super().__init__([])

            def set_points(self, points):
                self.args[0] = points  # preserve the trail on WebGL remount
                self.run_method("set_points", points)

        self.last_sample = None
        self.last_camera_mode = None
        self.state = None
        self.scene = ui.scene(width=480, height=520, grid=False,
                              camera=ui.scene.perspective_camera(fov=45, near=.1, far=5000),
                              background_color="#92b5cc", fps=20, click_events=[])
        self.scene.classes("w-full rounded-lg overflow-hidden").style("height:520px")
        self.scene.props('aria-label="FDM external aircraft view"')
        with self.scene:
            self.ground = Ground().with_name("Ground · 25 m grid")
            self.path = FlightPath().with_name("FDM flight path")
            self.aircraft = self.scene.group().with_name("FDM aircraft · schematic")
            with self.aircraft:
                self._airframe()
            self.aircraft.visible(False)
        with ui.row().classes("w-full items-center justify-between gap-2"):
            self.camera_mode = ui.toggle({"chase": "Follow", "side": "Side", "orbit": "Orbit"},
                                         value="chase", on_change=lambda: self.reset_camera())
            self.camera_mode.props("dense no-caps")
            ui.button("Reset view", icon="center_focus_strong", on_click=self.reset_camera).props("flat dense no-caps")
        self.readout = ui.label("Waiting for FDM").classes("font-mono text-xs text-slate-300")
        ui.label("Schematic aircraft · 25 m ground grid · Orbit: drag to rotate, scroll to zoom").classes("fd-muted text-xs")
        self.scene.move_camera(-7, 6, 103.2, look_at_x=0, look_at_y=0, look_at_z=100, duration=0)

    def _airframe(self):
        """Illustrative QuadPlane silhouette in body FRD, not a CAD model."""
        s = self.scene
        # A light upper skin, dark canopy and orange nose make attitude readable.
        s.sphere(1).scale(.84, .12, .14).material("#e2e8f0")
        s.sphere(1).scale(.28, .085, .065).move(.27, 0, -.12).material("#172d43")
        s.sphere(1).scale(.12, .125, .135).move(.72, 0, 0).material("#f59e0b")
        wing = [[.28, -.25], [.09, -1.975], [-.19, -1.975], [-.27, -.25],
                [-.27, .25], [-.19, 1.975], [.09, 1.975], [.28, .25]]
        s.extrusion(wing, .035).move(0, 0, -.04).material("#cbd5e1")
        for side in (-1, 1):
            s.box(.17, .22, .045).move(-.05, side*1.84, -.04).material(
                "#ef4444" if side < 0 else "#22c55e")
            s.box(1.8, .045, .045).move(-.15, side*.62, 0).material("#334155")
            # V-tail tips rise toward body -Z.
            s.box(.34, .75, .025).rotate(-side*math.pi/4, 0, 0).move(
                -1.02, side*.27, -.25).material("#cbd5e1")
            for station in (-.70, .55):
                s.cylinder(.055, .055, .10).move(station, side*.62, -.065).material("#172d43")
                s.box(.56, .025, .012).move(station, side*.62, -.13).material("#64748b")
        s.box(.025, .64, .035).move(-.9, 0, 0).material("#64748b")

    def reset_camera(self):
        self.last_camera_mode = None
        if self.state is not None:
            self._move_camera()

    def _move_camera(self):
        mode = self.camera_mode.value
        position = camera_pose(self.state, mode)
        self.scene.move_camera(*position, look_at_x=0, look_at_y=0,
                               look_at_z=float(-self.state[2]), duration=0)
        self.last_camera_mode = mode

    def update(self, snapshot):
        row = snapshot["truth"]
        if row is None:
            self.aircraft.visible(False)
            if self.last_sample is not None:
                self.path.set_points([])
            self.last_sample = self.state = None
            self.readout.text = "Waiting for FDM"
            return
        stamp = (snapshot["run_id"], row["time_s"])
        status = "LIVE" if snapshot["active"] else "LAST SAMPLE"
        if snapshot["stale_truth"]:
            status = "STALE · LAST SAMPLE"
        heading_deg = round(math.degrees(row["state"][8]) % 360, 1) % 360
        self.readout.text = (f"{status} · T+{row['time_s']:.2f} s · "
                             f"HDG {heading_deg:05.1f}°")
        if stamp == self.last_sample:
            return
        self.state = row["state"]
        pose = scene_pose(self.state)
        self.aircraft.visible(True).move(*pose["position"]).rotate_R(pose["rotation"])
        self.ground.move(*pose["ground_offset"])
        self.path.set_points(trail_points(snapshot["history"], row))
        if self.camera_mode.value != "orbit" or self.last_camera_mode != "orbit":
            self._move_camera()
        self.last_sample = stamp
