"""Commands executed by the runner that owns the local firmware connection."""

import json
import time

# Exact mode IDs in the pinned ArduPlane 4.7.1 mode.h (dbe792162d06).
# QAUTOTUNE, THERMAL, ADS-B and AUTOLAND are compiled into our SITL binary.
MODE_NAMES = {0: "MANUAL", 1: "CIRCLE", 2: "STABILIZE", 3: "TRAINING", 4: "ACRO",
              5: "FBWA", 6: "FBWB", 7: "CRUISE", 8: "AUTOTUNE", 10: "AUTO",
              11: "RTL", 12: "LOITER", 13: "TAKEOFF", 14: "AVOID_ADSB",
              15: "GUIDED", 16: "INITIALISING", 17: "QSTABILIZE", 18: "QHOVER",
              19: "QLOITER", 20: "QLAND", 21: "QRTL", 22: "QAUTOTUNE",
              23: "QACRO", 24: "THERMAL", 25: "LOITER_ALT_QLAND", 26: "AUTOLAND"}
VTOL_MODES = {17, 18, 19, 20, 21, 22, 23, 25}
PLANE_MODES = {key: value for key, value in MODE_NAMES.items() if key not in VTOL_MODES | {16}}
# AUTOLAND::_enter explicitly rejects QuadPlanes. INITIALISING is internal only.
QUADPLANE_MODES = {key: value for key, value in MODE_NAMES.items() if key not in (16, 26)}
MODE_HELP = {
    0: "Direct pilot control; no attitude or altitude hold.",
    1: "Circles with a gentle bank; altitude control uses throttle.",
    2: "Levels the wings when sticks are released; pilot controls throttle.",
    3: "Pilot control with roll and pitch limits; pilot controls throttle.",
    4: "Pilot commands rotation rates; pilot controls throttle.",
    5: "Pilot commands bank and pitch; pilot controls throttle.",
    6: "Pilot commands bank and altitude; automatic throttle controls speed.",
    7: "Holds course and altitude; automatic throttle controls speed.",
    8: "Fixed-wing gain learning; requires deliberate pilot manoeuvres. Gains may change in this run.",
    10: "Follows an uploaded mission. This deck does not upload missions yet; firmware may reject or fall back.",
    11: "Returns toward home. QuadPlane return behaviour depends on Q_RTL_MODE.",
    12: "Fixed-wing circles around a location; automatic speed and altitude control.",
    13: "Automatic fixed-wing takeoff; this simulator starts airborne and has no runway/contact model.",
    14: "Fixed-wing ADS-B avoidance; requires configured ADS-B avoidance and traffic. No traffic feed is supplied.",
    15: "Flies to guided targets. Target entry is not implemented here; changing mode alone does not set a destination.",
    17: "Levels in hover; direct lift throttle. The simulator supplies its calculated hover throttle, not altitude hold.",
    18: "Holds hover attitude and altitude; horizontal drift is free.",
    19: "Holds hover position and altitude using GPS/EKF. Centre sticks request position hold.",
    20: "Vertical landing; firmware controls descent. Run stops before the ground because contact is not modelled.",
    21: "Returns home and lands vertically; may use a fixed-wing approach depending on Q_RTL_MODE. No ground contact.",
    22: "VTOL gain learning; excites the aircraft and may change gains in this run. Available in our pinned SITL build.",
    23: "Hover rotation-rate control; direct lift throttle. Simulator supplies calculated hover throttle; no altitude hold.",
    24: "Circles in a thermal; requires enabled, active soaring. No thermal updraft model is supplied.",
    25: "Fixed-wing loiter down to landing altitude, then vertical landing. Requires forward propulsion and landing setup.",
    26: "Automatic fixed-wing landing; unavailable when QuadPlane is enabled and requires a captured takeoff direction.",
}


def session_modes(manifest):
    """Aircraft configuration, not the browser selection, determines command support."""
    return (QUADPLANE_MODES if (manifest or {}).get("parameter_overlay", {}).get("Q_ENABLE", 0)
            else PLANE_MODES)


RESULTS = {0: "Accepted", 1: "Temporarily rejected", 2: "Denied",
           3: "Unsupported", 4: "Failed", 5: "In progress", 6: "Cancelled"}
SEVERITIES = {0: "Emergency", 1: "Alert", 2: "Critical", 3: "Error",
              4: "Warning", 5: "Notice", 6: "Info", 7: "Debug"}


class ModeCommands:
    """One in-flight command, ACK outcome plus authoritative heartbeat confirmation."""

    def __init__(self, directory, emit, modes=None):
        self.directory, self.emit = directory, emit
        self.modes = dict(PLANE_MODES if modes is None else modes)
        self.seen = set()
        self.pending = None

    def update(self, link, message=None, *, ready=False):
        now = time.monotonic()
        if self.pending is not None:
            request = self.pending
            if message and message.get("mavpackettype") == "COMMAND_ACK" and message.get("command") == 176:
                result = message["result"]
                if result not in (0, 5):
                    self._event(request, "rejected", RESULTS.get(result, f"Result {result}"))
                    self.pending = None
                else:
                    self._event(request, "accepted", "Firmware accepted; waiting for mode heartbeat")
            elif message and message.get("mavpackettype") == "HEARTBEAT" and message.get("custom_mode") == request["mode"]:
                self._event(request, "confirmed", f"Mode confirmed: {MODE_NAMES[request['mode']]}")
                self.pending = None
            if self.pending is not None and now-request["sent_at"] > 8:
                self._event(request, "timeout", "No mode confirmation received within 8 seconds")
                self.pending = None
        for path in sorted(self.directory.glob("*.json")):
            if path.name in self.seen:
                continue
            self.seen.add(path.name)
            request = json.loads(path.read_text())
            if (request.get("kind") != "set_mode" or request.get("mode") not in self.modes
                    or now-request.get("created_at", -100) > 10 or not ready):
                self._event(request, "rejected", "Mode command expired, unsupported, or flight not ready")
                continue
            if self.pending is not None:
                self._event(request, "rejected", "Another mode change is awaiting confirmation")
                continue
            # MAV_CMD_DO_SET_MODE, CUSTOM_MODE_ENABLED, Plane custom mode.
            link.mav.command_long_send(link.target_system, link.target_component,
                                       176, 0, 1, request["mode"], 0, 0, 0, 0, 0)
            request["sent_at"] = now
            self.pending = request
            self._event(request, "sent", f"Requested {MODE_NAMES[request['mode']]}")

    def _event(self, request, status, text):
        self.emit(command_id=request.get("id"), mode=request.get("mode"),
                  status=status, text=text)
