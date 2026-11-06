"""Ground-station readouts from firmware telemetry; absent data stays unknown."""

from ..sitl.ground_control import SEVERITIES


class StatusTexts:
    """Reassemble MAVLink 2 chunks, retaining severity and original timestamp."""

    def __init__(self, feed):
        self.feed, self.parts = feed, {}

    def add(self, message, timestamp):
        text = message.get("text", "")
        if isinstance(text, bytes):
            text = text.decode(errors="replace")
        ident, sequence = message.get("id", 0), message.get("chunk_seq", 0)
        severity = message.get("severity", 6)
        if not ident:
            self._emit(text.rstrip("\0"), severity, timestamp)
            return
        if sequence == 0:
            self.parts[ident] = (text, severity, timestamp, 0)
        elif ident in self.parts and self.parts[ident][3]+1 == sequence:
            before, severity, timestamp, _ = self.parts[ident]
            self.parts[ident] = (before+text, severity, timestamp, sequence)
        else:
            self._emit("[incomplete message] " + text.rstrip("\0"), severity, timestamp)
            self.parts.pop(ident, None)
            return
        if len(text.encode()) < 50 or "\0" in text:
            joined, severity, timestamp, _ = self.parts.pop(ident)
            self._emit(joined.rstrip("\0"), severity, timestamp)
        # Bound incomplete messages even if a final chunk is lost.
        while len(self.parts) > 20:
            self.parts.pop(next(iter(self.parts)))

    def _emit(self, text, severity, timestamp):
        self.feed.append({"time_s": timestamp, "severity": severity,
                          "level": SEVERITIES.get(severity, "Unknown"),
                          "source": "ArduPilot", "text": text})


def firmware_health(messages, fdm_time, *, live=False):
    def recent(kind):
        row = messages.get(kind)
        if row is None or (live and (fdm_time is None or fdm_time-row["fdm_time_s"] > 3)):
            return {}
        return row["message"]

    gps, ekf, system = recent("GPS_RAW_INT"), recent("EKF_STATUS_REPORT"), recent("SYS_STATUS")
    flags = ekf.get("flags")
    ekf_ok = None if flags is None else flags & 63 == 63 and not flags & (1024 | 32768)
    present = system.get("onboard_control_sensors_present", 0)
    enabled = system.get("onboard_control_sensors_enabled", 0)
    healthy = system.get("onboard_control_sensors_health", 0)
    failed = present & enabled & ~healthy & 0x7fffffff
    sensor_names = {1: "Gyro", 2: "Accelerometer", 4: "Compass", 8: "Barometer",
                    16: "Airspeed", 32: "GPS", 64: "Optical flow", 128: "Vision position",
                    256: "Rangefinder position", 512: "External position", 1024: "Rate control", 2048: "Attitude control",
                    4096: "Yaw control", 8192: "Altitude control", 16384: "Position control",
                    32768: "Motor outputs", 65536: "RC receiver", 131072: "Second gyro",
                    262144: "Second accelerometer", 524288: "Second compass", 1048576: "Geofence",
                    2097152: "AHRS", 4194304: "Terrain data", 8388608: "Reverse motor",
                    16777216: "Logging", 33554432: "Battery", 67108864: "Proximity",
                    134217728: "Satellite link", 268435456: "Pre-arm checks",
                    536870912: "Collision avoidance", 1073741824: "Propulsion"}
    problems = [name for bit, name in sensor_names.items() if failed & bit]
    unknown = failed & ~sum(sensor_names)
    if unknown:
        problems.append(f"Other sensor flags: 0x{unknown:x}")
    voltage, current, remaining = (system.get(k) for k in ("voltage_battery", "current_battery", "battery_remaining"))
    satellites = gps.get("satellites_visible")
    return {"gps_fix": gps.get("fix_type"), "satellites": None if satellites == 255 else satellites,
            "ekf_ok": ekf_ok, "sensor_faults": problems, "system_known": bool(system),
            "battery_v": None if voltage in (None, 65535, 0) else voltage/1000,
            "battery_a": None if current in (None, -1) else current/100,
            "battery_pct": None if remaining in (None, -1) else remaining}
