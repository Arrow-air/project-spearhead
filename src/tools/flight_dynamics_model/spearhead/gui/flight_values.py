"""Numeric flight instruments from observed FDM state and firmware outputs."""

import math

VALUE_LABELS = {
    "airspeed": "True airspeed · m/s", "altitude": "Altitude · m", "climb": "Climb rate · m/s",
    "heading": "Heading · °", "roll_rate": "Roll rate · °/s", "pitch_rate": "Pitch rate · °/s",
    "yaw_rate": "Yaw rate · °/s", "pusher": "Pusher command · %",
    "lift_fl": "Lift FL · µs", "lift_fr": "Lift FR · µs", "lift_rl": "Lift RL · µs", "lift_rr": "Lift RR · µs",
    "elevator": "Elevator · °", "pusher_current": "Pusher current · A",
    "battery": "Lift pack · V / A", "battery_used": "Lift pack used · mAh", "wind": "Wind · m/s",
}
# Measured elevator coverage of the 11 October control delivery.
MEASURED_ELEVATOR_DEG = 15.
# Per-cell voltage under load treated as low for the 12S lift pack.
LOW_CELL_V = 3.5


def numeric_flight_values(view):
    values = dict.fromkeys(VALUE_LABELS, "—")
    row = view.get("truth")
    if row is None:
        return values
    x = row.get("state") or []

    def display(value, digits=2):
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            return "—"
        if abs(value) < .5 * 10**-digits:
            value = 0.
        return f"{value:.{digits}f}"

    if len(x) >= 12 and all(math.isfinite(v) for v in x):
        phi, theta = x[6:8]
        # State velocities are body-frame components; climb is earth vertical.
        down = -math.sin(theta)*x[3] + math.sin(phi)*math.cos(theta)*x[4] + math.cos(phi)*math.cos(theta)*x[5]
        # Air-relative speed when the run records it (wind runs); otherwise inertial speed.
        tas = row.get("true_airspeed_m_s")
        if not isinstance(tas, (int, float)):
            tas = math.sqrt(sum(v*v for v in x[3:6]))
        values.update(airspeed=display(tas), altitude=display(-x[2]),
                      climb=display(-down), heading=display(round(math.degrees(x[8]) % 360, 1) % 360, 1))
        for key, value in zip(("roll_rate", "pitch_rate", "yaw_rate"), x[9:12]):
            values[key] = display(math.degrees(value), 1)
    values["pusher"] = display(100*(row.get("control") or {}).get("throttle"), 1) if (row.get("control") or {}).get("throttle") is not None else "—"
    control = row.get("control") or {}
    if isinstance(control.get("de"), (int, float)):
        values["elevator"] = display(math.degrees(control["de"]), 1)
    pusher = row.get("pusher") or {}
    if isinstance(pusher.get("current_a"), (int, float)):
        values["pusher_current"] = display(pusher["current_a"], 1)
    battery = (row.get("sensor") or {}).get("battery") or {}
    if isinstance(battery.get("voltage"), (int, float)):
        values["battery"] = f"{display(battery['voltage'], 1)} / {display(battery.get('current'), 0)}"
    if isinstance(row.get("battery_consumed_mah"), (int, float)):
        values["battery_used"] = display(row["battery_consumed_mah"], 0)
    wind = row.get("wind_ned_m_s")
    values["wind"] = display(math.hypot(wind[0], wind[1]), 1) if wind else "0.0" if "wind_ned_m_s" in row else "—"
    pwm = row.get("pwm") or []
    for key, index in (("lift_fr", 4), ("lift_rr", 7), ("lift_rl", 5), ("lift_fl", 6)):
        if len(pwm) > index and pwm[index] != 0:
            values[key] = display(pwm[index], 0)
    return values


def flight_values_status(view):
    if view.get("truth") is None:
        return "WAITING FOR FLIGHT DATA", False
    if view.get("active"):
        if view.get("stale_truth") or view.get("stale_mavlink"):
            return "STALE · LAST RECEIVED VALUES", True
        if view["truth"].get("hold"):
            return "STARTUP · AIRCRAFT HELD", False
        return "LIVE", False
    return "RECORDED · LAST SAMPLE", False


def live_alerts(row, cells=12):
    """Short warnings a test engineer should see while flying, from one frame record."""
    if row is None or row.get("hold"):
        return []
    alerts = []
    if row.get("failure_active"):
        alerts.append("Failure injected")
    control = row.get("control") or {}
    if isinstance(control.get("de"), (int, float)) and abs(math.degrees(control["de"])) > MEASURED_ELEVATOR_DEG:
        alerts.append(f"Elevator beyond measured ±{MEASURED_ELEVATOR_DEG:g}° CFD data")
    elif row.get("control_coverage_exceeded"):
        alerts.append("Outside measured aerodynamic data")
    pusher = row.get("pusher") or {}
    if pusher.get("above_sustainable_current"):
        alerts.append(f"Pusher above {pusher.get('sustainable_current_a', 24):g} A sustainable current")
    voltage = ((row.get("sensor") or {}).get("battery") or {}).get("voltage")
    if isinstance(voltage, (int, float)) and voltage < LOW_CELL_V*cells:
        alerts.append(f"Lift pack below {LOW_CELL_V:g} V/cell")
    return alerts
