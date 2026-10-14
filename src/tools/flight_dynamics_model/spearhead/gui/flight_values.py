"""Numeric flight instruments from observed FDM state and firmware outputs."""

import math

VALUE_LABELS = {
    "airspeed": "True airspeed · m/s", "altitude": "Altitude · m", "climb": "Climb rate · m/s",
    "heading": "Heading · °", "roll_rate": "Roll rate · °/s", "pitch_rate": "Pitch rate · °/s",
    "yaw_rate": "Yaw rate · °/s", "pusher": "Pusher command · %",
    "lift_fl": "Lift FL · µs", "lift_fr": "Lift FR · µs", "lift_rl": "Lift RL · µs", "lift_rr": "Lift RR · µs",
}


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
        values.update(airspeed=display(math.sqrt(sum(v*v for v in x[3:6]))), altitude=display(-x[2]),
                      climb=display(-down), heading=display(round(math.degrees(x[8]) % 360, 1) % 360, 1))
        for key, value in zip(("roll_rate", "pitch_rate", "yaw_rate"), x[9:12]):
            values[key] = display(math.degrees(value), 1)
    values["pusher"] = display(100*(row.get("control") or {}).get("throttle"), 1) if (row.get("control") or {}).get("throttle") is not None else "—"
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
