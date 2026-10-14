"""Failure injection applied between the firmware's actuator outputs and the physics.

The firmware is not told about the failure; it sees only the aircraft's response,
as on the real aircraft. Times are seconds after release (ground starts: after the
lift motors spin up), matching the wind clock.
"""

from dataclasses import dataclass, replace
import math

# Lift order in the FDM: FR, RR, RL, FL (ArduPilot motors 1, 4, 2, 3).
LIFT_INDEX = {"lift_fr_out": 0, "lift_rr_out": 1, "lift_rl_out": 2, "lift_fl_out": 3}
KINDS = {*LIFT_INDEX, "pusher_out", "battery_low"}

# The firmware reacts to a low pack only with BATT_MONITOR and failsafe actions set.
# Those SITL parameters await approval (parameter proposal rows 7-12), so this stays off.
BATTERY_FAILSAFE_PARAMETERS_APPROVED = False
# State of charge the pack drops to when the battery_low failure starts.
BATTERY_LOW_STATE_OF_CHARGE = .15


@dataclass(frozen=True)
class Failure:
    kind: str
    start_s: float = 20.

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"Unknown failure {self.kind!r}")
        if not math.isfinite(self.start_s) or not 0 <= self.start_s <= 600:
            raise ValueError("Failure start must be 0-600 s after release")
        if self.kind == "battery_low" and not BATTERY_FAILSAFE_PARAMETERS_APPROVED:
            raise ValueError("Low-battery failure awaits approval of the SITL battery parameters")

    @classmethod
    def from_config(cls, config):
        if config is None:
            return None
        if not isinstance(config, dict) or set(config) - {"kind", "start_s"} or "kind" not in config:
            raise ValueError("Failure accepts kind and start_s")
        return cls(str(config["kind"]), float(config.get("start_s", 20.)))

    def as_config(self):
        return {"kind": self.kind, "start_s": self.start_s}

    def active(self, physics):
        return physics.release_time is not None and physics.time - physics.release_time >= self.start_s

    def apply(self, control, physics):
        """Return the control the physics actually receives."""
        if not self.active(physics):
            return control
        if self.kind in LIFT_INDEX:
            commands = control.lift_motor_throttles
            if commands is None:
                raise ValueError("Lift-motor failure requires per-motor lift commands")
            commands = list(commands)
            commands[LIFT_INDEX[self.kind]] = 0.
            return replace(control, lift_motor_throttles=tuple(commands))
        if self.kind == "pusher_out":
            return replace(control, throttle=0.)
        if self.kind == "battery_low" and physics.battery is not None:
            floor = physics.battery.capacity_ah*(1 - BATTERY_LOW_STATE_OF_CHARGE)
            physics.battery.consumed_ah = max(physics.battery.consumed_ah, floor)
        return control
