"""Configurable provisional surface travel and independent lift motor mapping."""

from dataclasses import dataclass
import math

from ..controls import Control


@dataclass(frozen=True)
class Channel:
    number: int
    minimum: int = 1000
    trim: int = 1500
    maximum: int = 2000
    reverse: bool = False
    travel_deg: float = 5.0

    def __post_init__(self):
        if not 1 <= self.number <= 32 or not 800 <= self.minimum < self.trim < self.maximum <= 2200:
            raise ValueError("Invalid PWM channel or limits")
        if not math.isfinite(self.travel_deg) or not 0 < self.travel_deg <= 45:
            raise ValueError("SITL surface travel must be in (0, 45] degrees")

    def surface(self, pwm):
        value = pwm[self.number - 1]
        if value == 0:
            value = self.trim
        span = self.maximum - self.trim if value >= self.trim else self.trim - self.minimum
        normalized = max(-1., min(1., (value - self.trim) / span))
        return math.radians(self.travel_deg) * normalized * (-1 if self.reverse else 1)

    def motor(self, pwm):
        return max(0., min(1., (pwm[self.number - 1] - self.minimum) / (self.maximum - self.minimum)))


@dataclass(frozen=True)
class ActuatorMap:
    aileron: Channel = Channel(1)
    elevator: Channel = Channel(2)
    pusher: Channel = Channel(3)
    rudder: Channel = Channel(4)
    # AP Quad-X motors 1,4,2,3 -> FDM FR,RR,RL,FL.
    lift: tuple[Channel, ...] = ()

    def __post_init__(self):
        channels = [self.aileron, self.elevator, self.pusher, self.rudder, *self.lift]
        if len({c.number for c in channels}) != len(channels):
            raise ValueError("Actuator channels must be unique")
        if len(self.lift) not in (0, 4):
            raise ValueError("Configure either zero or four lift motors")

    def control(self, pwm):
        channels = [self.aileron, self.elevator, self.pusher, self.rudder, *self.lift]
        if max(c.number for c in channels) > len(pwm):
            raise ValueError("Packet does not contain configured actuator channel")
        return Control(self.elevator.surface(pwm), self.aileron.surface(pwm),
                       self.rudder.surface(pwm), self.pusher.motor(pwm),
                       lift_motor_throttles=tuple(c.motor(pwm) for c in self.lift) or None)
