"""Steady wind with discrete 1-cosine gusts, in NED, for air-relative aerodynamics.

Directions follow the meteorological convention: ``from_deg`` is where the wind
blows from (0 = from north). Times are seconds after the SITL release, so the
stationary pitot calibration before release always sees still air. The field is
spatially uniform; there is no turbulence spectrum, wind shear or rotor wake.
"""

from dataclasses import dataclass
import math

import numpy as np

MAX_WIND_M_S = 20.


def _horizontal_ned(speed_m_s, from_deg):
    heading = math.radians(from_deg)
    # Wind blowing from north moves air toward the south (negative north).
    return np.array([-speed_m_s*math.cos(heading), -speed_m_s*math.sin(heading), 0.])


@dataclass(frozen=True)
class Gust:
    start_s: float
    duration_s: float
    speed_m_s: float
    from_deg: float
    vertical_m_s: float = 0.

    def __post_init__(self):
        values = (self.start_s, self.duration_s, self.speed_m_s, self.from_deg, self.vertical_m_s)
        if (not all(math.isfinite(v) for v in values) or self.start_s < 0 or self.duration_s <= 0
                or not 0 <= self.speed_m_s <= MAX_WIND_M_S or abs(self.vertical_m_s) > MAX_WIND_M_S):
            raise ValueError("Gust needs start >= 0, duration > 0 and speeds up to 20 m/s")

    def ned(self, t):
        if not self.start_s <= t <= self.start_s + self.duration_s:
            return np.zeros(3)
        shape = .5*(1 - math.cos(2*math.pi*(t - self.start_s)/self.duration_s))
        return shape*(_horizontal_ned(self.speed_m_s, self.from_deg) + np.array([0., 0., self.vertical_m_s]))


@dataclass(frozen=True)
class WindField:
    speed_m_s: float = 0.
    from_deg: float = 0.
    gusts: tuple = ()

    def __post_init__(self):
        if (not math.isfinite(self.speed_m_s) or not math.isfinite(self.from_deg)
                or not 0 <= self.speed_m_s <= MAX_WIND_M_S):
            raise ValueError("Steady wind must be 0-20 m/s with a finite direction")
        if not all(isinstance(gust, Gust) for gust in self.gusts):
            raise ValueError("Gusts must be Gust instances")

    @classmethod
    def from_config(cls, config):
        if config is None:
            return None
        if not isinstance(config, dict) or set(config) - {"speed_m_s", "from_deg", "gusts"}:
            raise ValueError("Wind accepts speed_m_s, from_deg and gusts")
        gusts = config.get("gusts", [])
        if not isinstance(gusts, list):
            raise ValueError("Wind gusts must be a list")
        return cls(float(config.get("speed_m_s", 0.)), float(config.get("from_deg", 0.)),
                   tuple(Gust(**{key: float(value) for key, value in gust.items()}) for gust in gusts))

    def as_config(self):
        return {"speed_m_s": self.speed_m_s, "from_deg": self.from_deg,
                "gusts": [{"start_s": g.start_s, "duration_s": g.duration_s, "speed_m_s": g.speed_m_s,
                           "from_deg": g.from_deg, "vertical_m_s": g.vertical_m_s} for g in self.gusts]}

    @property
    def calm(self):
        return self.speed_m_s == 0 and not self.gusts

    def ned(self, t):
        wind = _horizontal_ned(self.speed_m_s, self.from_deg)
        for gust in self.gusts:
            wind = wind + gust.ned(t)
        return wind
