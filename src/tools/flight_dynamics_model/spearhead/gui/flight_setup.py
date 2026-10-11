"""Validated, independent choices for an ArduPilot flight session."""

from dataclasses import asdict, dataclass
import math

import yaml

from ..sitl.backend import CONFIG_ROOT

AIRCRAFT = {"pt1_as_built": "PT1 · 21 kg", "reference_design": "Spearhead · 25 kg"}
POINTS = {
    "pt1_as_built": {"ground": "Ground · VTOL takeoff", "hover20": "Hover · 20 m", "hover30": "Hover · 30 m", "custom_hover": "Custom hover"},
    "reference_design": {"trim22": "Forward trim · 22 m/s", "trim24": "Forward trim · 24 m/s",
                         "trim26": "Forward trim · 26 m/s", "trim28": "Forward trim · 28 m/s", "custom_trim": "Custom forward trim"},
}
MANEUVERS = {
    "pt1_as_built": {"interactive": "Interactive flight", "hold": "Hover hold", "takeoff_hover_land": "Takeoff → hover → land"},
    "reference_design": {"interactive": "Interactive flight", "hold": "Straight & level hold",
                         "steps": "Pitch / roll steps", "plus": "+0.5° attitude disturbance",
                         "minus": "−0.5° attitude disturbance", "elevator": "Elevator doublet · ±1°",
                         "manual": "Surface direction checks"},
}
MIN_DURATION = {"steps": 118, "elevator": 45}


@dataclass(frozen=True)
class FlightSetup:
    aircraft: str = "pt1_as_built"
    initial_point: str = "hover20"
    maneuver: str = "interactive"
    duration_s: float = 60
    altitude_m: float = 20
    speed_m_s: float = 0

    def validate(self):
        if self.aircraft not in AIRCRAFT or self.initial_point not in POINTS[self.aircraft]:
            raise ValueError("Choose an initial point supported by the selected aircraft")
        if self.maneuver not in MANEUVERS[self.aircraft]:
            raise ValueError("Choose a maneuver supported by the selected aircraft")
        ground = self.initial_point == "ground"
        if ground != (self.maneuver == "takeoff_hover_land"):
            raise ValueError("Ground starts use Takeoff → hover → land; air starts use a flight maneuver")
        if not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
                   for v in (self.duration_s, self.altitude_m, self.speed_m_s)):
            raise ValueError("Duration, altitude and airspeed must be finite numbers")
        if ground and self.duration_s != int(self.duration_s):
            raise ValueError("VTOL mission hover time must be a whole number of seconds")
        minimum = MIN_DURATION.get(self.maneuver, 15)
        if not minimum <= self.duration_s <= 600:
            raise ValueError(f"This maneuver requires {minimum}–600 seconds of free flight")
        hover = self.aircraft == "pt1_as_built"
        if not ((3 if ground else 6) <= self.altitude_m <= (30 if ground else 35 if hover and self.maneuver == "hold" else 1000 if not hover else 110)):
            raise ValueError("Initial altitude is outside this runner's diagnostic envelope")
        if hover and self.speed_m_s != 0:
            raise ValueError("This PT1 hover runner starts at zero airspeed; forward starts need a separately validated trim setup")
        if not hover and not 22 <= self.speed_m_s <= 30:
            raise ValueError("Forward-flight initial airspeed must be between 22 and 30 m/s")
        # A named point has a precise meaning; arbitrary numeric edits use Custom.
        if self.initial_point.startswith("hover") and self.altitude_m != int(self.initial_point[5:]):
            raise ValueError("Select Custom hover to edit the initial altitude")
        if self.initial_point.startswith("trim") and self.speed_m_s != int(self.initial_point[4:]):
            raise ValueError("Select Custom forward trim to edit the initial airspeed")
        return self

    @property
    def scenario(self):
        if self.aircraft == "pt1_as_built":
            if self.initial_point == "ground":
                return "quad_ground"
            return "quad_interactive" if self.maneuver == "interactive" else "quad_hover"
        return "fbwa" if self.maneuver == "hold" else self.maneuver

    def config(self):
        self.validate()
        template = "quad_ground.yaml" if self.initial_point == "ground" else "quad_hover.yaml" if self.aircraft == "pt1_as_built" else "fixed_wing_control.yaml"
        data = yaml.safe_load((CONFIG_ROOT / template).read_text())
        data.update(aircraft=self.aircraft, altitude_m=float(self.altitude_m), speed_m_s=float(self.speed_m_s))
        return data

    def as_dict(self):
        return asdict(self)
