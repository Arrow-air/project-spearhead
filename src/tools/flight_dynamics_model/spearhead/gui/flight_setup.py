"""Validated, independent choices for an ArduPilot flight session."""

from dataclasses import asdict, dataclass
import math

import yaml

from ..aircraft import pusher_options
from ..sitl.backend import CONFIG_ROOT
from ..sitl import failures as failure_model

# Flight Deck flies PT1.5 only; the 25 kg reference design stays in Model Lab (old AeroDB and geometry).
AIRCRAFT = {"pt15": "PT1.5 · 21 kg"}
DOUBLET_LABELS = {"pitch_doublet": "Pitch doublet · elevator ±3°", "roll_doublet": "Roll doublet · aileron ±5°",
                  "yaw_doublet": "Yaw doublet · rudder ±5°", "throttle_doublet": "Throttle doublet · up to ±10 points"}
HOVER_DOUBLET_LABELS = {"pitch_doublet": "Pitch doublet · attitude ±8°", "roll_doublet": "Roll doublet · attitude ±8°",
                        "yaw_doublet": "Yaw doublet · stick ±20%", "throttle_doublet": "Lift throttle doublet · ±10 points"}
POINTS = {
    "pt15": {"ground": "Ground · VTOL takeoff", "hover20": "Hover · 20 m", "hover80": "Hover · 80 m (transition)", "custom_hover": "Custom hover",
                     "trim22": "Forward trim · 22 m/s", "trim24": "Forward trim · 24 m/s",
                     "trim26": "Forward trim · 26 m/s", "custom_trim": "Custom forward trim"},
    "reference_design": {"trim22": "Forward trim · 22 m/s", "trim24": "Forward trim · 24 m/s",
                         "trim26": "Forward trim · 26 m/s", "trim28": "Forward trim · 28 m/s", "custom_trim": "Custom forward trim"},
}
MANEUVERS = {
    "pt15": {"interactive": "Interactive flight", "hold": "Hover hold", "transition": "Hover → forward flight test", "takeoff_hover_land": "Takeoff → hover → land"},
    "reference_design": {"interactive": "Interactive flight", "hold": "Straight & level hold",
                         "steps": "Pitch / roll steps", "plus": "+0.5° attitude disturbance",
                         "minus": "−0.5° attitude disturbance", "elevator": "Elevator doublet · ±1°",
                         "manual": "Surface direction checks"},
}
MANEUVERS["reference_design"].update(DOUBLET_LABELS)
MANEUVERS["pt15"].update(HOVER_DOUBLET_LABELS)
MIN_DURATION = {"steps": 118, "elevator": 45, "transition": 60, "back_transition": 60,
                **{key: 45 for key in DOUBLET_LABELS}}
# Body x from the CFD origin (wing root quarter-chord), forward positive. The 14 October
# drawing puts the lift-rotor centroid at x = -37.5 mm; 22 m/s margins are October AeroDB
# estimates at z = +65 mm and move with any aerodynamic correction.
TRANSITION_CGS = {
    "target40": "PT1.5 target · x −40 mm (SM ≈ 6 %, default)",
    "pt1_log": "PT1 logged · x −47.5 mm (SM ≈ 5 %)",
    "aft60": "Aft · x −60 mm (SM ≈ 2 %)",
    "forward25": "Forward · x −25 mm (SM ≈ 10 %)",
    "forward0": "Forward · x 0 mm (SM ≈ 15 %, elevator-limited)",
}
TRANSITION_CG_X = {"target40": -.040, "pt1_log": -.0475, "aft60": -.060, "forward25": -.025, "forward0": 0.}
# Build-up estimate of the CG height below the origin; unmeasured.
TRANSITION_CG_Z = .065
FAILURES = {
    "none": "No failure",
    "lift_fr_out": "Lift motor FR (motor 1) out",
    "lift_rr_out": "Lift motor RR (motor 4) out",
    "lift_rl_out": "Lift motor RL (motor 2) out",
    "lift_fl_out": "Lift motor FL (motor 3) out",
    "pusher_out": "Pusher cut",
    "battery_low": "Low battery (awaits BATT parameter approval)",
}
# Pusher motor options of the PT1.5 profile (MAD V8013 150KV or KDE5215XF-220, both with the 20x12 propeller).
PUSHERS = pusher_options("pt15")
# A gust reaches its peak 1.5 s after it starts; same direction as the steady wind.
GUST_START_S, GUST_DURATION_S = 20., 3.
TRANSITION_PARAMETERS = {
    "first_transition_2026_10_09": "9 October proposal · frozen factory rate gains",
}


def forward_start(initial_point):
    return initial_point.startswith("trim") or initial_point == "custom_trim"


def maneuver_options(aircraft, initial_point):
    if initial_point == "ground":
        return {"takeoff_hover_land": MANEUVERS[aircraft]["takeoff_hover_land"]}
    if aircraft == "pt15" and forward_start(initial_point):
        return {**{key: MANEUVERS["reference_design"][key] for key in ("interactive", "hold")},
                "back_transition": "Forward → hover transition", **DOUBLET_LABELS}
    return {key: value for key, value in MANEUVERS[aircraft].items() if key != "takeoff_hover_land"}


@dataclass(frozen=True)
class FlightSetup:
    aircraft: str = "pt15"
    initial_point: str = "hover20"
    maneuver: str = "interactive"
    duration_s: float = 60
    altitude_m: float = 20
    speed_m_s: float = 0
    cg_case: str = "target40"
    transition_parameters: str = "first_transition_2026_10_09"
    wind_speed_m_s: float = 0
    wind_from_deg: float = 0
    gust_m_s: float = 0
    failure: str = "none"
    failure_time_s: float = 20
    pusher: str = "mad"

    def validate(self):
        if self.aircraft not in AIRCRAFT or self.initial_point not in POINTS[self.aircraft]:
            raise ValueError("Choose an initial point supported by the selected aircraft")
        if self.maneuver not in maneuver_options(self.aircraft, self.initial_point):
            raise ValueError("Choose a maneuver supported by the selected aircraft")
        if self.cg_case not in TRANSITION_CGS or self.transition_parameters not in TRANSITION_PARAMETERS:
            raise ValueError("Choose a supported transition CG and parameter set")
        conditions = (self.wind_speed_m_s, self.wind_from_deg, self.gust_m_s, self.failure_time_s)
        if not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in conditions):
            raise ValueError("Wind, gust and failure time must be finite numbers")
        if not (0 <= self.wind_speed_m_s <= 15 and 0 <= self.gust_m_s <= 10 and 0 <= self.wind_from_deg < 360):
            raise ValueError("Wind 0–15 m/s from 0–359°, gust 0–10 m/s")
        if self.failure not in FAILURES:
            raise ValueError("Choose a supported failure")
        if self.pusher not in PUSHERS:
            raise ValueError("Choose a supported pusher motor")
        if self.failure != "none":
            if self.aircraft != "pt15":
                raise ValueError("Failure injection is available for PT1.5")
            if self.failure == "battery_low" and not failure_model.BATTERY_FAILSAFE_PARAMETERS_APPROVED:
                raise ValueError("Low-battery failure awaits approval of the SITL battery parameters")
            if self.failure.startswith("lift_") and self.is_forward and self.maneuver != "back_transition":
                raise ValueError("Lift-motor failures need a hover, ground or forward → hover start")
            if self.failure == "pusher_out" and not (self.is_forward or self.maneuver == "transition"):
                raise ValueError("Pusher failures need a forward start or a hover → forward transition")
            if not 0 <= self.failure_time_s < self.duration_s:
                raise ValueError("Failure time must fall inside the flight")
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
        hover = not self.is_forward
        if not ((3 if ground else 6) <= self.altitude_m <= (30 if ground else 35 if hover and self.maneuver == "hold" else 1000 if not hover else 110)):
            raise ValueError("Initial altitude is outside this runner's diagnostic envelope")
        if hover and self.speed_m_s != 0:
            raise ValueError("This PT1.5 hover runner starts at zero airspeed; forward starts need a separately validated trim setup")
        if self.maneuver in ("transition", "back_transition") and not 40 <= self.altitude_m <= 110:
            raise ValueError("Choose Hover · 80 m or Custom hover at 40–110 m for the transition diagnostic")
        upper_speed = 26 if self.aircraft == "pt15" else 30
        if not hover and not 22 <= self.speed_m_s <= upper_speed:
            raise ValueError(f"Forward-flight initial airspeed must be between 22 and {upper_speed} m/s")
        # A named point has a precise meaning; arbitrary numeric edits use Custom.
        if self.initial_point.startswith("hover") and self.altitude_m != int(self.initial_point[5:]):
            raise ValueError("Select Custom hover to edit the initial altitude")
        if self.initial_point.startswith("trim") and self.speed_m_s != int(self.initial_point[4:]):
            raise ValueError("Select Custom forward trim to edit the initial airspeed")
        return self

    @property
    def is_forward(self):
        return forward_start(self.initial_point)

    @property
    def scenario(self):
        if self.aircraft == "pt15" and not self.is_forward:
            if self.initial_point == "ground":
                return "quad_ground"
            if self.maneuver == "transition":
                return "quad_transition"
            if self.maneuver in HOVER_DOUBLET_LABELS:
                return "quad_"+self.maneuver
            return "quad_interactive" if self.maneuver == "interactive" else "quad_hover"
        return "fbwa" if self.maneuver == "hold" else self.maneuver

    def config(self):
        self.validate()
        template = ("quad_ground.yaml" if self.initial_point == "ground" else
                    "quad_transition.yaml" if self.maneuver == "transition" else
                    "quad_transition.yaml" if self.aircraft == "pt15" and self.is_forward else
                    "quad_hover.yaml" if self.aircraft == "pt15" else "fixed_wing_control.yaml")
        data = yaml.safe_load((CONFIG_ROOT / template).read_text())
        data.update(aircraft=self.aircraft, altitude_m=float(self.altitude_m), speed_m_s=float(self.speed_m_s))
        if self.aircraft == "pt15" and self.is_forward:
            data["flight_regime"] = "fixed_wing"
            # This air start uses the fixed-wing runner; lift remains off.
            if self.maneuver != "back_transition":
                data["channels"].pop("lift")
        if self.wind_speed_m_s or self.gust_m_s:
            data["wind"] = {"speed_m_s": float(self.wind_speed_m_s), "from_deg": float(self.wind_from_deg),
                            "gusts": [{"start_s": GUST_START_S, "duration_s": GUST_DURATION_S,
                                       "speed_m_s": float(self.gust_m_s), "from_deg": float(self.wind_from_deg)}]
                            if self.gust_m_s else []}
        if self.failure != "none":
            data["failure"] = {"kind": self.failure, "start_s": float(self.failure_time_s)}
        if self.aircraft == "pt15":
            data["pusher"] = self.pusher
            data["cg_body_xyz_m"] = [TRANSITION_CG_X[self.cg_case], 0., TRANSITION_CG_Z]
            if not self.is_forward:
                data["hover_trim"] = "independent_rotors"
            else:
                data.pop("hover_trim", None)
        return data

    def as_dict(self):
        return asdict(self)
