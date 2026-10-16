"""Main (lift) battery state reported to ArduPilot through the JSON ``battery`` field.

The pack voltage and current are telemetry for ArduPilot's battery monitor and
failsafes. Motor and propeller thrust still use the profile's nominal voltage; the
sag is not fed back into the lift or pusher models. PT1.5's separate pusher pack is
not represented, so pusher current is excluded from this pack.
"""

from dataclasses import dataclass
import math

import numpy as np

# Generic lithium-polymer open-circuit voltage per cell against state of charge.
# Not a measured discharge curve for the installed pack.
SOC_POINTS = (0., .05, .10, .20, .30, .40, .50, .60, .70, .80, .90, 1.)
CELL_OCV_V = (3.30, 3.50, 3.60, 3.68, 3.73, 3.77, 3.80, 3.84, 3.89, 3.95, 4.05, 4.20)

# PT1 hover logs put the pack plus wiring resistance at about 40-55 mOhm.
DEFAULT_RESISTANCE_OHM = .045
# Flight controller, receivers, telemetry and servos; an allowance, not a measurement.
DEFAULT_AVIONICS_CURRENT_A = 1.5


@dataclass
class BatteryModel:
    cells: int
    capacity_ah: float
    resistance_ohm: float = DEFAULT_RESISTANCE_OHM
    avionics_current_a: float = DEFAULT_AVIONICS_CURRENT_A
    consumed_ah: float = 0.

    def __post_init__(self):
        values = (self.capacity_ah, self.resistance_ohm, self.avionics_current_a, self.consumed_ah)
        if (not isinstance(self.cells, int) or self.cells <= 0 or not all(math.isfinite(v) for v in values)
                or self.capacity_ah <= 0 or self.resistance_ohm < 0 or self.avionics_current_a < 0
                or not 0 <= self.consumed_ah <= self.capacity_ah):
            raise ValueError("Battery model needs positive cells/capacity and nonnegative resistance, load and use")

    @classmethod
    def from_params(cls, params):
        return cls(cells=int(params.battery.cells), capacity_ah=float(params.battery.capacity_ah))

    @property
    def state_of_charge(self):
        return 1. - self.consumed_ah/self.capacity_ah

    def open_circuit_voltage(self):
        return self.cells*float(np.interp(self.state_of_charge, SOC_POINTS, CELL_OCV_V))

    def pack_current(self, control, params):
        from ..propulsion import lift_propulsion
        _, _, samples = lift_propulsion(control, params)
        return self.avionics_current_a + sum(sample["current_a"] for sample in samples)

    def reading(self, current_a):
        if not math.isfinite(current_a) or current_a < 0:
            raise ValueError("Battery current must be finite and nonnegative")
        voltage = max(self.open_circuit_voltage() - current_a*self.resistance_ohm, 0.)
        return {"voltage": voltage, "current": float(current_a)}

    def consume(self, current_a, dt):
        if not math.isfinite(dt) or dt < 0:
            raise ValueError("Battery timestep must be finite and nonnegative")
        self.consumed_ah = min(self.capacity_ah, self.consumed_ah + current_a*dt/3600.)
