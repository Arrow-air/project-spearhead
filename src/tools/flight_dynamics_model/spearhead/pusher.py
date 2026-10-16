"""Provisional electric pusher estimate: motor equilibrium + APC propeller data.

APC coefficients are analytical, not a measured installed propulsion curve.
The equivalent DC resistance includes a fit to the manufacturer's motor/ESC
test sheet; it is not the published phase-to-phase winding resistance.
"""

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class PropellerMap:
    rpms: np.ndarray
    sections: tuple[np.ndarray, ...]  # J, Ct, Cp

    def coefficients(self, rpm: float, advance_ratio: float):
        high = int(np.clip(np.searchsorted(self.rpms, rpm), 1, len(self.rpms)-1))
        low = high-1
        fraction = float(np.clip((rpm-self.rpms[low])/(self.rpms[high]-self.rpms[low]), 0, 1))
        coefficients = []
        outside = rpm < self.rpms[0] or rpm > self.rpms[-1]
        for section in (self.sections[low], self.sections[high]):
            outside |= advance_ratio > section[-1, 0] or advance_ratio < section[0, 0]
            # Positive-thrust envelope only: beyond the last J, set thrust to
            # zero and hold endpoint power. Windmilling/braking is unmodelled.
            ct = np.interp(advance_ratio, section[:, 0], section[:, 1], right=0)
            cp = np.interp(advance_ratio, section[:, 0], section[:, 2])
            coefficients.append((max(0., ct), max(0., cp)))
        ct, cp = (1-fraction)*np.array(coefficients[0])+fraction*np.array(coefficients[1])
        return float(ct), float(cp), bool(outside)


@lru_cache(maxsize=8)
def _load_map(path: str, modified_ns: int, size: int) -> PropellerMap:
    sections = {}
    rpm = None
    for line in Path(path).read_text().splitlines():
        match = re.search(r"PROP RPM\s*=\s*(\d+)", line)
        if match:
            rpm = int(match[1])
            if rpm in sections:
                raise ValueError("Duplicate APC RPM section")
            sections[rpm] = []
        elif rpm is not None:
            try:
                row = [float(value) for value in line.split()]
            except ValueError:
                continue
            if len(row) == 15:
                sections[rpm].append([row[1], row[3], row[4]])
    if len(sections) < 2:
        raise ValueError("APC map needs at least two RPM sections")
    arrays = tuple(np.asarray(sections[key]) for key in sorted(sections))
    for section in arrays:
        if (section.ndim != 2 or section.shape[0] < 2 or section.shape[1] != 3
                or not np.all(np.isfinite(section)) or section[0, 0] != 0
                or np.any(np.diff(section[:, 0]) <= 0)):
            raise ValueError("Invalid APC advance-ratio grid")
    return PropellerMap(np.asarray(sorted(sections), dtype=float), arrays)


def load_propeller_map(path: Path) -> PropellerMap:
    path = Path(path).resolve()
    info = path.stat()
    return _load_map(str(path), info.st_mtime_ns, info.st_size)


def estimate_pusher(system, voltage_v: float, throttle: float, axial_speed_m_s: float, rho: float):
    """Solve steady RPM at commanded duty, reporting battery-side current.

    Maximum current, input power and RPM are model limits, not an emulation of
    ESC protection or ArduPilot commands. The sustainable-current rating is a
    reported operating limit; no invented thermal time constant is applied.
    """
    if not np.all(np.isfinite([voltage_v, throttle, axial_speed_m_s, rho])) or rho <= 0:
        raise ValueError("Pusher inputs must be finite and density positive")
    if not system.voltage_min_v <= voltage_v <= system.voltage_max_v:
        raise ValueError("Pusher voltage is outside the configured 12S envelope")
    # Duty is the voltage fraction applied to the motor; some ESCs (KDE) map throttle to duty nonlinearly.
    requested = system.duty(throttle, voltage_v)
    empty = dict(thrust_n=0., current_a=0., power_w=0., shaft_power_w=0., rpm=0.,
                 torque_nm=0., effective_duty=0., limited=False,
                 above_sustainable_current=False, outside_propeller_data=False)
    if requested == 0:
        return empty
    propeller = load_propeller_map(system.performance_path)
    diameter = system.propeller["diameter_in"]*.0254
    kv = system.motor["kv_rpm_per_v"]
    kt = 60/(2*np.pi*kv)
    resistance = system.equivalent_resistance_ohm
    idle_current = system.motor["no_load_current_a"]
    speed = max(0., axial_speed_m_s)

    def operating_point(rpm):
        n = rpm/60
        ct, cp, outside = propeller.coefficients(rpm, speed/(n*diameter))
        shaft_power = cp*rho*n**3*diameter**5
        torque = shaft_power/(2*np.pi*n)
        winding_current = idle_current+torque/kt
        motor_voltage = rpm/kv+winding_current*resistance
        duty = motor_voltage/voltage_v
        electrical_power = motor_voltage*winding_current/system.esc_efficiency
        return dict(thrust_n=float(ct*rho*n**2*diameter**4),
                    current_a=float(electrical_power/voltage_v), power_w=float(electrical_power),
                    shaft_power_w=float(shaft_power), rpm=float(rpm), torque_nm=float(torque),
                    effective_duty=float(duty), outside_propeller_data=bool(outside or axial_speed_m_s < 0))

    # A monotone envelope prevents extrapolation past the supported RPM range.
    # Below startup voltage there is no steady rotating solution.
    if requested*voltage_v <= idle_current*resistance:
        return empty
    maximum_rpm = min(system.propeller["max_rpm"], propeller.rpms[-1], requested*voltage_v*kv)
    def constraint(rpm):
        point = operating_point(rpm)
        return max(point["effective_duty"]/requested,
                   point["current_a"]/system.max_current_a,
                   point["power_w"]/system.max_input_power_w)-1
    rpm = maximum_rpm if constraint(maximum_rpm) <= 0 else brentq(constraint, 1e-6, maximum_rpm, xtol=1e-5)
    result = operating_point(rpm)
    result.update(limited=result["effective_duty"] < requested-1e-5,
                  above_sustainable_current=result["current_a"] > system.sustainable_current_a*(1+1e-9))
    return result
