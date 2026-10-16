"""Partial measured control increments; angles in source CSV are degrees.

Each surface is differenced against its own zero-deflection case. This avoids
adding a second base load when the independently exported databases differ.
Alpha/beta hold endpoints; deflections optionally extrapolate to assumed travel.
Surface interactions and aerodynamic rate derivatives are not measured here.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .aerodb import BODY_COEFFICIENT_COLUMNS
from .params import AeroDBReference

SURFACES = {"aileron": "da", "elevator": "de", "rudder": "dr"}
LOADS = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")


@dataclass(frozen=True)
class SurfaceIncrement:
    alpha_deg: np.ndarray
    deflection_deg: np.ndarray
    coefficients: np.ndarray

    def at(self, alpha_deg: float, deflection_deg: float, *, extrapolate=False) -> tuple[np.ndarray, bool]:
        alpha_used = float(np.clip(alpha_deg, self.alpha_deg[0], self.alpha_deg[-1]))
        deflection_used = deflection_deg if extrapolate else float(np.clip(deflection_deg, self.deflection_deg[0], self.deflection_deg[-1]))
        held = alpha_used != alpha_deg or deflection_used != deflection_deg
        return self._interpolator([[alpha_used, deflection_used]])[0], held

    def __post_init__(self):
        object.__setattr__(self, "_interpolator", RegularGridInterpolator(
            (self.alpha_deg, self.deflection_deg), self.coefficients, bounds_error=False, fill_value=None,
        ))


@dataclass(frozen=True)
class DeflectionDatabase:
    surfaces: dict[str, SurfaceIncrement]
    beta_deg: float
    reference: AeroDBReference
    case_count: int

    def increments_at(self, alpha_deg: float, beta_deg: float, control, *,
                      deflection_policy="hold", limits_deg=(5., 5., 5.), velocity_m_s=None):
        if deflection_policy not in ("hold", "linear"):
            raise ValueError("Deflection policy must be hold or linear")
        if len(limits_deg) != 3 or not np.isfinite(limits_deg).all() or not all(0 < v <= 45 for v in limits_deg):
            raise ValueError("Surface limits must contain three positive angles up to 45 degrees")
        limits = dict(zip(("elevator", "aileron", "rudder"), limits_deg))
        total = np.zeros(6)
        held = []
        extrapolated = []
        travel_limited = []
        active = []
        for name, channel in SURFACES.items():
            value_deg = float(np.rad2deg(getattr(control, channel)))
            if not np.isfinite(value_deg):
                raise ValueError("Control deflections must be finite")
            # A zero command has an exact zero increment, even outside coverage.
            if value_deg == 0:
                continue
            active.append(name)
            surface = self.surfaces[name]
            used = float(np.clip(value_deg, -limits[name], limits[name])) if deflection_policy == "linear" else value_deg
            if used != value_deg:
                travel_limited.append(name)
            if deflection_policy == "linear" and not surface.deflection_deg[0] <= used <= surface.deflection_deg[-1]:
                extrapolated.append(name)
            increment, outside = surface.at(alpha_deg, used, extrapolate=deflection_policy == "linear")
            total += increment
            if outside:
                held.append(name)
        beta_held = bool(active) and abs(beta_deg - self.beta_deg) > 1e-9
        return dict(zip(BODY_COEFFICIENT_COLUMNS, total.tolist())), {
            "held_surfaces": held,
            "extrapolated_surfaces": extrapolated,
            "travel_limited_surfaces": travel_limited,
            "deflection_policy": deflection_policy,
            "assumed_limits_deg": limits if deflection_policy == "linear" else None,
            "beta_held": beta_held,
            "beta_reference_deg": self.beta_deg,
            "combined_surfaces_assumed_additive": len(active) > 1,
            "within_measured_coverage": not held and not beta_held and not extrapolated and not travel_limited,
            "coverage": {name: {
                "alpha_deg": [float(s.alpha_deg[0]), float(s.alpha_deg[-1])],
                "deflection_deg": [float(s.deflection_deg[0]), float(s.deflection_deg[-1])],
            } for name, s in self.surfaces.items()},
        }


def load_deflection_database(path: str | Path, reference: AeroDBReference) -> DeflectionDatabase:
    if Path(path).suffix == ".json":
        from .expanded_deflections import load_expanded_deflections
        return load_expanded_deflections(path, reference)
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream)
        reader.fieldnames = [name.strip() for name in reader.fieldnames or []]
        required = {"alpha", "beta", "velocity", *SURFACES, *LOADS}
        if required - set(reader.fieldnames):
            raise ValueError("Deflection CSV is missing required angles, velocity or body loads")
        rows = []
        for row in reader:
            case = (row.get("case_id") or "unnamed case").strip()
            if None in row or any(row[key] is None for key in required):
                raise ValueError(f"Deflection CSV line {reader.line_num} ({case}) has missing or extra columns")
            values = {}
            for key in sorted(required):
                try:
                    values[key] = float(row[key])
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Deflection CSV line {reader.line_num} ({case}): invalid {key}={row[key]!r}; "
                        "use one decimal point and no thousands separators. Decimal positions must not be guessed."
                    ) from exc
            rows.append(values)
    if not rows or not all(np.all(np.isfinite(list(row.values()))) and row["velocity"] > 0 for row in rows):
        raise ValueError("Deflection rows must be finite with positive velocity")
    betas = {row["beta"] for row in rows}
    if len(betas) != 1:
        raise ValueError("Partial deflection model currently requires a single beta plane")
    baselines = {}
    samples = {name: {} for name in SURFACES}
    seen = set()
    for row in rows:
        key = tuple(row[name] for name in ("alpha", "beta", "velocity", *SURFACES))
        if key in seen:
            raise ValueError("Duplicate deflection case")
        seen.add(key)
        moved = [name for name in SURFACES if row[name] != 0]
        condition = (row["alpha"], row["beta"], row["velocity"])
        if not moved:
            baselines[condition] = row
        elif len(moved) != 1:
            raise ValueError("Independent increment model requires one moved surface per case")
    for row in rows:
        moved = [name for name in SURFACES if row[name] != 0]
        if not moved:
            continue
        condition = (row["alpha"], row["beta"], row["velocity"])
        if condition not in baselines:
            raise ValueError("Deflection case has no matching zero-deflection baseline")
        base = baselines[condition]
        scale = 0.5 * reference.rho * row["velocity"]**2 * reference.S
        denominators = scale * np.array([1, 1, 1, reference.b, reference.c, reference.b])
        delta = np.array([row[key] - base[key] for key in LOADS]) / denominators
        name = moved[0]
        sample_key = (row["alpha"], row[name])
        if sample_key in samples[name]:
            raise ValueError("Multiple velocities at the same deflection point are not supported yet")
        samples[name][sample_key] = delta
        samples[name][(row["alpha"], 0.0)] = np.zeros(6)
    surfaces = {}
    for name, values in samples.items():
        if not values:
            raise ValueError(f"No measured deflection cases for {name}")
        alphas = np.array(sorted({key[0] for key in values}))
        angles = np.array(sorted({key[1] for key in values}))
        table = np.full((len(alphas), len(angles), 6), np.nan)
        for i, alpha in enumerate(alphas):
            for j, angle in enumerate(angles):
                if (alpha, angle) in values:
                    table[i, j] = values[(alpha, angle)]
        if not np.all(np.isfinite(table)):
            raise ValueError(f"Incomplete independent deflection grid for {name}")
        surfaces[name] = SurfaceIncrement(alphas, angles, table)
    return DeflectionDatabase(surfaces, betas.pop(), reference, len(rows))


@lru_cache(maxsize=16)
def _cached(path: str, mtime_ns: int, size: int, reference_key: tuple) -> DeflectionDatabase:
    rho, velocity, area, span, chord, center = reference_key
    return load_deflection_database(path, AeroDBReference(
        rho=rho, V_ref=velocity, S=area, b=span, c=chord,
        moment_center_body_xyz=np.array(center),
    ))


def load_deflection_database_cached(path: str | Path, reference: AeroDBReference) -> DeflectionDatabase:
    path = Path(path).resolve()
    stat = path.stat()
    key = (reference.rho, reference.V_ref, reference.S, reference.b, reference.c,
           tuple(reference.moment_center_body_xyz.tolist()))
    return _cached(str(path), stat.st_mtime_ns, stat.st_size, key)
