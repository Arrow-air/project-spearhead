"""Sparse alpha/deflection planes, speed dependence and measured interactions.

The sidecar explicitly opts into this model. Historical rectangular CSVs keep
their original loader. Neutral sideslip samples are audit data, not a second
base database. All increments use matching alpha/beta/speed zero baselines.
"""

import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull

from .aerodb import BODY_COEFFICIENT_COLUMNS
from .deflections import LOADS, SURFACES


@dataclass(frozen=True)
class SparseSurface:
    points: np.ndarray
    coefficients: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "alpha_deg", np.unique(self.points[:, 0]))
        object.__setattr__(self, "deflection_deg", np.unique(self.points[:, 1]))
        object.__setattr__(self, "_interpolator", LinearNDInterpolator(self.points, self.coefficients))
        object.__setattr__(self, "hull", self.points[ConvexHull(self.points).vertices])

    def bounds_at(self, alpha):
        cuts = []
        for first, second in zip(self.hull, np.roll(self.hull, -1, axis=0)):
            if first[0] == second[0]:
                if abs(alpha-first[0]) < 1e-9:
                    cuts.extend([first[1], second[1]])
            elif min(first[0], second[0])-1e-9 <= alpha <= max(first[0], second[0])+1e-9:
                cuts.append(first[1]+(alpha-first[0])/(second[0]-first[0])*(second[1]-first[1]))
        return min(cuts), max(cuts)

    def contains(self, alpha, angle):
        return bool(self._interpolator.tri.find_simplex([alpha, angle]) >= 0)

    def at(self, alpha_deg, deflection_deg, *, extrapolate=False):
        alpha = float(np.clip(alpha_deg, self.alpha_deg[0], self.alpha_deg[-1]))
        direct = np.asarray(self._interpolator([[alpha, deflection_deg]])[0])
        if np.isfinite(direct).all():
            return direct, alpha != alpha_deg
        lo, hi = self.bounds_at(alpha)
        angle = float(np.clip(deflection_deg, lo, hi))
        value = np.asarray(self._interpolator([[alpha, angle]])[0])
        if not np.isfinite(value).all():
            raise ValueError("Sparse control lookup left the measured convex hull")
        outside = not self.contains(alpha_deg, deflection_deg)
        if extrapolate and angle != deflection_deg:
            # Continue the boundary triangle's deflection slope at fixed alpha.
            inside = angle + min(.001, (hi-lo)/100) * (1 if angle == lo else -1)
            inner = self._interpolator([[alpha, inside]])[0]
            value = value + (deflection_deg-angle)*(value-inner)/(angle-inside)
        return value, outside


@dataclass(frozen=True)
class ExpandedDeflections:
    planes: dict
    interactions: list
    reference: object
    case_count: int
    audit_cases: list
    beta_deg: float = 0.

    @property
    def surfaces(self):
        return {name: planes[25.] for name, planes in self.planes.items()}

    def increments_at(self, alpha_deg, beta_deg, control, *, deflection_policy="hold",
                      limits_deg=(5., 5., 5.), velocity_m_s=None):
        if deflection_policy not in ("hold", "linear"):
            raise ValueError("Deflection policy must be hold or linear")
        if len(limits_deg) != 3 or not np.isfinite(limits_deg).all() or not all(0 < v <= 45 for v in limits_deg):
            raise ValueError("Surface limits must contain three positive angles up to 45 degrees")
        speed = 25. if velocity_m_s is None else float(velocity_m_s)
        if not np.isfinite([alpha_deg, beta_deg, speed]).all() or speed < 0:
            raise ValueError("Control lookup coordinates must be finite with nonnegative speed")
        for measured_speed in (10., 15., 25.):
            if abs(speed-measured_speed) < 1e-9:
                speed = measured_speed
        total = np.zeros(6)
        held, extrapolated, limited, velocity_held, active = [], [], [], [], []
        limits = dict(zip(("elevator", "aileron", "rudder"), limits_deg))
        used_angles = {}
        for name, field in SURFACES.items():
            angle = float(np.rad2deg(getattr(control, field)))
            if not np.isfinite(angle):
                raise ValueError("Control deflections must be finite")
            used = float(np.clip(angle, -limits[name], limits[name])) if deflection_policy == "linear" else angle
            used_angles[name] = used
            if angle == 0:
                continue
            active.append(name)
            if used != angle:
                limited.append(name)
            planes = self.planes[name]
            speeds = sorted(planes)
            clipped = float(np.clip(speed, speeds[0], speeds[-1]))
            upper = int(np.searchsorted(speeds, clipped))
            upper = min(upper, len(speeds)-1)
            lower = upper if speeds[upper] == clipped else max(0, upper-1)
            selected = [speeds[lower]] if lower == upper else [speeds[lower], speeds[upper]]
            # Evaluate every participating plane with the same explicit hold /
            # bounded extrapolation policy before interpolating coefficients.
            # Switching abruptly to 25 m/s outside a sparse hull would create
            # an artificial force jump at alpha=0 or deflection=5 degrees.
            if clipped != speed:
                velocity_held.append(name)
            values = []
            for velocity in selected:
                surface = planes[velocity]
                value, outside = surface.at(alpha_deg, used, extrapolate=deflection_policy == "linear")
                values.append(value)
                if outside and name not in held:
                    held.append(name)
                if outside:
                    a = float(np.clip(alpha_deg, surface.alpha_deg[0], surface.alpha_deg[-1]))
                    lo, hi = surface.bounds_at(a)
                    if deflection_policy == "linear" and not lo-1e-9 <= used <= hi+1e-9 and name not in extrapolated:
                        extrapolated.append(name)
            total += values[0] if len(values) == 1 else values[0]+(clipped-selected[0])/(selected[1]-selected[0])*(values[1]-values[0])
        interaction_applied = False
        interaction_exact = False
        de, dr = used_angles["elevator"], used_angles["rudder"]
        if de < 0 and dr != 0:
            sign = np.sign(dr)
            cases = sorted([r for r in self.interactions if np.sign(r["rudder"]) == sign], key=lambda r: r["alpha"])
            alpha = np.array([r["alpha"] for r in cases])
            residual = np.array([r["residual"] for r in cases])
            value = np.array([np.interp(alpha_deg, alpha, residual[:, j]) for j in range(6)])
            measured_de = np.interp(alpha_deg, alpha, [r["elevator"] for r in cases])
            total += value*min(1., de/measured_de)*min(1., abs(dr)/15.)
            interaction_applied = True
            interaction_exact = any(abs(alpha_deg-r["alpha"]) < 1e-9 and abs(de-r["elevator"]) < 1e-9
                                    and abs(dr-r["rudder"]) < 1e-9 and speed == 25. and beta_deg == 0.
                                    for r in cases)
        beta_held = bool(active) and abs(beta_deg) > 1e-9
        coverage = {name: {"alpha_deg": [float(s.alpha_deg[0]), float(s.alpha_deg[-1])],
                           "deflection_deg": [float(s.deflection_deg[0]), float(s.deflection_deg[-1])],
                           "convex_hull_alpha_deflection_deg": s.hull.tolist(),
                           "measured_velocities_m_s": sorted(self.planes[name])}
                    for name, s in self.surfaces.items()}
        return dict(zip(BODY_COEFFICIENT_COLUMNS, total.tolist())), {
            "held_surfaces": held, "extrapolated_surfaces": extrapolated,
            "travel_limited_surfaces": limited, "velocity_held_surfaces": velocity_held,
            "deflection_policy": deflection_policy, "assumed_limits_deg": limits if deflection_policy == "linear" else None,
            "beta_held": beta_held, "beta_reference_deg": 0., "velocity_m_s": speed,
            "combined_surfaces_assumed_additive": len(active) > 1 and not interaction_exact,
            "interaction_residual_applied": interaction_applied, "interaction_at_measured_case": interaction_exact,
            "interaction_model": "alpha interpolation; bounded product of negative elevator and signed rudder fractions",
            "within_measured_coverage": not held and not beta_held and not extrapolated and not limited
                                        and not velocity_held and (not interaction_applied or interaction_exact),
            "coverage": coverage}


def load_expanded_deflections(path, reference):
    path = Path(path)
    metadata = json.loads(path.read_text())
    if metadata["model"] != "sparse_controls_v1":
        raise ValueError("Unknown expanded deflection model")
    source = path.parent/metadata["source"]
    if hashlib.sha256(source.read_bytes()).hexdigest() != metadata["source_sha256"]:
        raise ValueError("Expanded control source hash mismatch")
    with source.open(newline="") as stream:
        raw = list(csv.DictReader(stream))
    fields = ("alpha", "beta", "velocity", *SURFACES, *LOADS)
    rows, seen = [], set()
    for row in raw:
        if None in row or any(row.get(k) is None for k in fields):
            raise ValueError(f"Malformed expanded control case {row.get('case_id')}")
        r = {k: float(row[k]) for k in fields}
        if not np.isfinite(list(r.values())).all() or r["velocity"] <= 0:
            raise ValueError("Expanded control cases must be finite with positive speed")
        r["case_id"] = row["case_id"]
        key = tuple(r[k] for k in fields[:6])
        if key in seen:
            raise ValueError("Duplicate expanded control case")
        seen.add(key)
        rows.append(r)
    baselines = {(r["alpha"], r["beta"], r["velocity"]): r for r in rows
                 if all(r[k] == 0 for k in SURFACES)}
    samples = {name: {} for name in SURFACES}
    combinations, audit = [], []
    for r in rows:
        moved = [k for k in SURFACES if r[k] != 0]
        if r["beta"] != 0:
            if moved:
                raise ValueError("Nonzero-beta controls require a measured beta-dependent model")
            audit.append(r["case_id"])
            continue
        if not moved:
            continue
        condition = (r["alpha"], r["beta"], r["velocity"])
        if condition not in baselines:
            raise ValueError("Expanded case has no matching zero-deflection baseline")
        base = baselines[condition]
        scale = .5*reference.rho*r["velocity"]**2*reference.S*np.array([1, 1, 1, reference.b, reference.c, reference.b])
        delta = np.array([r[k]-base[k] for k in LOADS])/scale
        if len(moved) == 1:
            plane = samples[moved[0]].setdefault(r["velocity"], {})
            plane[(r["alpha"], r[moved[0]])] = delta
            plane[(r["alpha"], 0.)] = np.zeros(6)
        elif set(moved) == {"elevator", "rudder"} and r["velocity"] == 25. and r["elevator"] < 0 and abs(r["rudder"]) == 15:
            combinations.append((r, delta))
        else:
            raise ValueError("Unsupported measured surface interaction")
    planes = {name: {speed: SparseSurface(np.array(list(points)), np.array(list(points.values())))
                     for speed, points in speeds.items()} for name, speeds in samples.items()}
    if any(25. not in speeds for speeds in planes.values()):
        raise ValueError("Each surface requires a 25 m/s reference plane")
    interactions = []
    for r, measured in combinations:
        independent = sum((planes[name][25.].at(r["alpha"], r[name])[0] for name in ("elevator", "rudder")), np.zeros(6))
        interactions.append({**r, "residual": measured-independent})
    if len(interactions) != 4 or {r["alpha"] for r in interactions} != {0., 18.}:
        raise ValueError("Interaction model requires both rudder signs at alpha 0 and 18")
    return ExpandedDeflections(planes, interactions, reference, len(rows), audit)
