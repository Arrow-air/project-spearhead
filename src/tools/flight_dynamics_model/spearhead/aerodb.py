"""Aerodynamic database loading, interpolation, and dimensionalization."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .params import AeroDBReference

BODY_COEFFICIENT_COLUMNS = ("cfx", "cfy", "cfz", "cmx", "cmy", "cmz")
WIND_COEFFICIENT_COLUMNS = ("cd", "cs", "cl", "cr", "cm", "cn")
COEFFICIENT_COLUMNS = BODY_COEFFICIENT_COLUMNS + WIND_COEFFICIENT_COLUMNS
_CACHE: dict[tuple, "AeroDatabase"] = {}


def _normalize_column_name(name: str) -> str:
    return name.strip().lower().lstrip("#").strip()


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def _find_header_row(rows: list[list[str]]) -> int:
    for i, row in enumerate(rows):
        normalized = [_normalize_column_name(value) for value in row]
        has_alpha = any(name in normalized for name in ("alpha_deg", "alpha", "alpha(deg)"))
        has_beta = any(name in normalized for name in ("beta_deg", "beta", "beta(deg)"))
        has_body_coefficients = all(column in normalized for column in BODY_COEFFICIENT_COLUMNS)
        if has_alpha and has_beta and has_body_coefficients:
            return i
    raise ValueError("Could not find ADB header row with alpha, beta, and body coefficients")


def _find_column(column_index: dict[str, int], candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        if candidate in column_index:
            return candidate
    raise ValueError(f"Could not find any of these ADB columns: {candidates}")


def _metadata_from_rows(rows: list[list[str]]) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for row in rows:
        clean = [value.strip() for value in row if value.strip()]
        if len(clean) < 2:
            continue
        key = clean[0].lstrip("#").strip()
        if not key or _is_number(key):
            continue
        metadata[key] = ",".join(clean[1:])
    return metadata


@dataclass(frozen=True)
class AeroDatabase:
    alpha_grid_deg: np.ndarray
    beta_grid_deg: np.ndarray
    coefficient_tables: dict[str, np.ndarray]
    metadata: dict[str, str]
    reference: AeroDBReference

    def __post_init__(self):
        interpolators = {
            column: RegularGridInterpolator(
                (self.alpha_grid_deg, self.beta_grid_deg),
                table,
                bounds_error=True,
            )
            for column, table in self.coefficient_tables.items()
        }
        object.__setattr__(self, "_interpolators", interpolators)

    @property
    def alpha_bounds_deg(self) -> tuple[float, float]:
        return float(self.alpha_grid_deg[0]), float(self.alpha_grid_deg[-1])

    @property
    def beta_bounds_deg(self) -> tuple[float, float]:
        return float(self.beta_grid_deg[0]), float(self.beta_grid_deg[-1])

    def coefficients_at(
        self,
        alpha_deg: float,
        beta_deg: float,
        *,
        clamp: bool = False,
    ) -> dict[str, float]:
        alpha_query = float(alpha_deg)
        beta_query = float(beta_deg)
        if clamp:
            alpha_query = np.clip(alpha_query, *self.alpha_bounds_deg)
            beta_query = np.clip(beta_query, *self.beta_bounds_deg)
        point = np.array([[alpha_query, beta_query]])
        try:
            return {
                column: float(interpolator(point)[0])
                for column, interpolator in self._interpolators.items()
            }
        except ValueError as exc:
            raise ValueError(
                "ADB interpolation requested outside alpha/beta bounds: "
                f"alpha={alpha_deg} deg in {self.alpha_bounds_deg}, "
                f"beta={beta_deg} deg in {self.beta_bounds_deg}"
            ) from exc


def load_aerodynamic_database(
    path: str | Path,
    reference: AeroDBReference | None = None,
) -> AeroDatabase:
    """Load a Nondimit-style coefficient CSV with optional metadata rows."""
    db_path = Path(path)
    if reference is None:
        reference = AeroDBReference()

    with db_path.open(newline="") as csv_file:
        rows = list(csv.reader(csv_file))

    header_index = _find_header_row(rows)
    metadata = _metadata_from_rows(rows[:header_index])
    header = [_normalize_column_name(value) for value in rows[header_index]]
    column_index = {name: i for i, name in enumerate(header)}
    alpha_name = _find_column(column_index, ("alpha_deg", "alpha", "alpha(deg)"))
    beta_name = _find_column(column_index, ("beta_deg", "beta", "beta(deg)"))
    available_coefficients = [column for column in COEFFICIENT_COLUMNS if column in column_index]

    numeric_rows = []
    for row in rows[header_index + 1 :]:
        if not row or not any(value.strip() for value in row):
            continue
        first = row[0].strip()
        if first.startswith("#"):
            continue
        if not _is_number(row[column_index[alpha_name]].strip()):
            continue
        numeric_rows.append(row)

    alpha_values = np.array([float(row[column_index[alpha_name]]) for row in numeric_rows])
    beta_values = np.array([float(row[column_index[beta_name]]) for row in numeric_rows])
    alpha_grid = np.unique(alpha_values)
    beta_grid = np.unique(beta_values)

    tables = {
        column: np.full((alpha_grid.size, beta_grid.size), np.nan)
        for column in available_coefficients
    }
    alpha_lookup = {value: i for i, value in enumerate(alpha_grid)}
    beta_lookup = {value: i for i, value in enumerate(beta_grid)}

    for row in numeric_rows:
        alpha = float(row[column_index[alpha_name]])
        beta = float(row[column_index[beta_name]])
        alpha_i = alpha_lookup[alpha]
        beta_i = beta_lookup[beta]
        for column in available_coefficients:
            tables[column][alpha_i, beta_i] = float(row[column_index[column]])

    missing = [
        column
        for column, table in tables.items()
        if np.any(~np.isfinite(table))
    ]
    if missing:
        raise ValueError(f"ADB grid is incomplete for columns: {missing}")

    return AeroDatabase(alpha_grid, beta_grid, tables, metadata, reference)


def _reference_cache_key(reference: AeroDBReference) -> tuple:
    return (
        reference.rho,
        reference.V_ref,
        reference.q,
        reference.S,
        reference.b,
        reference.c,
        tuple(reference.moment_center_body_xyz.tolist()),
    )


def load_aerodynamic_database_cached(path: str | Path, reference: AeroDBReference) -> AeroDatabase:
    db_path = Path(path).resolve()
    key = (str(db_path), _reference_cache_key(reference))
    if key not in _CACHE:
        _CACHE[key] = load_aerodynamic_database(db_path, reference)
    return _CACHE[key]


def dimensional_body_force_moment(
    coefficients: dict[str, float],
    qbar: float,
    reference: AeroDBReference,
    *,
    shift_to_cg: bool = True,
    cg_body_xyz: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert body-axis coefficients into dimensional loads.

    If requested, moments are shifted from the ADB moment center to the CG with
    M_CG = M_db + cross(old_center - CG, F_body).
    """
    force_b = qbar * reference.S * np.array(
        [coefficients["cfx"], coefficients["cfy"], coefficients["cfz"]]
    )
    moment_b = qbar * reference.S * np.array(
        [
            reference.b * coefficients["cmx"],
            reference.c * coefficients["cmy"],
            reference.b * coefficients["cmz"],
        ]
    )

    if shift_to_cg:
        if cg_body_xyz is None:
            cg_body_xyz = np.zeros(3)
        r_new_to_old = reference.moment_center_body_xyz - np.asarray(cg_body_xyz, dtype=float)
        moment_b = moment_b + np.cross(r_new_to_old, force_b)

    return force_b, moment_b
