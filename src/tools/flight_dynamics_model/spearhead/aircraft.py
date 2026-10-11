"""Load aircraft hardware profiles without duplicating the simulation engine."""

from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
import hashlib
from pathlib import Path
from typing import Any

import numpy as np
import yaml

FDM_ROOT = Path(__file__).resolve().parents[1]
AIRCRAFT_CONFIGS = FDM_ROOT / "aircraft_configs"
AERODYNAMIC_CONFIGS = FDM_ROOT / "aerodynamic_configs"


def aerodynamic_options() -> dict[str, str]:
    return {data["id"]: data["label"] for path in sorted(AERODYNAMIC_CONFIGS.glob("*.yaml"))
            for data in [yaml.safe_load(path.read_text())]}


def select_aerodynamics(params, dataset_id: str):
    """Select a versioned aerodynamic dataset independently of the hardware."""
    from .aerodb import load_aerodynamic_database_cached
    from .deflections import load_deflection_database_cached

    if dataset_id not in aerodynamic_options():
        raise ValueError(f"Unknown aerodynamic dataset: {dataset_id}")
    path = AERODYNAMIC_CONFIGS / f"{dataset_id}.yaml"
    data = yaml.safe_load(path.read_text())
    database_path = (path.parent / data["path"]).resolve()
    database = load_aerodynamic_database_cached(database_path)
    deflection_path = data.get("deflection_path")
    deflection_path = None if deflection_path is None else (path.parent / deflection_path).resolve()
    if deflection_path is not None:
        load_deflection_database_cached(deflection_path, database.reference)
    return replace(params, aero_dataset_id=dataset_id,
                   aero_database_version=data["version"], aero_database_path=database_path,
                   deflection_database_path=deflection_path, aerodb_reference=database.reference,
                   rho=database.reference.rho, Sref=database.reference.S,
                   bref=database.reference.b, cref=database.reference.c,
                   geometry_version=data["geometry_version"],
                   geometry=replace(params.geometry, **data.get("geometry", {})))


@dataclass(frozen=True)
class BatterySpec:
    cells: int = 12
    voltage_v: float = 44.4
    capacity_ah: float = 22.0
    mass_kg: float = 3.709

    def __post_init__(self) -> None:
        if self.cells <= 0 or not np.isfinite(self.voltage_v) or self.voltage_v <= 0:
            raise ValueError("Battery cells and voltage must be positive")


@dataclass(frozen=True)
class LiftSystem:
    motor: dict[str, Any] = field(default_factory=dict)
    esc: dict[str, Any] = field(default_factory=dict)
    propeller: dict[str, Any] = field(default_factory=dict)
    positions_body_m: tuple[tuple[float, float, float], ...] = ()
    reaction_torque_signs: tuple[int, ...] = ()
    performance_path: Path | None = None

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions_body_m, dtype=float)
        if positions.ndim != 2 or positions.shape[1] != 3 or not np.all(np.isfinite(positions)):
            raise ValueError("Lift rotor positions must have shape (N, 3) and be finite")
        if len(self.reaction_torque_signs) != len(positions) or any(
            sign not in (-1, 1) for sign in self.reaction_torque_signs
        ):
            raise ValueError("Each lift rotor needs a reaction torque sign of -1 or +1")


@dataclass(frozen=True)
class PusherSystem:
    motor: dict[str, Any]
    esc: dict[str, Any]
    propeller: dict[str, Any]
    performance_path: Path
    equivalent_resistance_ohm: float
    esc_efficiency: float = 0.97
    voltage_min_v: float = 36.0
    voltage_max_v: float = 50.4
    max_current_a: float = 77.8
    max_input_power_w: float = 3569.0
    sustainable_current_a: float = 24.0
    # No measured thrust line or rotation direction supplied: use CG-aligned
    # thrust and omit reaction torque until installation data is available.
    position_body_m: tuple[float, float, float] | None = None
    reaction_torque_sign: int = 0

    def __post_init__(self):
        positive = [self.equivalent_resistance_ohm, self.voltage_min_v, self.voltage_max_v,
                    self.max_current_a, self.max_input_power_w, self.sustainable_current_a,
                    self.motor["kv_rpm_per_v"], self.propeller["diameter_in"], self.propeller["max_rpm"]]
        if not np.all(np.isfinite(positive)) or min(positive) <= 0:
            raise ValueError("Pusher constants must be finite and positive")
        if (not 0 < self.esc_efficiency <= 1 or self.voltage_min_v >= self.voltage_max_v
                or not 0 <= self.motor["no_load_current_a"] < self.max_current_a
                or not self.sustainable_current_a <= self.max_current_a <= self.esc["continuous_current_a"]):
            raise ValueError("Invalid pusher voltage, efficiency or current limits")
        if self.position_body_m is not None and (np.shape(self.position_body_m) != (3,)
                                                 or not np.all(np.isfinite(self.position_body_m))):
            raise ValueError("Pusher position must contain three finite values")
        if self.reaction_torque_sign not in (-1, 0, 1):
            raise ValueError("Pusher reaction torque sign must be -1, 0 or +1")


def aircraft_options() -> dict[str, str]:
    """Return the bundled profiles for UI/CLI selection."""
    options = {}
    for path in sorted(AIRCRAFT_CONFIGS.glob("*.yaml")):
        data = yaml.safe_load(path.read_text())
        options[data["id"]] = data["label"]
    return options


def load_aircraft_config(name_or_path: str | Path = "reference_design", *, aerodynamics: str | None = None):
    """Build independent parameters from a bundled profile or a YAML file.

    Relative AeroDB and propulsion paths are resolved against the profile's
    directory. Geometry changes and new databases remain explicit profile inputs.
    """
    from .params import AircraftGeometry, AircraftParams

    path = Path(name_or_path)
    if not path.is_file():
        path = AIRCRAFT_CONFIGS / f"{name_or_path}.yaml"
    data = yaml.safe_load(path.read_text())
    allowed = {"id", "label", "description", "aircraft", "geometry", "aerodynamics", "battery", "propulsion"}
    if not isinstance(data, dict) or set(data) - allowed:
        raise ValueError(f"Invalid aircraft profile: {path}")
    values = dict(data.get("aircraft", {}))
    allowed_params = {"mass", "inertia", "cg_body_xyz", "rho"}
    if set(values) - allowed_params:
        raise ValueError("Unsupported aircraft parameter in profile")
    for key in ("inertia", "cg_body_xyz"):
        if key in values:
            values[key] = np.asarray(values[key], dtype=float)
    if "mass" in values and (not np.isfinite(values["mass"]) or values["mass"] <= 0):
        raise ValueError("Aircraft mass must be positive")
    if "cg_body_xyz" in values and (values["cg_body_xyz"].shape != (3,) or not np.all(np.isfinite(values["cg_body_xyz"]))):
        raise ValueError("Aircraft CG must contain three finite values")
    if "inertia" in values:
        inertia = values["inertia"]
        if inertia.shape != (3, 3) or not np.all(np.isfinite(inertia)) or not np.allclose(inertia, inertia.T) or np.min(np.linalg.eigvalsh(inertia)) <= 0:
            raise ValueError("Aircraft inertia must be symmetric positive definite with shape (3, 3)")
    geometry = data.get("geometry", {})
    if set(geometry) - {item.name for item in fields(AircraftGeometry)}:
        raise ValueError("Unsupported aircraft geometry field")
    aero = data.get("aerodynamics", {})
    propulsion = data.get("propulsion", {})
    pusher = propulsion.get("pusher", {})
    pusher_system = dict(pusher.get("estimate", {}))
    if pusher_system:
        pusher_system["performance_path"] = (path.parent / pusher_system["performance_path"]).resolve()
        if pusher_system.get("position_body_m") is not None:
            pusher_system["position_body_m"] = tuple(pusher_system["position_body_m"])
    lift = dict(propulsion.get("lift", {}))
    curve = lift.pop("performance_path", None)
    if curve is not None:
        lift["performance_path"] = (path.parent / curve).resolve()
    for key in ("positions_body_m", "reaction_torque_signs"):
        if key in lift:
            lift[key] = tuple(tuple(row) if isinstance(row, list) else row for row in lift[key])
    values.update(
        aircraft_id=str(data["id"]),
        aircraft_label=str(data["label"]),
        aircraft_description=str(data.get("description", "")),
        geometry=AircraftGeometry(**geometry),
        geometry_version=str(aero.get("geometry_version", "preliminary_tail_chord_0.28m")),
        aero_database_version=str(aero.get("version", "adb_v1.1")),
        battery=BatterySpec(**data.get("battery", {})),
        pusher_installed=bool(pusher.get("installed", True)),
        pusher_system=PusherSystem(**pusher_system) if pusher_system else None,
        Tmax=float(pusher.get("max_thrust_n", 120.0)),
        lift_system=LiftSystem(**lift) if lift else None,
    )
    if "path" in aero:
        values["aero_database_path"] = (path.parent / aero["path"]).resolve()
        if not values["aero_database_path"].is_file():
            raise ValueError(f"Profile AeroDB file is missing: {values['aero_database_path']}")
    params = AircraftParams(**values)
    dataset_id = aerodynamics or aero.get("dataset")
    params = params if dataset_id is None else select_aerodynamics(params, str(dataset_id))
    if params.pusher_system is not None:
        from .pusher import estimate_pusher
        params = replace(params, Tmax=estimate_pusher(params.pusher_system, params.battery.voltage_v,
                                                     1, 0, params.rho)["thrust_n"])
    return params


def aircraft_summary(params) -> dict[str, Any]:
    """Record enough profile identity to distinguish saved simulation runs."""
    lift = params.lift_system
    control_data = None
    if params.deflection_database_path is not None:
        from .deflections import load_deflection_database_cached
        database = load_deflection_database_cached(params.deflection_database_path, params.aerodb_reference)
        from .controls import Control
        _, control_data = database.increments_at(0, 0, Control(0, 0, 0, 0))
        control_data = {"coverage": control_data["coverage"], "beta_deg": database.beta_deg,
                        "outside_coverage_policy": "hold endpoint and report",
                        "combined_surfaces_policy": "add independent increments; interactions unmeasured",
                        "rate_derivatives": "provisional"}
    return {
        "id": params.aircraft_id,
        "label": params.aircraft_label,
        "description": params.aircraft_description,
        "mass_kg": params.mass,
        "cg_body_xyz_m": params.cg_body_xyz.tolist(),
        "inertia_kg_m2": params.inertia.tolist(),
        "tail_chord_m": params.geometry.tail_chord,
        "geometry_version": params.geometry_version,
        "aero_database_version": params.aero_database_version,
        "aero_database_path": str(params.aero_database_path),
        "aero_dataset_id": params.aero_dataset_id,
        "rho_kg_m3": params.rho,
        "moment_center_body_xyz_m": params.aerodb_reference.moment_center_body_xyz.tolist(),
        "deflection_database_path": None if params.deflection_database_path is None else str(params.deflection_database_path),
        "control_data": control_data,
        "battery_voltage_v": params.battery.voltage_v,
        "pusher_installed": params.pusher_installed,
        "pusher_max_thrust_n": params.Tmax,
        "pusher_estimate": None if params.pusher_system is None else {
            "motor": params.pusher_system.motor, "esc": params.pusher_system.esc,
            "propeller": params.pusher_system.propeller,
            "performance_path": str(params.pusher_system.performance_path),
            "performance_sha256": hashlib.sha256(params.pusher_system.performance_path.read_bytes()).hexdigest(),
            "equivalent_resistance_ohm": params.pusher_system.equivalent_resistance_ohm,
            "esc_efficiency": params.pusher_system.esc_efficiency,
            "sustainable_current_a": params.pusher_system.sustainable_current_a,
            "max_current_a": params.pusher_system.max_current_a,
            "max_input_power_w": params.pusher_system.max_input_power_w,
            "voltage_range_v": [params.pusher_system.voltage_min_v, params.pusher_system.voltage_max_v],
            "position_body_m": params.pusher_system.position_body_m,
            "reaction_torque_sign": params.pusher_system.reaction_torque_sign,
            "model": "steady motor equilibrium + APC analytical Ct/Cp; provisional",
        },
        "lift_motor": None if lift is None else lift.motor,
        "lift_esc": None if lift is None else lift.esc,
        "lift_propeller": None if lift is None else lift.propeller,
        "lift_performance_path": None if lift is None or lift.performance_path is None else str(lift.performance_path),
    }
