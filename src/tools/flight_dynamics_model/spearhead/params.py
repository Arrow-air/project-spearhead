"""Aircraft, geometry, and aerodynamic reference parameters."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


def default_aero_database_path() -> Path:
    project_root = Path(__file__).resolve().parents[4]
    return (
        project_root
        / "docs/information-note/phase-1/Aerodynamics/"
        "0002-Preliminary-Design-Aerodynamic-Database/assets/"
        "adb_v1_1_drag_correction.csv"
    )


@dataclass(frozen=True)
class AeroDerivatives:
    CL0: float = 0.25
    CLalpha: float = 4.8
    CLde: float = 0.35
    CD0: float = 0.04
    k: float = 0.08
    CYbeta: float = -0.6
    CYdr: float = 0.18
    Clbeta: float = -0.08
    Clda: float = 0.08
    Clp: float = -0.45
    Clr: float = 0.08
    Cm0: float = 0.02
    Cmalpha: float = -1.2
    Cmde: float = -0.9
    Cmq: float = -8.0
    Cnbeta: float = 0.12
    Cndr: float = -0.08
    Cnp: float = -0.02
    Cnr: float = -0.25


@dataclass(frozen=True)
class AircraftGeometry:
    layout: str = "high-wing quadplane pusher"
    fuselage_length: float = 1.095
    fuselage_max_width: float = 0.240
    fuselage_max_height: float = 0.275
    wing_airfoil: str = "Clark Z"
    wing_gross_span: float = 3.95
    wing_effective_span: float = 3.65
    wing_chord: float = 0.457
    wing_gross_area: float = 1.805
    wing_effective_area: float = 1.668
    wing_incidence: float = np.deg2rad(1.0)
    wing_ac_x_cg: float = -0.15
    tail_type: str = "inverse V-tail"
    tail_airfoil_design_note: str = "NACA 0015"
    tail_airfoil_config_py: str = "NACA 0018"
    tail_airfoil_validation_note: str = (
        "TODO: PR #7 design note says NACA 0015, but assets/config.py says NACA 0018."
    )
    tail_projected_span: float = 1.10
    tail_geometric_span: float = 1.343
    tail_chord: float = 0.28
    tail_area: float = 0.376
    tail_dihedral: float = np.deg2rad(35.0)
    tail_incidence: float = np.deg2rad(2.0)
    tail_ac_x_cg: float = 1.15
    boom_lateral_station: float = 0.55
    vtol_lift_station_x_abs: float = 0.70


@dataclass(frozen=True)
class AeroDBReference:
    rho: float = 1.225
    V_ref: float = 25.0
    q: float = 382.8125
    S: float = 1.805
    b: float = 3.95
    c: float = 0.457
    qS: float = 690.9766
    qSb: float = 2729.3574
    qSc: float = 315.7763
    moment_center_body_xyz: np.ndarray = field(default_factory=lambda: np.array([-0.15, 0.0, 0.0]))
    alpha_min_deg: float = -90.0
    alpha_max_deg: float = 90.0
    beta_min_deg: float = -90.0
    beta_max_deg: float = 90.0


@dataclass(frozen=True)
class AircraftParams:
    g: float = 9.80665
    rho: float = 1.225
    mass: float = 25.0
    inertia: np.ndarray = field(default_factory=lambda: np.diag([2.5, 4.0, 5.0]))
    Sref: float = 1.805
    bref: float = 3.95
    cref: float = 0.457
    Vmin: float = 0.1
    alpha_limit: float = np.deg2rad(90.0)
    beta_limit: float = np.deg2rad(90.0)
    geometry: AircraftGeometry = field(default_factory=AircraftGeometry)
    aerodb_reference: AeroDBReference = field(default_factory=AeroDBReference)
    aero_database_path: Path = field(default_factory=default_aero_database_path)
    use_aero_database: bool = True
    aero_database_clamp: bool = True
    shift_adb_moments_to_cg: bool = True
    aero: AeroDerivatives = field(default_factory=AeroDerivatives)
    Tmax: float = 120.0
