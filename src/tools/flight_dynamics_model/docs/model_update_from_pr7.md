# Model Update From PR #7

This note records the simulation-model changes needed after Arrow-air/project-spearhead PR #7 merged the preliminary geometry and ADB v1.1 documentation.

## Changed Assumptions

- The model reference geometry has been updated from the initial placeholder values to the preliminary high-wing quadplane pusher layout.
- Mass remains `25 kg`.
- Main reference values are now `S = 1.805 m^2`, `b = 3.95 m`, and `c = 0.457 m`.
- The preferred fixed-wing aerodynamic source is ADB v1.1 drag-corrected CSV:
  `docs/information-note/phase-1/Aerodynamics/0002-Preliminary-Design-Aerodynamic-Database/assets/adb_v1_1_drag_correction.csv`.
- Body-axis coefficients are used for equations of motion:
  `Fx = cfx*qS`, `Fy = cfy*qS`, `Fz = cfz*qS`,
  `Mx = cmx*qSb`, `My = cmy*qSc`, `Mz = cmz*qSb`.
- ADB moments are referenced to `[-0.15, 0.0, 0.0] m`; the simulation shifts them to CG using `M_CG = M_db + cross(r, F_body)`.

## Remaining Uncertainties

- The latest `project-spearhead` repository includes the large ADB CSV and the flight dynamics tool loads it by default. The fallback derivative model is preserved for cases where the CSV is unavailable.
- The inertia matrix is still a placeholder and has not been updated from a mass-properties model.
- The wing/tail/fuselage force split is still bookkeeping only, not validated component physics.

## Tail Airfoil Inconsistency

- The preliminary design note says the tail airfoil is `NACA 0015`.
- `assets/config.py` in PR #7 says `tail_airfoil = "NACA 0018"`.
- This is visible in `AircraftGeometry` and needs validation before freezing geometry assumptions.

## Control-Surface Database Missing

ADB v1.1 only covers `alpha` and `beta`; it does not include elevator, ruddervator, aileron, flap, or other control-surface dimensions. The model keeps provisional derivative increments separate from the ADB interpolation so they can be replaced cleanly when control-surface tables exist.

## Propulsion / ADM Limitation

The ADB workflow includes an actuator-disk correction based on pusher thrust cases at `30 N`. This is not a throttle-dependent propulsion model. The simulation keeps pusher thrust and pusher moment in the propulsion model rather than treating the ADM correction as throttle-varying propulsion.
