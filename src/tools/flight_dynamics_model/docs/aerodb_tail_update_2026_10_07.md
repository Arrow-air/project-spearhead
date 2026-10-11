# October 2026 AeroDB comparison and integration

Received `database.csv` and `deflection.csv` on 2026-10-07. Both use user-confirmed body-axis moment center (0,0,0) m and density 1.021 kg/m³. Original source bytes are preserved in `sources/`, with SHA-256 hashes in `source.json`. The runtime database changes only the contradictory moment-center metadata; coefficients and dimensional loads are unchanged.

The later `results.csv` control expansion is integrated as of 9 October. See the [update assessment and multi-speed validation](../../../../artifacts/control_results_2026_10_09/assessment.md) for current results, provenance and limitations.

## Quick comparison

Both base databases contain the complete 181 × 181 alpha/beta grid, −90..90° at 1° steps. The table compares coefficients directly. Moments are evaluated at the same origin, density 1.021 kg/m³ and speed 25 m/s; the old moments are transferred from their documented (−0.15,0,0.15) m center. Density alone does not explain coefficient differences. Updated tail geometry is expected to change loads, but these files do not isolate its effect from export/correction differences.

| α (β=0) | Old CL | New CL | Old CD | New CD | Old My at origin (N·m) | New My at origin (N·m) |
|---:|---:|---:|---:|---:|---:|---:|
| -5° | -0.0531 | -0.0980 | 0.0340 | 0.0427 | -16.44 | -15.48 |
| 0° | 0.3759 | 0.3587 | 0.0343 | 0.0426 | -35.70 | -13.21 |
| 5° | 0.7976 | 0.8173 | 0.0612 | 0.0729 | -46.76 | -5.00 |

At α=0°, CL decreases 4.6% and CD increases 24.2%. Moment comparisons without a common center would mix geometry effects with lever-arm effects.

## Deflection coverage and model

The active 9 October control update contains 49 cases at β=0° and V=25 m/s. All three surfaces have α=−5,0,5,10,15,20,30° coverage, with ±5° and matching zero-deflection baselines. The original 17-case delivery remains preserved in `sources/deflection.csv`; the complete 54-case `results.csv` export is preserved in `sources/results_2026_10_09.csv`. Five additional zero-control sideslip cases are retained for audit, not injected into the full base table. Metadata columns are omitted from the normalized runtime control CSV. All 17 overlapping cases match exactly; source CFD convergence is marked `not_assessed`. Units: angles in degrees, forces in N, moments in N·m; runtime controls remain radians.

Measured increments are computed by subtracting the matching zero-deflection row, then normalized by ρV²S/2 and the wing reference span/chord. They are interpolated linearly in alpha and surface angle and added to the full alpha/beta base database before one moment transfer to the configured CG. This preserves the base database at zero control and reproduces all supplied load deltas. Measured controls replace the provisional control terms; provisional rate damping remains.

Outside measured alpha bounds, hold the nearest endpoint and flag the approximation. The default deflection policy holds measured endpoints; the configured linear policy extrapolates each measured side up to the assumed mechanical travel limit and flags it. At nonzero beta, use the measured beta=0 plane and flag it. Multiple controls add independent increments; interactions have not been measured. No additional 1.3 drag correction is applied to these measured increments. GUI labels and exported aircraft metadata expose coverage and model assumptions. This is a partial local surrogate, not a complete deflection database.

## Export discrepancy

The base CSV header says moment center (−0.1,0,0.15), contrary to the user-confirmed (0,0,0). At α=0°, the base table gives My≈−13.21 N·m, while the raw zero-deflection row gives −36.57 N·m; similar discrepancies occur at ±5°. The header center and force lever arm could account for much of this difference, but no unrequested moment transformation is applied. The user-confirmed origin is authoritative. Baseline subtraction prevents the deflection file from introducing a second base moment. Reconcile this export discrepancy before treating the new database as a validated stability baseline.

## Trim and geometry status

The original 17-case intake study at V=25 m/s, mass=25 kg and CG=(−0.15,0,0.15) m gave a bounded preliminary trim at alpha=-0.557°, elevator=5.000°. It did **not** satisfy force/moment balance within measured elevator bounds ±5°. Residual [Fx/m, Fz/m, My/Iyy, theta−alpha] = [0.0, 1.826671, 1.789204, 0.358784]. This remains historical intake evidence. The later update assessment documents successful diagnostic-reference trims at 18–28 m/s using a different CG and propulsion setup, including elevator extrapolation at 18 and 20 m/s. A solver termination flag alone is insufficient for trim success; initialization checks actual equilibrium residuals.

The numeric updated tail chord/area were not supplied; they are marked pending instead of assigning the old 0.28 m value to the new geometry. Wing reference chord remains 0.457 m, area 1.805 m² and span 3.95 m. Use `aero_database_check.yaml` for a 0.4 s untrimmed response and `pt1_hover.yaml` for the KDE static-lift setup. Neither is a validated equilibrium/flight or transition model. Longer untrimmed runs with this CG depart rapidly and can approach the existing Euler-angle singularity; the diagnostic scenario is deliberately limited to 0.4 s.

## Selecting a version

Hardware and aerodynamic versions are independent. Reference Spearhead preserves original MAD hardware and ADB v1.1 by default; PT1 KDE defaults to the new database. Flight Deck has an Aerodynamics selector. Both simulation and stability CLIs accept `--aerodynamics adb_v1_1` or `--aerodynamics tail_update_2026_10_07`. Scenario YAML accepts `aerodynamics`; explicit CLI/API selection wins. Original pitch-doublet scenario explicitly retains the old database.

The control revision is identified as `tail_update_2026_10_07_controls_2026_10_09` within the existing October dataset selector. Its path and hashes are recorded in `source.json`. Future deliveries should add new versioned configurations and preserved sources. The current loader rejects missing zero baselines, incomplete independent-surface grids, duplicate cases and combined-surface cases; additional beta planes or interaction data require extending the interpolation model, not silently flattening new dimensions.
