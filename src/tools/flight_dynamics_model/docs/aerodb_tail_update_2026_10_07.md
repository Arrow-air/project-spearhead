# October 2026 AeroDB comparison and integration

Dated validation reports and recordings referenced below are retained locally and excluded from Git.

**Current control revision: 11 October 2026.** `deflection-3.csv` supersedes the
49-case control revision discussed below. Its preserved source and explicit sparse
lookup sidecar are versioned in `data/aerodynamics/tail_update_2026_10_07/`.
See the new intake, transition and flight mechanics review (local report).
The historical 9 October selector is `tail_update_2026_10_07_controls_2026_10_09`.
The following numerical comparisons describe their original deliveries and historical intake interpretation.

**14 October confirmation:** the Aero team confirmed that the Nondimit export center of `database.csv` was (-0.1,0,+0.15) and that the deflection tables are referenced to (0,0,0). The base table carries the 1.3 wind-axis drag multiplier and the deflection tables do not (Fx at alpha=0: -24.5 vs -18.9 N). Nondimit applies force corrections before the moment-center shift (Nondimit README, "Force and Moment Corrections"), so restoring `database.csv` to (0,0,0) with its own corrected forces is an exact inverse and the multiplier adds no pitching moment at the origin; the extra drag is in effect applied through (0,0,0). The 0.97 N·m neutral residual at alpha=0 is therefore not explained by the multiplier; base-surrogate versus control-run differences at other alphas are 0.5 to 3.3 N·m (alpha=-5, 5, 10°).

**12 October correction:** Nondimit writes the final processed moment center in the export header. The active base is now numerically restored from (-0.1,0,+0.15) to the user-confirmed raw CFD origin (0,0,0). Earlier header-only comparisons and flight results are historical; see the [current validation summary](validation_summary.md).

Received `database.csv` and `deflection.csv` on 2026-10-07. The user confirms raw CFD body-axis moment center (0,0,0) m and density 1.021 kg/m³. The delivered base is a processed Nondimit export at (-0.1,0,+0.15); the control loads inherit the raw origin. Original source bytes are preserved in `sources/`, with SHA-256 hashes in `source.json`. The active base changes dimensional moments and derived moment coefficients by the documented reference transfer; dimensional forces are unchanged, and already-applied drag correction/zeroing is not repeated.

The later `results.csv` control expansion is integrated as of 9 October. See the update assessment and multi-speed validation (local report) for current results, provenance and limitations.

## Historical quick comparison

Both base databases contain the complete 181 × 181 alpha/beta grid, −90..90° at 1° steps. The table compares coefficients directly. This historical table treated the new moments as origin-referenced without restoring the processed export first. That assumption was corrected on 12 October. Both datasets use density 1.021 kg/m³ and speed 25 m/s here; the old moments are transferred from their documented (−0.15,0,0.15) m center. Density alone does not explain coefficient differences. Updated tail geometry is expected to change loads, but these files do not isolate its effect from export/correction differences.

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

## Processed export reference

The base header records a final Nondimit export center (−0.1,0,+0.15). Native Nondimit documentation and source confirm that its header describes the output reference. Restoring the same loads to the raw CFD origin gives `M_origin = M_export + r_export × F`. At V=25 m/s and α=β=0, My changes from −13.2069 to −37.540345 N·m; the control neutral is −36.5681 N·m. The residual is about 0.972 N·m, rather than 23.36 N·m. A reference shift explains most of this discrepancy, but the original CFD/postprocessing history is not independently available and the remaining difference is unresolved. SakeDB supplies dimensional loads; the inspected template does not perform a CG-reference transfer. Baseline subtraction still prevents adding a second neutral base moment. `scripts/recenter_october_database.py` reproduces the 32761-row conversion and verifies it independently.

## Trim and geometry status

The original 17-case intake study at V=25 m/s, mass=25 kg and CG=(−0.15,0,0.15) m gave a bounded preliminary trim at alpha=-0.557°, elevator=5.000°. It did **not** satisfy force/moment balance within measured elevator bounds ±5°. Residual [Fx/m, Fz/m, My/Iyy, theta−alpha] = [0.0, 1.826671, 1.789204, 0.358784]. This remains historical intake evidence. The later update assessment documents successful diagnostic-reference trims at 18–28 m/s using a different CG and propulsion setup, including elevator extrapolation at 18 and 20 m/s. A solver termination flag alone is insufficient for trim success; initialization checks actual equilibrium residuals.

The numeric updated tail chord/area were not supplied; they are marked pending instead of assigning the old 0.28 m value to the new geometry. Wing reference chord remains 0.457 m, area 1.805 m² and span 3.95 m. Use `aero_database_check.yaml` for a 0.4 s untrimmed response and `pt1_hover.yaml` for the KDE static-lift setup. Neither is a validated equilibrium/flight or transition model. Longer untrimmed runs with this CG depart rapidly and can approach the existing Euler-angle singularity; the diagnostic scenario is deliberately limited to 0.4 s.

## Selecting a version

Hardware and aerodynamic versions are independent. Reference Spearhead preserves original MAD hardware and ADB v1.1 by default; PT1 KDE defaults to the new database. Flight Deck has an Aerodynamics selector. Both simulation and stability CLIs accept `--aerodynamics adb_v1_1` or `--aerodynamics tail_update_2026_10_07`. Scenario YAML accepts `aerodynamics`; explicit CLI/API selection wins. Original pitch-doublet scenario explicitly retains the old database.

The control revision is identified as `tail_update_2026_10_07_controls_2026_10_09` within the existing October dataset selector. Its path and hashes are recorded in `source.json`. Future deliveries should add new versioned configurations and preserved sources. The current loader rejects missing zero baselines, incomplete independent-surface grids, duplicate cases and combined-surface cases; additional beta planes or interaction data require extending the interpolation model, not silently flattening new dimensions.
