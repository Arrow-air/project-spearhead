# Flight Deck validation summary

**14 October 2026 · Provisional engineering baseline.** This is the single public validation report. Recordings, logs and dated diagnostic folders remain local and excluded from Git.

## Configuration and moment reference

Native **ArduPlane 4.7.1** (`dbe792162d06cab66c3475fd5556bf7a120f119e`) controls the Spearhead 6-DOF FDM. Runs use normal arming, verified parameters and fresh mode evidence. PT1.5 uses 21 kg, KDE lift hardware, an estimated MAD 150KV/APC 20×12 pusher and proposed ±30° surface travel. SITL home is the Ankara test field (1125 m above sea level). Controller gains and flight acceptance limits remain frozen.

The aero team confirms raw CFD origin `(0,0,0)` and density 1.021 kg/m³. The delivered Nondimit base export records its **processed** center `(-0.1,0,+0.15)`. On 12 October the active table was numerically restored to the raw origin, replacing the former header-only handling. Dimensional forces are unchanged; previously applied corrections are not repeated. All 32761 rows match an independent wrench transfer within 5.7×10⁻¹⁴ N·m. The preserved delivery and reproducible conversion are hashed in `source.json`.

Runtime then applies `M_CG = M_source + (r_source − r_CG) × F` once, after adding origin-referenced control increments. At V=25 m/s, α=β=0, corrected base My is −37.540345 N·m versus control-neutral −36.5681 N·m. The Aero team confirmed the (−0.1,0,+0.15) export center on 14 October. Nondimit applies the 1.3 drag multiplier before the moment shift, so the restore is exact and the multiplier does not cause the remaining 0.972 N·m; base-surrogate and control-run differences at other alphas are 0.5–3.3 N·m.

## Geometry, CG and trim (14 October 2026)

Body axes are x forward, y right, z down, from the CFD origin (wing root quarter-chord). Lift-rotor and pusher stations follow the 14 October drawing: rotors at x = +0.675/−0.750 m, y = ±0.55 m, z = −0.043 m (centroid x = −0.0375 m); pusher at (−0.374, 0, +0.050) m, the lowest of the open mount heights. The former nominal rectangle (+0.55/−0.85 m) is withdrawn.

PT1.5 (`pt15`, 21 kg) uses the target CG **[−0.040, 0, +0.065] m**: inside the 22 m/s window (static margin about 6 %). The PT1.5 build-up estimate without ballast is about −0.08 m, so the built aircraft must be weighed and ballasted to the target. PT1 as flown (`pt1_as_built`, 17.9 kg, wings on, no pusher) uses x = −0.0475 m from the logged 49.3 % front hover share; separate hover segments vary by about ±10 mm. CG z (+0.065 m), inertia and rate damping are build-up and vortex-lattice estimates.

| TAS | Alpha | Elevator | Pusher command | Estimated pusher current | Max real part, longitudinal |
|---|---:|---:|---:|---:|---:|
| 20 m/s | 2.541° | −21.085° | 77.23 % | 15.03 A | −0.052 |
| 22 m/s | 1.494° | −18.206° | 83.09 % | 18.00 A | −0.062 |
| 24 m/s | 0.699° | −15.824° | 89.51 % | 21.78 A | −0.073 |
| 26 m/s | 0.079° | −14.552° | 96.20 % | 26.31 A | −0.090 |

All four satisfy six-axis trim and level-flight residual checks. Elevator trim at every stable CG lies beyond the measured ±15° coverage, so it rests on linear extrapolation until the AQ-12 cases (−20°, −30°) arrive. Above about 24 m/s the pusher exceeds its 24 A sustainable-current guideline. Forward trim now tries several elevator starting points; the former single start failed between solvable neighboring points.

## Current verification

| Check | Result and limits |
|---|---|
| Source-only regression | **539 passed**, from the repository root |
| PT1.5 hover hold, target CG | **Passed:** normal arm; 60.005 s free flight; max tilt 0.61°; settled altitude error 0.39 m; QHOVER drifted 71 m horizontally without position hold |
| PT1.5 hover → forward, 9 October parameters | **Passed:** 60.010 s free flight; 30.35 s sustained accepted FBWA with lift idle; final TAS 24.16 m/s; peak pitch 2.11°; altitude loss 2.99 m |
| PT1.5 hover → forward, PT1.5 parameter set | **Passed:** 60.010 s free flight; 26.06 s sustained accepted FBWA; final TAS 24.20 m/s; peak pitch 0.75°; altitude loss 0.94 m |
| Lift-battery telemetry | ArduPilot read the JSON pack voltage and current (46.6 V, 83.2 A in hover) with a temporary `BATT_MONITOR = 4` overlay; the repository SITL parameters still have the monitor off |
| Performance | About 340–445 s wall time for 107 s of simulation with three runs sharing two CPU cores |

Earlier recordings used the former rotor rectangle, a +10 mm CG and z = +150 mm; they do not validate the current geometry. Forward→hover and doublet runs were not repeated in this update.

## Remaining limits

Measured PT1.5 CG and inertia, the pusher mount height and motor/propeller, and servo/PWM calibration remain open. Sparse control/sideslip coverage, provisional rate damping, absent rotor inflow/propwash/transients and a generic battery discharge curve (telemetry only; thrust uses nominal voltage) limit predictions. Physical placement must also account for the wing, not merely the sign relative to an arbitrary origin.

The display retains the disabled 3D implementation, numeric FL/FR/RL/RR lift values and logged 1100–1940 µs limits. Documentation diagrams and local links are checked. Broader formatting/spelling findings remain; provisioned Docker checks were not run because its daemon was unavailable. This is a reviewable provisional simulation baseline, not hardware flight readiness. See the [capabilities](capabilities.md), [model reference](model_reference.md) and [setup README](../README.md).
