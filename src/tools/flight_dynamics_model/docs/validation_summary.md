# Flight Deck validation summary

**12 October 2026 · Provisional engineering baseline.** This is the single public validation report. Recordings, logs and dated diagnostic folders remain local and excluded from Git.

## Configuration and moment reference

Native **ArduPlane 4.7.1** (`dbe792162d06cab66c3475fd5556bf7a120f119e`) controls the Spearhead 6-DOF FDM. Runs use normal arming, verified parameters and fresh mode evidence. PT1.5 uses 21 kg, KDE lift hardware, an estimated MAD 150KV/APC 20×12 pusher and proposed ±30° surface travel. Controller gains and flight acceptance limits remain frozen.

The aero team confirms raw CFD origin `(0,0,0)` and density 1.021 kg/m³. The delivered Nondimit base export records its **processed** center `(-0.1,0,+0.15)`. On 12 October the active table was numerically restored to the raw origin, replacing the former header-only handling. Dimensional forces are unchanged; previously applied corrections are not repeated. All 32761 rows match an independent wrench transfer within 5.7×10⁻¹⁴ N·m. The preserved delivery and reproducible conversion are hashed in `source.json`.

Runtime then applies `M_CG = M_source + (r_source − r_CG) × F` once, after adding origin-referenced control increments. At V=25 m/s, α=β=0, corrected base My is −37.540345 N·m versus control-neutral −36.5681 N·m. The remaining 0.972 N·m difference is unresolved; the export-center transfer explains most of the previous 23.36 N·m discrepancy. Original postprocessing history is not independently confirmed.

## Forward CG and trim

Body axes are x forward, y right, z down. The provisional default is **[+0.010,0,+0.150] m**, 10 mm forward of the October origin. All prototype starts use this selection. Hover initialization solves independent rotor loads; actual flight motor outputs still come from ArduPilot.

The [preliminary drawing](../../../../docs/information-note/phase-1/Aerodynamics/0001-Preliminary-Design/information-note.md) shows a 275 mm tall high-wing fuselage, but cannot establish the October datum or mass CG. The existing +150 mm z offset is retained explicitly as an unmeasured assumption. A +80 mm sensitivity case also trims; it is not a CAD-derived CG measurement.

| TAS | Alpha | Elevator | Pusher command | Estimated current |
|---|---:|---:|---:|---:|
| 22 m/s | 1.653° | −24.796° | 84.21% | 18.82 A |
| 24 m/s | 0.826° | −20.552° | 90.73% | 22.80 A |
| 26 m/s | 0.181° | −18.092° | 97.29% | 27.34 A |

All three satisfy six-axis trim and level-flight residual checks. Fixed-control Cmα is −1.01 to −1.09 per radian and the four local longitudinal eigenvalues have negative real parts. This supports local model stability, not a physical flight envelope. Elevator trim uses substantial extrapolation beyond measured samples; 26 m/s has little throttle margin and exceeds the estimated 24 A sustainable-current guideline. +25 mm leaves less elevator margin; +50 mm fails the 22 m/s trim within existing bounds. +100 mm failed the sampled 24 m/s balance checks.

## Current verification

| Check | Result and limits |
|---|---|
| Source-only regression | **521 passed in 93.60 s**, from repository root with both artifact folders absent |
| +10 mm hover → forward | **Passed:** normal arm; 60.010 s free flight; 32.455 s sustained accepted FBWA with lift idle; final TAS 23.813 m/s; peak pitch 14.951°; altitude loss 3.636 m |
| +10 mm forward → hover | **Passed:** 60.015 s free flight; 37.725 s accepted QLOITER hover; final TAS 0.014 m/s; pusher off; braking peak pitch 20.182° |
| Startup and performance | Frozen 9 October proposal; initial hover pitch transient reached 14.951° and firmware reported potential VTOL thrust loss. 158.94 s wall time for 107.13 s total simulation; real-time performance is not established |
| Corrected legacy aft CG | Failed: 60.836° peak pitch, only 2.715 s accepted forward interval; stopped by the existing envelope |

Earlier recordings used the former reference handling. The two transition directions above are newly verified; the eight earlier doublet passes remain historical and do not validate this corrected model. Forward and hover doublets are available; their input delivery, recovery and actuator headroom are checked separately.

## Remaining limits

CG/datum, inertia, final motor coordinates, thrust line and servo/PWM calibration remain unmeasured. Sparse control/sideslip coverage, provisional rate damping, absent rotor inflow/propwash/transients and unvalidated battery behavior limit predictions. Physical placement must also account for the wing, not merely the sign relative to an arbitrary origin.

The display retains the disabled 3D implementation, numeric FL/FR/RL/RR lift values and logged 1100–1940 µs limits. Documentation diagrams and local links are checked. Broader formatting/spelling findings remain; provisioned Docker checks were not run because its daemon was unavailable. This is a reviewable provisional simulation baseline, not hardware flight readiness. See the [capabilities](capabilities.md), [model reference](model_reference.md) and [setup README](../README.md).
