# Aircraft profiles and propulsion data

The simulation, realtime stepping, GUI and stability tools share one engine.
Hardware differences live in `aircraft_configs/reference_design.yaml`,
`aircraft_configs/pt1_as_built.yaml` (PT1 as flown, no pusher) and
`aircraft_configs/pt15.yaml` (PT1.5 with pusher; used by SITL and the Flight Deck).
Scenario YAML accepts `aircraft_config:`;
legacy `aircraft_case: nominal` retains its original parameters. An explicit
`params` argument or CLI `--aircraft` takes precedence over the scenario profile.
All profile-relative file paths resolve from the profile directory.

## Hardware baseline

| Item | Reference Spearhead | PT1.5 KDE |
| --- | --- | --- |
| Lift motors, four | MAD V8013 PRO IPE 150KV | KDE5215XF-220, 220 RPM/V |
| Lift ESC | MAD AMPX 80A V2 DroneCAN | KDE-UAS75HVC+ (product ID KDE-UAS75HVC) |
| Lift propeller | MAD FLUXER PRO 26×7.8, two blades | KDE-CF215-DP 21.5×7.3, two blades |
| Battery | 12S, 44.4 V nominal, 22 Ah | 12S, 44.4 V nominal, 22 Ah |
| Pusher | Existing 120 N × throttle placeholder | MAD V8013 150KV / AMPX 80A / APC Sport 20×12; provisional speed-dependent estimate |
| Lift curve | Pending manufacturer/bench data | Official KDE static data, provisional for PT1.5 |

The planned hardware is preserved from `docs/pt1-parts-list.md`; the KDE hardware
request supersedes those lift components for PT1.5. The original 12S/two-blade/pusher-absent
first-hover setup is retained as local development history; its diagnostic report is excluded from Git.
The user added the electric pusher and confirmed 150 KV on 11 October 2026. Confirm
motor edition and the installed hardware against the real prototype.

## 14 October 2026 drawing update

Lift-rotor and pusher stations come from the 14 October drawing table (origin at the
CFD origin, X aft, Z up). In body axes (x forward, z down) the rotors are at
x = +0.675 m (front) and −0.750 m (rear), y = ±0.55 m, z = −0.043 m; their centroid is
x = −0.0375 m. The pusher is at x = −0.374 m; its height is open between Z = +40 and
−50 mm, and `pt15` uses the lowest option (body z = +0.050 m), closest to the estimated
CG height. The former nominal rectangle (+0.55/−0.85 m, z = +0.15 m) was a placeholder
and also set the logged-hover CG interpretation; it is withdrawn.

| Item | `pt1_as_built` (PT1 as flown) | `pt15` (PT1.5) |
| --- | --- | --- |
| Mass | 17.9 kg (log-derived, wings on) | 21.0 kg build estimate |
| CG x | −0.0475 m from the logged 49.3 % front hover share (±10 mm between hover segments) | −0.040 m target inside the 22 m/s window; the build-up estimate without ballast is about −0.08 m |
| CG z | +0.065 m build-up estimate | +0.065 m build-up estimate |
| Inertia Ixx/Iyy/Izz/Ixz (kg·m²) | 4.30 / 3.73 / 7.66 / −0.29 | 4.42 / 4.07 / 8.11 / −0.31 |
| Pusher | none | MAD V8013 150KV / AMPX 80A / APC 20×12 estimate at (−0.374, 0, +0.050) m |

Inertia is a component build-up (±20 %). Profiles may set
`aerodynamics.rate_derivatives` (`Clp`, `Clr`, `Cmq`, `Cnp`, `Cnr`); both prototype
profiles use wing + inverse-V-tail vortex-lattice values without fuselage or booms:
Cmq −9.6/−9.7, Clp −0.58, Clr +0.10, Cnp −0.04, Cnr −0.04. The low Cnr is the
lifting-surface value only; the earlier −0.25 had no source.

At the PT1.5 target CG the October AeroDB gives a static margin of about 6 % and
elevator trim of about −18° at 22 m/s (−21° at 20 m/s). That is beyond the measured
±15° elevator coverage; the AQ-12 cases (−20°, −30°) are needed to confirm it.

PT1.5 mass is 21 kg and the planned finished reference is 25 kg, as supplied by
the user. The PT1.5 CG must be set by weighing and ballasting the built aircraft;
the profile value is the simulation target, not a measurement.
Component masses are recorded as hardware metadata and are not silently
added to the total mass. This avoids double counting and inventing inertia/CG.
This profile is not a reconstructed wingless vehicle or a flight controller.

## PT1.5 electric pusher update — 11 October 2026

`pt1_as_built.yaml` enables the user-selected MAD V8013 PRO IPE **150 KV** motor,
AMPX 80A (5–14S) ESC and **APC LP20012 Sport 20×12** propeller. `pusher.py`
solves quasi-steady motor RPM using a fit to all 15 official motor test-sheet
RPM/torque/throttle points and APC's exact analytical Ct/Cp map. It uses body
axial speed and FDM density; it does not reuse MAD 26-inch static thrust.

At 44.4 V and density 1.021, estimated full-duty thrust is 62.9 N static and
32.8 N at 26 m/s. These are short-run estimates: the corresponding currents
are 29.0 A and 30.0 A, above the sheet's 24 A sustainable-current guideline.
The operating flag does not impose an invented thermal derating model.
The pusher thrust acts at the drawing station (see the 14 October update);
rotation direction is undecided, so propeller reaction torque is omitted. The
installed motor (KDE or MAD) and the wooden 20×12 propeller are not yet characterized.

The exact manufacturer map and test-sheet image, transcription, fitted constants,
assumptions and hashes live under `data/propulsion/pt1_pusher_source.json`.
Generated reports and performance plots remain in local, ignored `artifacts/` folders.
Reproduce them with `scripts/estimate_pt1_pusher.py`. PT1.5 hover starts still keep
the pusher off; ground takeoff and transition sequences require separate validation.

## KDE source and interpolation

The official [motor page](https://www.kdedirect.com/products/kde5215xf-220),
[ESC page](https://www.kdedirect.com/products/kde-uas75hvc) and
[propeller page](https://www.kdedirect.com/products/kde-cf215-dp) supply the
hardware specifications. The motor page links the HE edition performance data.

The preserved source workbook is
`data/propulsion/sources/KDE5215-220kV_rev2.xlsx`. Its URL, SHA256, source rows,
test conditions and column mapping are in `data/propulsion/kde5215xf_220_source.json`.
The runtime CSV contains 28 samples for the requested two-blade propeller:
25.2, 33.6, 42.0 and 50.4 V, each with throttle from 25% to 100%.

The engine linearly interpolates voltage and throttle, including 44.4 V nominal
between the 42.0 and 50.4 V tests. Values outside the measured voltage envelope
are rejected. Throttle saturates at 0–1. Below 25%, an assumed linear ramp from
motor-off fills the unmeasured startup region.

Runtime thrust uses the measured **THRUST OUTPUT [N]** column. The manufacturer's
estimated sea-level thrust is preserved for reference and is not substituted
for measured thrust. No density, advance-speed/inflow, rotor wake, battery sag,
ESC loss, spool dynamics or thermal-duration corrections are applied. Electrical
readouts are bench estimates at the selected voltage and throttle.

The source test used a **KDE-UAS125UVC** on firmware **D460337**, with synchronous
rectification deactivated, rather than the PT1.5 UAS75HVC+. The resulting map is
a static surrogate until the exact installed combination is bench validated.
Some source points exceed the motor's published 44 A / 1955 W 180-second rating;
the full test map is not a continuous-duty operating envelope. Ratings are
recorded in the profile, not converted into an invented thermal controller.

## Lift forces, moments and scenarios

`Control.throttle` continues to command the forward pusher.
`Control.lift_throttle` commands all four lift rotors collectively, and optional
`lift_motor_throttles` provides one normalized command per rotor. A lift rotor
acts along body **−z**. Its moment is `(rotor_position − CG) × force`, plus the
signed reaction torque from the table. Motor positions follow the 14 October
drawing; their list order is front-right, rear-right, rear-left, front-left
(ArduPilot motors 1, 4, 2, 3).

The `pt1_hover.yaml` example starts level at 20 m and estimates collective that
supports the selected mass on the static curve, then applies a 3% collective
doublet. It is an open-loop response experiment, without an altitude/attitude
controller, ground contact or validated hover stability. It does not implement
VTOL transition. A pusher-absent profile cannot use fixed-wing trim/stability.
At rest, aerodynamic dynamic pressure is zero. The existing `Vmin` guard is
used only for rate normalization, preventing artificial aerodynamic loads from
disturbing the static lift initial condition.
The reference MAD lift system rejects nonzero commands until its own curve is
provided, rather than reusing the KDE data.

JSON exports identify the aircraft, AeroDB, geometry and hardware; CSV appends
aircraft/AeroDB identity after the original state columns. GUI selection resets
simulation history and dependent analysis results.

## AeroDB selection and changed tail chord

Reference Spearhead retains **ADB v1.1** and its **0.28 m tail chord**. PT1.5 KDE
defaults to **tail_update_2026_10_07**, using the delivered alpha/beta table and
the 103-case `deflection-3.csv` control delivery received 11 October. Sparse
beta-zero samples cover all three surfaces at 25 m/s, including wider travel,
two additional elevator speed planes and four elevator/rudder combinations.
Six neutral sideslip cases remain audit data. The 49-case 9 October revision is
preserved as `tail_update_2026_10_07_controls_2026_10_09`. Numeric updated tail dimensions are pending.
The wing reference chord remains 0.457 m. Both datasets are selectable with
either hardware profile through GUI, CLI and API. See the
[comparison and integration report](aerodb_tail_update_2026_10_07.md) for coverage,
the corrected moment-center metadata and the unresolved export discrepancy.

Keep the old CSV and versioned `aerodynamic_configs/`. New deliveries should
add a dataset configuration with paths, geometry version and supplied geometry
values. The current expanded loader triangulates sparse alpha/deflection planes,
interpolates normalized elevator coefficients across speed with explicit endpoint
assumptions, and adds a bounded residual for the measured elevator/rudder pairs.
Control beta dependence and other surface interactions remain unmeasured. Preserve the database's
reference-area/chord/moment-centre metadata rather than substituting tail chord
for the aircraft reference chord. Database-provided control effects must replace
the corresponding provisional increments to avoid counting them twice.

Before adopting the new baseline, validate coefficient curves against Aero
source data, zero-deflection differences, trim, static derivatives, eigenvalues,
CG sensitivity and control effectiveness. Zero-deflection differences may be
intended geometry effects from the changed tail, so equality with the old
database is not the acceptance criterion.
