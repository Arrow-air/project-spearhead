# Aircraft profiles and propulsion data

The simulation, realtime stepping, GUI and stability tools share one engine.
Hardware differences live in `aircraft_configs/reference_design.yaml` and
`aircraft_configs/pt1_as_built.yaml`. Scenario YAML accepts `aircraft_config:`;
legacy `aircraft_case: nominal` retains its original parameters. An explicit
`params` argument or CLI `--aircraft` takes precedence over the scenario profile.
All profile-relative file paths resolve from the profile directory.

## Hardware baseline

| Item | Reference Spearhead | PT1 KDE |
| --- | --- | --- |
| Lift motors, four | MAD V8013 PRO IPE 150KV | KDE5215XF-220, 220 RPM/V |
| Lift ESC | MAD AMPX 80A V2 DroneCAN | KDE-UAS75HVC+ (product ID KDE-UAS75HVC) |
| Lift propeller | MAD FLUXER PRO 26×7.8, two blades | KDE-CF215-DP 21.5×7.3, two blades |
| Battery | 12S, 44.4 V nominal, 22 Ah | 12S, 44.4 V nominal, 22 Ah |
| Pusher | Existing 120 N × throttle placeholder | Absent |
| Lift curve | Pending manufacturer/bench data | Official KDE static data, provisional for PT1 |

The planned hardware is preserved from `docs/pt1-parts-list.md`; the KDE hardware
request supersedes those lift components for PT1. The 12S/two-blade/pusher-absent
setup follows `artifacts/firmware/spearhead-kde-first-hover/README.md`. Confirm
motor edition and the installed hardware against the real prototype.

Aircraft mass (25 kg), CG, inertia, full reference geometry and nominal motor
coordinates are inherited engineering assumptions until measured PT1 values
arrive. Component masses are recorded as hardware metadata and are not silently
added to the total mass. This avoids double counting and inventing inertia/CG.
This profile is not a reconstructed wingless vehicle or a flight controller.

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
rectification deactivated, rather than the PT1 UAS75HVC+. The resulting map is
a static surrogate until the exact installed combination is bench validated.
Some source points exceed the motor's published 44 A / 1955 W 180-second rating;
the full test map is not a continuous-duty operating envelope. Ratings are
recorded in the profile, not converted into an invented thermal controller.

## Lift forces, moments and scenarios

`Control.throttle` continues to command the forward pusher.
`Control.lift_throttle` commands all four lift rotors collectively, and optional
`lift_motor_throttles` provides one normalized command per rotor. A lift rotor
acts along body **−z**. Its moment is `(rotor_position − CG) × force`, plus the
signed reaction torque from the table. Motor positions and reaction signs remain
nominal assumptions; their list order is front-right, rear-right, rear-left,
front-left. They are not claimed to be flight-controller output ordering.

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

Reference Spearhead retains **ADB v1.1** and its **0.28 m tail chord**. PT1 KDE
defaults to **tail_update_2026_10_07**, using the delivered alpha/beta table and
17 measured deflection cases. Numeric updated tail dimensions are pending.
The wing reference chord remains 0.457 m. Both datasets are selectable with
either hardware profile through GUI, CLI and API. See the
[comparison and integration report](aerodb_tail_update_2026_10_07.md) for coverage,
the corrected moment-center metadata and the unresolved export discrepancy.

Keep the old CSV and versioned `aerodynamic_configs/`. New deliveries should
add a dataset configuration with paths, geometry version and supplied geometry
values. The current deflection loader interpolates independent surface increments
within one beta plane; more beta planes or combined-surface interactions need
an explicit model extension. Preserve the database's
reference-area/chord/moment-centre metadata rather than substituting tail chord
for the aircraft reference chord. Database-provided control effects must replace
the corresponding provisional increments to avoid counting them twice.

Before adopting the new baseline, validate coefficient curves against Aero
source data, zero-deflection differences, trim, static derivatives, eigenvalues,
CG sensitivity and control effectiveness. Zero-deflection differences may be
intended geometry effects from the changed tail, so equality with the old
database is not the acceptance criterion.
