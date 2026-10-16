# Validation summary

**15 October 2026 · provisional engineering baseline.** What has been checked, with which
result, and what is not validated. Engineering conclusions for PT1.5 are in
[Flight Dynamics note 0003](../../../../docs/information-note/phase-1/Flight-Dynamics/0003-PT1-5-Hover-and-Fixed-Wing-Assessment/information-note.md);
the model itself is described in the [model reference](model_reference.md).

## Configuration

- ArduPlane 4.7.1 (`dbe79216`) in JSON SITL lockstep with the FDM; normal arming, frozen gains
  (PT1 logged VTOL gains, 9 October fixed-wing proposal), Ankara home (1125 m above sea level).
- PT1.5 (`pt15`): 21 kg, KDE lift rotors at the 14 October drawing stations, MAD 150KV pusher
  estimate (KDE5215XF-220 pusher as an option) with the 20×12 propeller, proposed ±30° surface
  travel, October AeroDB at 1.021 kg/m³.
- Simulation CG (−0.040, 0, +0.065) m: the target inside the 22 m/s window, not a measurement.
  PT1.5 without ballast is estimated at about −0.08 m.

## Trim at the target CG

| TAS | Alpha | Elevator | Pusher command | Pusher current | Max real part, longitudinal |
|---|---:|---:|---:|---:|---:|
| 20 m/s | 2.54° | −21.1° | 77 % | 15.0 A | −0.052 |
| 22 m/s | 1.49° | −18.2° | 83 % | 18.0 A | −0.062 |
| 24 m/s | 0.70° | −15.8° | 90 % | 21.8 A | −0.073 |
| 26 m/s | 0.08° | −14.6° | 96 % | 26.3 A | −0.090 |

All four meet the six-axis trim and level-flight residual checks. Elevator trim is beyond the
measured ±15° at every stable CG (linear extrapolation until the AQ-12 cases arrive), and above
about 24 m/s the MAD pusher exceeds its 24 A sustainable current.

With the KDE5215XF-220 pusher the attitude and elevator are the same; the command and current
change:

| TAS | Pusher command | Pusher current |
|---|---:|---:|
| 20 m/s | 62 % | 15.5 A |
| 22 m/s | 69 % | 18.7 A |
| 24 m/s | 76 % | 22.9 A |
| 26 m/s | 85 % | 28.0 A |
| 28 m/s | 95 % | 33.7 A |

The KDE option trims to 28 m/s, where the MAD option has no trim, and stays below the 44 A
180 s motor rating; KDE publishes no continuous rating, and it draws slightly more current than
the MAD motor for the same thrust.

## Checks

| Check | Result |
|---|---|
| Model tests | 568 passed, from the repository root and from the FDM folder |
| Fresh setup | The [setup guide](setup_guide.md) followed from a clean clone on Ubuntu 24.04 / Python 3.13: packages, firmware build, tests, a 60 s hover and the Flight Deck page |
| Default model on `main` | The default 20 s open-loop run, `run_trim.py` and the pitch-doublet scenario give results identical to `main` |
| Moment reference | Restored October base matches an independent wrench transfer for all 32 761 rows within 6×10⁻¹⁴ N·m |
| Hover hold, calm | Passed: 60 s free flight, max tilt 0.61°, settled altitude error 0.39 m; QHOVER drift 71 m (no position hold) |
| Hover hold, 6 m/s wind + 3 m/s gust | Passed: max tilt 0.8°, altitude error 1.0 m |
| Hover → FBWA, calm | Passed: 30.4 s sustained FBWA without lift, 24.2 m/s, peak pitch 2.1°, altitude loss 3.0 m |
| Hover → FBWA, 6 m/s headwind | Passed on airspeed: 32.5 s sustained FBWA, 24.2 m/s |
| KDE pusher option | Hover → FBWA passed (33.6 s sustained FBWA, 26.4 m/s, pusher up to 35.6 A); 60 s FBWA hold at 22 m/s passed (speed within 0.04 m/s, altitude within 0.3 m, 69 % and 18.7 A) |
| Hover doublets (pitch, roll, yaw, collective) | All passed; pitch/roll well damped, yaw sluggish and lightly damped (note 0003) |
| 22 m/s doublets (elevator, aileron, rudder, pusher) | All passed; Dutch roll period 2.50–2.53 s and damping ratio 0.22–0.24 against the linear 0.41 Hz / 0.24 |
| Front-right lift motor out at 20 s | Left the attitude envelope 0.4 s later, as expected for a quad-X |
| Lift-pack telemetry | The firmware read 46.6 V and 83 A in hover with a temporary `BATT_MONITOR = 4` overlay |
| Speed | About 0.3× real time on a two-core cloud machine running several simulations or a build at once |

Rows without a pusher named use the MAD pusher. Recordings made before 14 October used the former rotor rectangle and CG and do not validate
the current geometry. Forward → hover and ground-start runs were not repeated after the drawing
update.

## Physical consistency review

Checked against the September aerodynamic questions, the delivery specification and the
30 September CFD request:

- Control signs match the agreed convention, and the firmware's elevator (reversed) and rudder
  mappings match them.
- The base table is antisymmetric in sideslip (stable Cn_beta and Cl_beta), so zeroing Fy, Mx and
  Mz at beta = 0 adds no kink.
- Control-case drag stays near CD 0.03, so the 30 N actuator disk's direct thrust is not in the
  integrated loads and the pusher thrust is not counted twice.
- Open: the speed dependence of the neutral cases (AQ-18), lateral loads in the elevator
  increments, and the 0.97 N·m neutral pitching-moment residual. See the
  [model reference](model_reference.md#known-data-issues).
- The reference design is statically unstable with ADB v1.1 at its stated CG (static margin
  about −17 %); its fixed-wing SITL configurations use a CG 0.45 m below the origin and are for
  software regression only.

## Not validated

- Measured PT1.5 mass, CG and inertia; the pusher motor, propeller and mount height; servo and PWM
  calibration and the ±30° travel; the V-tail mixer outputs.
- Sparse control and sideslip coverage, provisional rate damping, the low-speed aerodynamic data.
- Rotor inflow, propwash, ground effect, actuator dynamics, battery sag on thrust, detailed wind.
- Any physical flight envelope. Simulation passes support test planning; they are not flight
  clearance.
