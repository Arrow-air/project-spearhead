# Validation summary

**21 November 2026 · provisional engineering baseline.** What has been checked, with which
result, and what is not validated. The PT1.5 hardware changed on 21 November (MAD lift units,
HOBBYWING pusher; see the [hardware update](update_2026_11_21.md)). Engineering conclusions for
the earlier 21 kg KDE configuration are in
[Flight Dynamics note 0003](../../../../docs/information-note/phase-1/Flight-Dynamics/0003-PT1-5-Hover-and-Fixed-Wing-Assessment/information-note.md);
the model itself is described in the [model reference](model_reference.md).

## Configuration

- ArduPlane 4.7.1 (`dbe79216`) in JSON SITL lockstep with the FDM; normal arming, frozen gains
  (PT1 logged VTOL gains, 9 October fixed-wing proposal), Ankara home (1125 m above sea level).
- PT1.5 (`pt15`): 22 kg; MAD V8013 150KV / AMPX 80A / FLUXER PRO 26×7.8 lift units at the
  14 October drawing stations; HOBBYWING 5230SL 380KV / AMPX 80A / APC 20×12 pusher limited to
  45 % of the ESC range (`SERVO3_MAX` 1450 µs); proposed ±30° surface travel; October AeroDB at
  1.021 kg/m³.
- Simulation CG (−0.040, 0, +0.065) m: the target inside the 22 m/s window, not a measurement.
  PT1.5 without ballast is estimated at about −0.085 m.

## Trim at the target CG

| TAS | Alpha | Elevator | Pusher command | Pusher current | Phugoid real part |
|---|---:|---:|---:|---:|---:|
| 20 m/s | 2.83° | −21.3° | 71 % | 15.9 A | −0.045 |
| 22 m/s | 1.73° | −18.4° | 76 % | 18.9 A | −0.056 |
| 24 m/s | 0.90° | −16.0° | 82 % | 22.8 A | −0.067 |
| 26 m/s | 0.25° | −14.7° | 88 % | 27.4 A | −0.081 |
| 28 m/s | −0.26° | −14.2° | 94 % | 32.7 A | −0.085 |

All five meet the six-axis trim and level-flight residual checks. The pusher command is a
percentage of the limited range (100 % = 1450 µs = 45 % of the ESC range). Elevator trim is beyond
the measured ±15° up to 24 m/s (linear extrapolation until the AQ-12 cases arrive). The spiral
mode stays mildly divergent (about +0.08 1/s), as before.

## Checks, 21 November hardware

| Check | Result |
|---|---|
| Model tests | 572 passed (1 skipped: PT1 DataFlash sample not in the repository) |
| Lift curve against the MAD sheet | Thrust within 5 %, current within 7.5 %, RPM within 2 % at 48 V and sea level |
| Pusher fit against the HOBBYWING sheet | Duty within 4.8 % of throttle, battery current RMS 2.3 %, up to 80 A |
| Hover hold, calm | Passed: 60 s free flight, max tilt 0.57°, settled altitude error 0.68 m; lift PWM 1629–1651 µs; lift pack 79 A at 45.7–46.6 V; QHOVER drift 72 m (no position hold) |
| Hover → FBWA, 9 October parameters | Passed: 32.4 s sustained FBWA without lift, 26.0 m/s, peak pitch 1.4°, altitude loss 0.4 m; pusher up to 90 % of the limited range, 33 A and 55 N |
| 22 m/s FBWA hold, 60 s | Passed: speed 21.94–22.04 m/s, altitude within 0.74 m, attitude RMS error below 0.05°; pusher at 76 % (1354 µs) and 18.7 A; firmware VFR_HUD throttle 76 % |

Not yet repeated with the new hardware: wind and gust hover, hover and forward doublets, lift
motor failure, forward → hover and ground start.

## Checks, October baseline (21 kg, KDE lift; superseded hardware)

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

Rows without a pusher named used the MAD pusher. Recordings made before 14 October used the
former rotor rectangle and CG and do not validate the current geometry.

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

- Measured PT1.5 mass, CG and inertia; lift and pusher bench data with the installed propellers;
  the pusher mount height; servo and ESC PWM calibration and the ±30° travel; the V-tail mixer
  outputs.
- Sparse control and sideslip coverage, provisional rate damping, the low-speed aerodynamic data.
- Rotor inflow, propwash, ground effect, actuator dynamics, battery sag on thrust, detailed wind.
- Any physical flight envelope. Simulation passes support test planning; they are not flight
  clearance.
