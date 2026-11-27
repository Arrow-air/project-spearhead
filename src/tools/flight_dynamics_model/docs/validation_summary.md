# Validation summary

**26 November 2026 · provisional engineering baseline.** What has been checked, with which
result, and what is not validated. The PT1.5 hardware changed on 21 November (MAD lift units,
HOBBYWING pusher) and the AeroDB reference density was corrected on 26 November; see the
[hardware update](update_2026_11_21.md). Engineering conclusions for the earlier 21 kg KDE
configuration are in
[Flight Dynamics note 0003](../../../../docs/information-note/phase-1/Flight-Dynamics/0003-PT1-5-Hover-and-Fixed-Wing-Assessment/information-note.md);
the model itself is described in the [model reference](model_reference.md).

## Configuration

- ArduPlane 4.7.1 (`dbe79216`) in JSON SITL lockstep with the FDM; normal arming, frozen gains
  (PT1 logged VTOL gains, 9 October fixed-wing proposal), Ankara home (1125 m above sea level).
- PT1.5 (`pt15`): 22 kg; MAD V8013 150KV / AMPX 80A / FLUXER PRO 26×7.8 lift units at the
  14 October drawing stations; HOBBYWING 5230SL 380KV / AMPX 80A / APC 20×12 pusher limited to
  45 % of the ESC range (`SERVO3_MAX` 1450 µs); proposed ±30° surface travel.
- October AeroDB with reference density 1.092 kg/m³ (CFD run condition), flown at 1.021 kg/m³.
- Simulation CG (−0.040, 0, +0.065) m: the target inside the 22 m/s window, not a measurement.
  PT1.5 without ballast is estimated at about −0.085 m.

## Trim at the target CG

| TAS | Alpha | Elevator | Pusher command | Pusher current | Static margin | Phugoid real part |
|---|---:|---:|---:|---:|---:|---:|
| 20 m/s | 3.27° | −21.7° | 70 % | 15.5 A | 4.4 % | −0.044 |
| 22 m/s | 2.10° | −18.8° | 75 % | 18.3 A | 5.5 % | −0.053 |
| 24 m/s | 1.20° | −16.4° | 81 % | 21.8 A | 6.9 % | −0.062 |
| 26 m/s | 0.51° | −15.0° | 87 % | 26.0 A | 8.3 % | −0.077 |
| 28 m/s | −0.04° | −14.4° | 93 % | 31.1 A | 9.8 % | −0.081 |

All five meet the six-axis trim and level-flight residual checks. The pusher command is a
percentage of the limited range (100 % = 1450 µs = 45 % of the ESC range). Elevator trim is beyond
the measured ±15° up to 24 m/s and at the limit at 26 m/s (linear extrapolation until the AQ-12 cases arrive). The spiral
mode stays mildly divergent (+0.06 to +0.09 1/s).

The density correction lowers the aerodynamic loads by 6.5 %, so every trim needs more alpha.
The static margin falls with alpha in this database, so the 22 m/s window (margin at least 5 %,
elevator no lower than −20°) narrows from −0.046…−0.027 m to about −0.042…−0.033 m. The
−0.040 m target is still inside it. The untrimmed CLmax of the 25 m/s table is 1.30, a stall
speed of 13.4 m/s at 22 kg (12.9 m/s before the correction).

## Checks, 26 November (new hardware, corrected density)

| Check | Result |
|---|---|
| Model tests | 572 passed (1 skipped: PT1 DataFlash sample not in the repository) |
| Corrected AeroDB | Coefficients from the preserved dimensional loads at 1.092 kg/m³ match the revised note 0003 table within its 10⁻⁴ rounding; dimensional loads unchanged |
| Lift curve against the MAD sheet | Thrust within 5 %, current within 7.5 %, RPM within 2 % at 48 V and sea level |
| Pusher fit against the HOBBYWING sheet | Duty within 4.8 % of throttle, battery current RMS 2.3 %, up to 80 A |
| Hover hold, calm | Passed: 60 s free flight, max tilt 0.58°, settled altitude error 0.68 m, lift PWM 1630–1656 µs; QHOVER drift 73 m (no position hold) |
| Hover hold, 6 m/s headwind + 3 m/s gust | Passed: max tilt 0.64°, altitude error 0.87 m, lift PWM 1593–1648 µs |
| Front-right lift motor out at 20 s | Left the attitude envelope 0.46 s after the failure, as expected for a quad-X |
| Ground start (takeoff, 10 m, 20 s hover, landing) | Passed: normal arming, max tilt 0.31°, landed on the gear and disarmed |
| Hover → FBWA, calm | Passed: 32.2 s sustained FBWA without lift, 26.3 m/s, peak pitch 1.7°, altitude loss 0.4 m; pusher up to 90 % of the limited range, 33 A and 55 N |
| Hover → FBWA, 6 m/s headwind | Passed: 34.2 s sustained FBWA, 26.3 m/s, peak pitch 1.7°, no altitude loss |
| 22 m/s FBWA hold, 60 s | Passed: speed 22.00–22.02 m/s, altitude within 0.05 m, pusher 76 % and 18.3 A |
| Forward → hover (QLOITER), 9 October parameters | Passed: 38.8 s sustained hover, altitude loss 0.4 m, peak pitch 20.9° during the deceleration |
| Hover doublets, ±8° pitch / roll | Passed: rise to 90 % in 0.60 / 0.61 s, overshoot on the reversal 16 % / 13 % (KDE: 36 % / 24 %), within 1° 0.6 s after neutral |
| Hover doublet, yaw stick ±20 % | Passed: 20.6 / −20.8°/s achieved against 38 / −51°/s demanded (KDE: about 15°/s), lag 1.1 s; heading oscillation period 5.0 s, damping ratio 0.11, within 2° after 15 s (KDE: 28 s) |
| Hover doublet, collective ±10 points | Passed: +1.8 / −2.5 m/s², climb 3.3 m/s, QHOVER holds the new height (+2.9 m) |
| 22 m/s doublets, elevator ±3° | Passed: q 9.2°/s (2.3°/s per degree), pitch 4.6°, load factor 0.34 g; rates within 1°/s 1.9 s after the input |
| 22 m/s doublets, aileron ±5° / rudder ±5° | Passed: Dutch roll period 2.75 / 2.70 s, damping ratio 0.24 / 0.22; adverse yaw, sideslip 5.1° from aileron; rudder trim 2.2° |
| 22 m/s doublet, pusher ±10 points | Passed: speed 0.8 m/s, height 2.1 m; pusher up to 27 A |

No run reached the pusher output limit; the highest command was 90 % of the limited range, in
the transitions. The 45 % limit is unchanged.

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
