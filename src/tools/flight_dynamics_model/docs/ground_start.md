# PT1 ground-start simulation

Choose **PT1 → Ground · VTOL takeoff → Takeoff → hover → land** in the Flight Deck.
Set the takeoff altitude (3–30 m) and hover time (15–600 whole seconds), then select **Start flight**.
Calibration and normal arming happen on the ground. Climb and landing are additional to the hover time.
The display reports ON GROUND, TAKEOFF, HOVER, LANDING and the completed validation result.
ArduPlane 4.7.1 owns the motor outputs throughout. This automatic test does not enable interactive mode changes.

The uploaded mission contains HOME, NAV_VTOL_TAKEOFF, DO_VTOL_TRANSITION(MC),
NAV_LOITER_TIME and NAV_VTOL_LAND. The explicit MC command keeps the aircraft in VTOL between takeoff and loiter.
Mission acceptance, firmware identity, parameter readback, EKF navigation validity and normal arming are checked.
Completion requires physical gear contact and firmware disarm; the runner never force-arms or force-disarms.

The rigid-body FDM starts freely at rest on four spring/damper feet. Their reaction forces and moments
enter the equations of motion and accelerometer specific force. There is no airborne hold, imposed trajectory or landing teleport.
Feet are ±0.40 m longitudinally, ±0.35 m laterally and 0.30 m below CG. Each foot has stiffness 20,000 N/m,
vertical damping 650 Ns/m, tangential damping 100 Ns/m and friction coefficient 0.5. These are assumptions, not measured PT1 gear.
At 21 kg, the resting CG is about 0.297 m above the ground datum; the horizon altitude shows CG height.
ArduPilot's default VTOL takeoff altitude is a climb above the takeoff point; the subsequent loiter uses altitude relative to home.

Ground-session parameters added to the existing logged-gain hover baseline:

| Parameter | Value | Purpose |
| --- | ---: | --- |
| Q_GUIDED_MODE | 1 | VTOL loiter in AUTO |
| Q_OPTIONS | 65536 | Disable fixed-wing landing approach (bit 16) |
| Q_NAVALT_MIN | 2 m | Keep early takeoff level before position navigation |
| Q_WP_SPD_UP | 1 m/s | Moderate climb |
| Q_WP_SPD_DN | 1 m/s | Moderate upper descent |
| Q_LAND_FINAL_SPD | 0.5 m/s | Final landing descent |
| TERRAIN_ENABLE | 0 | This scenario uses a flat local datum, with no terrain service |

Existing VTOL controller gains, PT1 mass (21 kg), four KDE lift systems, MAD pusher estimate,
October AeroDB, density 1.021 and moment origin (0,0,0) remain the selected model.
Conventional aero coefficients still have limited measured coverage during vertical flight;
the model has no rotor ground effect, propwash, battery sag, actuator dynamics or validated ground handling.
Forward transition and physical PT1 forward trim require separate validation.

Validation records and plots are in `../artifacts/ground_start_2026_10_11/`.

Sources: [ArduPilot QuadPlane AUTO missions](https://ardupilot.org/plane/docs/quadplane-auto-mode.html),
[pinned 4.7.1 mission handling](https://github.com/ArduPilot/ardupilot/blob/Plane-4.7.1/ArduPlane/commands_logic.cpp),
[pinned 4.7.1 QuadPlane implementation](https://github.com/ArduPilot/ardupilot/blob/Plane-4.7.1/ArduPlane/quadplane.cpp).
