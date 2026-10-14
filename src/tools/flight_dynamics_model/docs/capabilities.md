# Flight Deck capabilities and scope

**Status: provisional engineering/research tool · 12 October 2026.**

Flight Deck runs the actual **ArduPlane 4.7.1** autopilot as a native Software-in-the-Loop (SITL) executable, connected to the Spearhead six-degree-of-freedom aircraft model. Use it to examine trims, controller response, VTOL transitions and model assumptions. Start with the [installation README](../README.md).

## What controls the aircraft?

ArduPilot owns navigation estimation, normal arming, flight modes, controllers, motor mixing and servo outputs. The FDM converts those outputs to forces/moments, integrates aircraft motion, and supplies simulated sensor data. Flight Deck observes the same run; opening a browser does not create another physics model.

![ArduPilot and FDM exchange motor/servo commands and simulated sensor feedback.](diagrams/firmware_model.svg)

The firmware source commit is pinned to `dbe79216`; its executable identity and MAVLink version are checked. Each run records its binary hash, source/configuration hashes, selected CG and parameter overlay. This is computer-based simulation; no flight controller is flashed.

## Available experiments

| Capability | Current behavior |
|---|---|
| PT1.5 hover starts | 20 m, 80 m or custom altitude; normal arming and checked motor startup before release |
| PT1.5 ground mission | Automatic takeoff, timed hover and landing with provisional flat-ground gear contact |
| Forward starts | Computed six-axis trim; PT1.5 named speeds 22/24/26 m/s and custom 22–26 m/s; reference profile also offers 28 m/s |
| Hover → forward | 10 s hover, FBWA request, then a 5 s pusher ramp anchored to confirmed mode; minimum 60 s |
| Forward → hover | Starts at selected forward trim, settles 10 s, requests QLOITER and checks sustained low-speed hover; minimum 60 s |
| Hover doublets | Pitch/roll attitude ±8°, yaw stick ±20%, lift collective ±10 percentage points; 2 s per lobe; minimum 45 s |
| Forward doublets | Elevator ±3° for 0.75 s, aileron ±5° for 0.5 s, rudder ±5° for 0.75 s, pusher up to ±10 points for 2 s per lobe; minimum 45 s |
| Interactive flight | Firmware-confirmed mode requests; neutral pilot sticks, not joystick/manual-stick flying |
| Reference fixed-wing checks | Holds, pitch/roll steps, small attitude disturbances and physical surface-direction checks |
| Recorded review | Reopen completed managed runs and local CLI diagnostics with read-only controls |

Doublets settle for 10 s, deliver positive and negative pulses, restore neutral input and observe recovery with frozen gains. Forward inputs use MANUAL during excitation and FBWA during recovery; hover attitude/yaw use QHOVER. Lift collective briefly uses QSTABILIZE, because QHOVER throttle commands climb rate instead of direct collective. Pusher doublets are limited symmetrically by available trim headroom; the recorded amplitude may be smaller at a high-throttle trim.

QHOVER holds altitude/attitude and allows horizontal drift; QLOITER adds position holding. Ailerons and elevator can still move in hover under ArduPilot; the rudder is centered in QHOVER/QSTABILIZE. Lift motors provide the dominant hover control in the saved tests. The FDM has no rotor-propwash model, so aerodynamic surface authority vanishes at exactly zero airspeed.

## Follow and interpret a run

![Four stages from flight planning to recorded review, with early-stop handling.](diagrams/flight_lifecycle.svg)

Startup occurs before the free-flight clock. Progress distinguishes wall time from simulated time; a measured speed below 1× means the computer is taking longer than real time. The server owns one shared flight, and closing a tab does not stop it. The 3D implementation is retained but disabled; hidden charts do not render/update their browser plots.

Beside the horizon, numeric values show true airspeed, altitude, earth-vertical climb rate, heading, body p/q/r rates, pusher actuator command and four lift PWM outputs arranged FL/FR above RL/RR. Values use the FDM's physical state and actual firmware outputs. Unknown/no-signal values remain dashes; saved and stale samples are labeled. The motor plot shows the logged 1100–1940 µs limits. Messages expose arming, mode, GPS/EKF and other firmware events.

**A passing result has a defined scope.** Transition acceptance requires a sustained interval with fresh recorded mode evidence, suitable motion and the expected motor state. Firmware “Transition done” alone is insufficient. Doublet input delivery is checked separately from response/recovery. Interactive completion does not establish stability. Peaks, altitude offsets, drift and aerodynamic extrapolation remain visible even when input checks pass.

## Aircraft model, trim and CG studies

The shared model provides nonlinear rigid-body dynamics, measured database interpolation/control increments, provisional rate damping, four independent lift rotors, a pusher estimate and optional ground loads. The PT1.5 profile is 21 kg with KDE lift hardware; the reference profile is 25 kg. Hardware and aerodynamic profiles can be selected independently in the model tools; the current SITL integration requires the October dataset.

![Aerodynamic loads, one CG moment transfer, trim and local analysis.](diagrams/cg_trim.svg)

Body axes are x forward, y right, z down; inertial axes are NED. The runtime moment transfer is `M_CG = M_source + (r_source − r_CG) × F`. The active October base is first numerically restored from its documented Nondimit export center (-0.1,0,+0.15) to the CFD origin (0,0,0). Changing CG then recomputes its trim and moment arm; it does not require regenerating or double-shifting the database. The provisional PT1.5 default is [+0.010,0,+0.150] m, 10 mm forward of the origin. All prototype starts use the selected CG, including hover and ground; z remains unmeasured. Hover starts use checked collective or explicitly selected independent-rotor trim. After release, the firmware owns all actual actuator outputs.

The model/API tools also provide open-loop scenarios, batch simulation/export, longitudinal trim, numerical linearization, longitudinal/lateral modes, damping/growth metrics, static derivative estimates and CG sweeps. These analyses are local to the selected trim and inherit the model's assumptions. A balanced trim does not prove stability; an unreachable trim is reported as a failure. See the [model/API reference](model_reference.md).

## Evidence and current boundaries

The [final validation summary](validation_summary.md) records the tested configurations, key results and remaining limits.

Earlier forward-transition, reverse-transition and doublet recordings used the former header-only moment-reference handling and do not validate the corrected model. Current trims and flight checks are in the [validation summary](validation_summary.md). At the forward +10 mm candidate, 22/24/26 m/s trim and local longitudinal modes are stable in this model; trim uses elevator extrapolation and 26 m/s leaves little pusher headroom. These calculations do not establish a physical CG envelope.

Recordings, numeric reports, plots and run manifests remain local under `.sitl/` or ignored `artifacts/` folders. They are not distributed with the repository; each new run generates its own evidence.

Physical CG/datum, inertia, final motor geometry, pusher mounting/thrust line and servo calibration remain unconfirmed. ±30° travel is a proposal, not measured linkage calibration. Control data are sparse and predominantly β=0; extrapolation/endpoint assumptions are explicit. Base/control neutral pitching-moment disagreement remains unresolved. Rotor inflow/propwash/transients, battery sag, detailed wind and broader physical-aircraft validation are absent or unvalidated. Reference MAD lift performance is not yet available for a validated reference hover.

The current package is suitable as a **provisional research baseline on GitHub** with these limits stated. It is not a hardware flight recommendation, complete flight-envelope validation, or hosted multi-user simulator. macOS is tested; a Linux setup recipe is supplied but not verified end to end. Firmware upgrades require explicit compatibility checks and new maneuver verification.
