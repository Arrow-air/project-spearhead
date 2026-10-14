# Flight Deck test tools (14 October 2026)

These tools prepare and check PT1.5 flight tests. The Flight Deck flies PT1.5 (`pt15`) only;
the 25 kg reference design remains available in Model Lab with its older AeroDB and geometry.

## Wind and gusts

**Starting conditions & model settings → Wind and failures** sets a steady wind (0–15 m/s,
meteorological "from" direction) and an optional 1-cosine gust (0–10 m/s, same direction,
starting 20 s after release and lasting 3 s).

- Aerodynamics, the pusher inflow and the simulated pitot use the air-relative velocity;
  the state keeps the inertial velocity, so drift and groundspeed are real.
- Wind starts at release (ground starts: when the lift motors spin up after arming), so the
  firmware's stationary pitot calibration still sees still air.
- The field is uniform: no turbulence spectrum, wind shear or rotor wake. Forward starts are
  trimmed in still air and meet the wind at release.

Runs record the wind vector, true airspeed, alpha and beta per frame, and the manifest records
the wind settings.

## Failure injection

The same panel injects one failure at a chosen time after release:

| Failure | Effect in the FDM | Available for |
|---|---|---|
| Lift motor FR / RR / RL / FL out | That rotor's command is zeroed | Hover, ground, forward → hover |
| Pusher cut | Forward throttle is zeroed | Forward starts, hover → forward transition |
| Low battery | Lift pack drops to 15 % charge | Disabled until the SITL battery parameters are approved |

The firmware is not told about the failure; it sees only the aircraft's response. Records keep
both the commanded and the applied control.

## Live alerts and readouts

The Fly panel shows alerts while flying:

- elevator beyond the measured ±15° CFD data (or any other measured-coverage limit);
- pusher above its 24 A sustainable current;
- lift pack below 3.5 V per cell;
- failure injected.

Flight values add elevator angle, pusher current, lift-pack voltage and current, lift-pack
charge used and wind speed. Airspeed is air-relative.

The lift pack is reported to ArduPilot through the JSON `battery` field (generic LiPo curve,
45 mΩ, lift motors plus 1.5 A avionics). ArduPilot reads it only with `BATT_MONITOR = 4`,
which awaits approval for the repository SITL parameters. Thrust still uses nominal voltage.

## CG summary and CG window

Below the CG choice the Plan panel shows, for that CG: static margin at 22 m/s, elevator trim
at 20/22/24 m/s with ±30° travel (marked when beyond the measured ±15°) and the busiest lift
rotor's thrust as a fraction of its 95 % spin thrust. **Analyze → CG window** tabulates the same
quantities over a list of CG positions. Both are October-AeroDB estimates, not a measured
envelope.

## Compare with aircraft

**Review → Compare with aircraft** takes a DataFlash `.BIN` from the aircraft (path on this
computer or upload) and compares it with the open SITL run:

- Lift outputs, roll/pitch and altitude/pack voltage, both from takeoff (log) and release or
  spin-up (SITL), at 10 Hz. A log with several flights lists them; choose one by number
  (0 = last).
- Median lift PWM per rotor after 5 s and the front share of PWM above idle.
- Parameter differences between the aircraft and the run's SITL parameters (read from the run's
  own SITL DataFlash log), leaving out calibration, telemetry, hardware and simulator settings.

A `.param` file gives the parameter comparison only.

## Removed or simplified

- One frozen transition parameter set remains: the 9 October proposal. The later "PT1.5"
  comparison overlay was removed; the aircraft's own parameters replace it after alignment.
- The disabled 3D view was deleted. It remains in Git history.
- **Analyze → CG Sweep** became the PT1.5 CG window above.
- "Diagnostic CG" now reads "Simulation CG (target, not measured)".
- Recordings made with the former rotor geometry are labelled "legacy geometry" in Recorded
  flights.
