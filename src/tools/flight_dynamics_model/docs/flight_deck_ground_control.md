# Flight Deck ground-control functions

The main **Flight Deck** at `/` runs the ArduPilot-backed simulation; `/sitl` remains a compatible link to the same page. The model-only console is **Model Lab**, under **Analyze**, at `/model-lab`. **Recordings** opens the history picker. The attitude and external aircraft views stay together above the **Messages**, **Flight plots**, and **Model & run details** panels. Messages default to warnings and errors; choose **All messages** for the full recent feed. Recording status and controls are explicitly read-only. Model assumptions and run provenance remain available in the details panel.

Open `/`, select **PT1**, **Hover · 20 m**, **Interactive flight** and a duration, then press **Start flight**. The 21 kg PT1 uses four KDE lift motors, an estimated MAD 150KV / APC 20×12 pusher, October AeroDB and official ArduPlane 4.7.1. The air-start runner performs stationary estimator initialization, normal arming, motor spool-up, airborne release and QHOVER entry. It verifies HOME at simulated ground (500 m MSL). Mode controls become available only after the QHOVER heartbeat confirms this bootstrap. Choose a **Requested mode**, read its description, and press **Change mode**. The current mode is taken from the firmware heartbeat.

For a ground start, select **Ground · VTOL takeoff** and **Takeoff → hover → land**. Set takeoff altitude and hover dwell separately. This automatic sequence uses free ground physics, a normally armed AUTO mission, timed VTOL loiter and physical landing contact followed by firmware disarm. Mode commands are disabled while the sequence runs. Gear geometry/compliance/friction are assumptions; see [ground-start setup and validation](ground_start.md).

## Modes

The QuadPlane selector includes all 24 user-selectable modes in our pinned build: MANUAL, CIRCLE, STABILIZE, TRAINING, ACRO, FBWA, FBWB, CRUISE, AUTOTUNE, AUTO, RTL, LOITER, TAKEOFF, AVOID_ADSB, GUIDED, QSTABILIZE, QHOVER, QLOITER, QLAND, QRTL, QAUTOTUNE, QACRO, THERMAL and LOITER_ALT_QLAND. Each choice has a brief description and prerequisite/coverage notes. QAUTOTUNE is compiled into this local 4.7.1 binary, although it is omitted from many downloadable builds and deprecated in current documentation.

INITIALISING (16) is displayed if reported during startup but is never a command choice. AUTOLAND (26) is available in the fixed-wing selector only: the pinned firmware explicitly rejects it when QuadPlane is enabled. Mode 9 is unused. Spearhead **Interactive flight** uses Q_ENABLE=0; its selector includes all 17 user-selectable fixed-wing modes and excludes Q modes. Aircraft configuration determines the allowed requests at both the manager and runner, independent of the browser selection.

After bootstrap, interactive runs do not force the aircraft back to QHOVER/FBWA. Firmware fallbacks and automatic mode changes remain authoritative. Direct-throttle QSTABILIZE and QACRO receive the calculated hover-throttle pilot input, accounting for RC3 deadzone (currently 1505 µs at 21 kg). Altitude-controlled Q modes receive centre throttle (neutral climb request). Roll, pitch and yaw sticks remain neutral. These are pilot inputs, not replacement motor outputs: firmware drives the four independent PWM outputs into the same physical FDM used by the horizon and external 3D view.

## Model and command limits

Mode selection is integrated; it does not validate every mode's flight performance. PT1's pusher uses an analytical motor/propeller estimate; transition and sustained forward flight remain unvalidated, and PT1's physical CG has no validated forward trim. The reference profile retains its planned propulsion model. The dedicated ground-start sequence includes provisional gear contact; air-start interactive flights still stop before the ground. Battery sag, rotor transients, ground effect and propwash are not modelled. The interactive QuadPlane envelope is roll/pitch within 60 degrees and altitude between 5 and 120 m. Scripted hover tests retain tighter limits and QHOVER acceptance criteria.

AUTO requires a suitable uploaded mission. The ground-start runner uploads its built-in takeoff/hover/landing mission; a general mission editor/uploader is not implemented. Air-start interactive runs do not upload a mission. GUIDED target entry and ADS-B traffic are not supplied. THERMAL requires an active soaring configuration; no thermal updraft model is supplied. AUTOTUNE needs pilot manoeuvres; QAUTOTUNE excites the aircraft and can change gains within that isolated run. Requesting a mode does not set its auxiliary parameters or invent the missing mission/targets. Firmware can reject a request or select a fallback; the deck reports that outcome.

Command progress shows queued, sent, accepted and confirmed states; denial and timeout are visible. Acceptance alone never claims that the mode changed. A second request waits for the first to finish. Commands are disabled during startup, stale connections, recordings and scripted tests. The browser targets only its server-owned local simulation, not physical aircraft.

Arming is automatic through the established bootstrap; **Stop run** ends the simulation process group rather than sending a disarm command. Interactive sessions last at most ten simulated minutes and retain estimator-validity checks. Normal completion is **COMPLETED**, not a passed flight-validation test. Logs include manifest/configuration provenance, truth, firmware telemetry, command outcomes and a report.

## Status and messages

The status row shows the link, actual armed state and mode, GPS fix/satellite count, EKF navigation validity, enabled/present sensor-health failures, and simulated battery voltage/remaining charge when reported. Missing or stale values remain unknown. Recorded flights label their last status as recorded and never enable controls.

**Messages & warnings** retains the latest 200 ArduPilot text messages and command results, with timestamps and severity badges. Select **Warnings and errors** to filter the feed. MAVLink 2 multi-chunk messages are assembled before display; missing chunks are marked incomplete when detected. Historical warnings do not by themselves mean a fault is currently active; current sensor/EKF health is shown separately. Runner failures and aerodynamic-coverage warnings remain visible beside the physical-state display.

The single native runner owns the firmware connection and physics clock. Browser requests are atomic files in that run's `commands/` directory. The runner sends MAV_CMD_DO_SET_MODE over its existing connection and records ACK/rejection and heartbeat confirmation in `control-events.jsonl` and `events.json`. Requests expire after ten wall-clock seconds; confirmation times out after eight seconds. The passive browser observer never advances physics or creates a second firmware connection.

PT1 **Hover hold** remains a scripted air-start test with selectable altitude/duration. It uses normal arming and logged VTOL gains, then QSTABILIZE → QHOVER. Mode commands are disabled. QHOVER holds altitude and attitude, not horizontal position, so horizontal drift is expected. Current 21 kg static hover demand is about 67% of the motor PWM range; the firmware thrust setting is separately calculated through motor linearization. Historical recordings retain their original mass and hover demand.

The FDM horizon and external 3D view continue to display physical truth. Estimator status and firmware messages remain separate diagnostics.

## Sources

- [Pinned 4.7.1 mode enum](https://github.com/ArduPilot/ardupilot/blob/dbe792162d06cab66c3475fd5556bf7a120f119e/ArduPlane/mode.h)
- [Pinned AUTOLAND entry restrictions](https://github.com/ArduPilot/ardupilot/blob/dbe792162d06cab66c3475fd5556bf7a120f119e/ArduPlane/mode_autoland.cpp)
- [QuadPlane modes](https://ardupilot.org/plane/docs/quadplane-flight-modes.html)
- [Plane modes](https://ardupilot.org/plane/docs/flight-modes.html)
- [Flight-mode commands and heartbeat confirmation](https://ardupilot.org/dev/docs/mavlink-get-set-flightmode.html)
- [MAVLink STATUSTEXT](https://mavlink.io/en/messages/common.html#STATUSTEXT)
