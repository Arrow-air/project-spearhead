# Flight Deck guide

The Flight Deck is the browser interface to the ArduPilot simulation. It plans a test, starts
the firmware and the FDM, shows the flight live and keeps the recording. Install and start it
with the [setup guide](setup_guide.md); how the firmware connection works is in
[SITL integration](sitl_integration.md).

The Flight Deck flies **PT1.5** (`pt15`, 21 kg) only. The 25 kg reference design stays available
in Model Lab and on the command line. Results are simulations of a provisional model, not
flight clearance.

## Layout

| Area | Content |
|---|---|
| **Plan** | Aircraft, starting point, test, duration; collapsible starting conditions (CG, wind, failure) |
| **Fly** | Attitude view, flight values, live alerts, progress and simulation speed |
| **Review** | Overview, Messages, Flight plots, Technical details, Compare with aircraft; Recorded flights |
| **Analyze** | Model Lab (model only, no firmware) and the CG window |

The page is at `/` (`/sitl` and `/embed` show the same page). One server runs one flight at a
time; closing the tab does not stop it. **Stop flight** ends the firmware and FDM processes.

## Starting points and tests

| Starting point | Tests |
|---|---|
| Ground · VTOL takeoff | Takeoff → hover → land |
| Hover · 20 m, Custom hover | Hover hold, Interactive flight, hover doublets |
| Hover · 80 m (transition) | Hover → forward flight test, plus the hover tests |
| Forward trim · 22 / 24 / 26 m/s, Custom forward trim (22–26 m/s) | Straight & level hold, Interactive flight, Forward → hover transition, forward doublets |

- **Hover doublets:** pitch or roll attitude ±8°, yaw stick ±20 %, lift collective ±10 points;
  2 s per half after 10 s of settling; at least 45 s. Pitch, roll and yaw run in QHOVER;
  collective uses QSTABILIZE for the input and returns to QHOVER. The input is checked from the
  firmware's own RC readback.
- **Forward doublets:** elevator ±3° for 0.75 s per half, aileron ±5° for 0.5 s, rudder ±5° for
  0.75 s, pusher up to ±10 points for 2 s (limited symmetrically by the trim headroom). Input in
  MANUAL, recovery in FBWA; at least 45 s.
- **Hover → forward:** 10 s hover, FBWA request, then a 5 s pusher ramp that starts at confirmed
  FBWA; at least 60 s. Losing FBWA returns the pusher to neutral.
- **Forward → hover:** settle 10 s at the forward trim, then QLOITER braking and position hold;
  at least 60 s.
- **Takeoff → hover → land:** an AUTO mission (VTOL takeoff, VTOL-mode hold, timed loiter, VTOL
  land) on spring/damper landing gear. Choose takeoff height (3–30 m) and hover time (15–600 s);
  climb and landing add to it. The run ends with gear contact and firmware disarm.
- **Interactive flight:** after startup, request any firmware mode from the selector and watch
  what the firmware does. Pilot sticks stay neutral (hover throttle in QSTABILIZE/QACRO, center
  throttle in altitude-holding Q modes); there is no joystick. At most 10 simulated minutes.

All controller gains are frozen: the VTOL gains logged on PT1 and the 9 October fixed-wing
proposal. Pass criteria for each test are in [SITL integration](sitl_integration.md#runners).
A pass covers the stated checks only. Peaks, drift, altitude offsets and use of extrapolated
aerodynamic data stay visible in the result even when the checks pass.

## Starting conditions

**CG.** The simulation CG (target, not measured) defaults to x −40 mm, z +65 mm from the CFD
origin. Alternatives: PT1 logged x −47.5 mm, aft −60 mm, forward −25 mm and 0 mm. Below the
choice the panel shows, for that CG, the static margin at 22 m/s, elevator trim at 20/22/24 m/s
with ±30° travel (marked when beyond the measured ±15°) and the busiest lift rotor's thrust as a
fraction of its 95 % spin thrust. Forward trim is recomputed for the chosen CG; gains do not
change.

**Wind and gusts.** Steady wind 0–15 m/s from a chosen direction (meteorological "from") and an
optional 1-cosine gust up to 10 m/s in the same direction, starting 20 s after release and
lasting 3 s. Aerodynamics, pusher inflow and the simulated pitot use the air-relative velocity.
The field is uniform: no turbulence, shear or rotor wake. Forward starts are trimmed in still
air and meet the wind at release.

**Failure injection.** One failure at a chosen time after release:

| Failure | Effect | Available for |
|---|---|---|
| Lift motor FR / RR / RL / FL out | That rotor's command is zeroed | Hover, ground, forward → hover |
| Pusher cut | Pusher command is zeroed | Forward starts, hover → forward |
| Low battery | Lift pack drops to 15 % charge | Disabled until the SITL battery parameters are approved |

The firmware is not told; it sees only the aircraft's response. Recordings keep both the
commanded and the applied control.

## While flying

Startup (calibration, EKF, arming, motor start) takes about 1–2 minutes before the flight
clock starts. The page shows the stage, the simulator clock and the measured simulation speed;
below 1× the computer is slower than real time, which does not change the physics.

**Flight values** beside the attitude view: true airspeed, altitude, climb rate, heading, body
rates, elevator angle, pusher command and current, lift-pack voltage, current and charge used,
wind, and the four lift PWM values (FL/FR above RL/RR). Missing values show a dash; stale and
recorded values are labelled.

**Live alerts:** elevator beyond the measured ±15° data (or any other control beyond measured
coverage), forward flight between 3 and 20 m/s where the aerodynamic data are least reliable,
pusher above its 24 A sustainable current, lift pack below 3.5 V per cell, and an injected
failure.

**Status row:** link, armed state and mode from the firmware heartbeat, GPS fix and satellites,
EKF navigation validity, sensor-health failures and the reported battery.

## Interactive mode commands

The selector lists the QuadPlane modes of the pinned build: MANUAL, CIRCLE, STABILIZE,
TRAINING, ACRO, FBWA, FBWB, CRUISE, AUTOTUNE, AUTO, RTL, LOITER, TAKEOFF, AVOID_ADSB, GUIDED,
QSTABILIZE, QHOVER, QLOITER, QLAND, QRTL, QAUTOTUNE, QACRO, THERMAL and LOITER_ALT_QLAND, each
with a short description. Commands show queued, sent, accepted and confirmed; acceptance alone
does not count as a mode change, only the heartbeat does. The firmware may refuse a mode or fall
back to another one.

Some modes need inputs the Flight Deck does not supply: AUTO needs a mission (only the
ground-start runner uploads one), GUIDED needs targets, THERMAL needs a soaring setup. Commands
are disabled during startup, on a stale link, in scripted tests and in recordings.

## Review

- **Messages:** the latest 200 firmware and Flight Deck messages with severity, timed from
  release (T−/T+).
- **Flight plots:** roll, pitch, body rates, heading and sideslip, alpha (hidden below 1 m/s),
  speed and altitude, firmware outputs or lift PWM with the 1100–1940 µs limits, elevator against
  its trim, surface travel against the assumed limits, pusher command and, for doublets, the
  requested and measured input. Time is flight time: 0 at release, startup negative.
- **Technical details:** firmware and database identity, recorded CG, control coverage, forward
  trim, parameter overlay.
- **Recorded flights:** reopens finished runs and command-line runs under `.sitl/`, read-only.

### Compare with aircraft

Takes a DataFlash `.BIN` from the aircraft (a path on this computer or an upload) and compares it
with the open SITL run:

- lift outputs, roll and pitch, altitude and pack voltage, from takeoff (log) and release or
  spin-up (SITL), at 10 Hz; a log with several flights lists them, choose one by number
  (0 = last);
- median lift PWM per rotor after 5 s and the front share of PWM above idle;
- parameter differences between the aircraft and the run's SITL parameters (read from the run's
  own DataFlash log), leaving out calibration, telemetry, hardware and simulator settings.

A `.param` file gives the parameter comparison only.

## Analyze

- **CG window:** static margin, elevator trim and hover headroom over a list of CG positions,
  the same quantities as the CG summary. October-AeroDB estimates, not a measured envelope.
- **Model Lab** (`/model-lab`): the model without firmware. Choose aircraft, aerodynamic dataset,
  scenario and maneuver; play, pause, step and reset at 0.5–5× speed with attitude, telemetry,
  control and response plots. Batch simulation and stability analysis use the same model.

## Limits

- No joystick or manual stick flying; no general parameter editor or mission editor.
- One local session per server; no hosted multi-user mode.
- Rotor inflow, propwash, ground effect, actuator dynamics and battery sag on thrust are not
  modelled; landing gear properties are assumptions.
- Control surfaces have no aerodynamic effect at zero airspeed (no propwash model). In QHOVER and
  QSTABILIZE the firmware still moves aileron and elevator; the rudder is centered.
- QHOVER holds attitude and height, not position: hover runs drift horizontally.
