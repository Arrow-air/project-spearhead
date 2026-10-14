# Model reference and development notes

Dated validation reports and recordings referenced below are retained locally and excluded from Git.

For installation and the current Flight Deck, start with the [setup README](../README.md).
These notes retain the detailed model/API documentation and dated validation history.

**12 October reference correction:** earlier recorded results below used header-only origin handling. The active base now numerically restores the Nondimit export from (-0.1,0,+0.15) to (0,0,0), then transfers moments once to the selected CG. Those earlier results do not validate this corrected model. See the [current validation summary](validation_summary.md).


Python 6-DOF flight dynamics, stability analysis and ArduPilot simulation for Project Spearhead.

The prototype is displayed as **PT1.5**; its compatible configuration ID remains `pt1_as_built`.
The **Hover → forward flight test** starts at 80 m and uses the prototype’s calculated 24 m/s servo trim.
The 11 October `deflection-3.csv` delivery is now the active October control revision.
Its integration and flight mechanics review (local report)
records the trim checks and the frozen earlier transition proposal, without controller tuning.
The CG sensitivity campaign (local report)
compares explicit CG variants with the same model and older controller proposal.
The CG/reference audit (local report)
verifies runtime moment transfer against native Nondimit exports and distinguishes
the October source origin from the inherited provisional aircraft CG.
Historical forward-CG cases passed with the former reference handling. The corrected model defaults to x=+10 mm forward; z=+150 mm remains unmeasured. Physical CG and the flight envelope remain unconfirmed.
See the PT1.5 transition report (local report) for both trials, parameter changes and plots.
The pitch diagnosis (local report)
replays the recorded handover, separates pitch-moment contributions and quantifies
the current model's unstable mode and limited elevator recovery authority.

The external ArduPilot JSON backend is pinned to **ArduPlane 4.7.1** and
requires the October tail-update AeroDB. See the
[SITL quickstart](ardupilot_sitl_quickstart.md) for build/run commands
and the [integration plan](ardupilot_sitl_plan.md) for staged acceptance.
The bootstrap smoke test verifies the real firmware connection. Normal arming,
MANUAL direction checks, nominal FBWA, small steps and combined disturbances
are exercised by `python -m spearhead.sitl.flight`; the
fixed-wing control report (local report)
records results and model limits. The main Flight Deck at `/` starts/stops
owned runner processes, displays the FDM horizon alongside numeric flight values,
and opens
completed recordings. `/sitl` is a compatible link to the same console. See the
web integration and PWM replay report (local report).
The interactive session adds real firmware mode commands, link/GPS/EKF/sensor
status, and a timestamped warning/message feed. See the
[ground-control guide](flight_deck_ground_control.md) and
verification results (local report).
The **Elevator ±1° doublet** preset applies a physical MANUAL elevator input
through ArduPlane at 32 m/s; its surface-angle trace shows the actual FDM input.
See the input and trim assessment (local report).
New PT1.5 hover, ground and transition flights use the aero team's proposed **±30° surface travel** over 1000–2000 µs (11 October update). The reference fixed-wing presets retain **±20°**, matching their earlier checks.
with linear continuation of measured control slopes outside the local measured coverage, bounded
at the assumed travel. Alpha/beta endpoint assumptions remain. Travel values
are configurable and marked uncalibrated in the UI and run manifest. The proposal does not establish measured linkage or aerodynamic coverage, and older ±20° transition passes do not validate ±30°. CG awaits assembly; inertia is unavailable; the latest motor-coordinate drawing is pending. The physical pusher is not yet installed and its thrust-line placement is undecided; CG-aligned thrust remains a simulation assumption. See the
travel assessment (local report).
Pitch/roll step and disturbance acceptance with this wider travel is tracked
in the control revalidation report (local report),
including revised provisional firmware gains and the unsuccessful tuning trials.
The isolated-axis follow-up (local report)
checks both signs of roll and pitch disturbances at 28 m/s with those same
gains and records local mass/CG sensitivity. Higher-speed firmware tests are
deferred at the user's request. The earlier malformed `deflection-2.csv` remains
inactive; the corrected `deflection-3.csv` is integrated separately with provenance.
The trim sweep (local report)
examines speed/mass sensitivity at the same diagnostic CG. Broader
envelope sensitivity, real-time performance, repeated-run determinism and
long-duration hover remain separate acceptance checks. PT1.5 ground start now
supports **takeoff → timed hover → landing**, with normal arming and physical
gear contact. The [ground-start guide](ground_start.md) describes the
provisional gear and mission; the validation report (local report)
records successful 5, 10 and 20 m tests with ArduPlane 4.7.1.

The model uses NED inertial axes and aircraft body axes with `x` forward, `y` right, and `z` down. Reference Spearhead retains the original ADB v1.1. PT1.5 KDE uses the October 2026 tail-update database with partial measured controls. Hardware and aerodynamic versions can be selected independently. The preserved original database is:

```text
docs/information-note/phase-1/Aerodynamics/0002-Preliminary-Design-Aerodynamic-Database/assets/adb_v1_1_drag_correction.csv
```

The selected ADB CSV metadata defines aerodynamic reference values. Explicit user corrections are recorded in a versioned runtime copy, with original source files preserved. For the October delivery, the user confirmed raw CFD rho=1.021 kg/m³ and source moment center=(0,0,0). The preserved Nondimit export header records its final processed moment center (-0.1,0,+0.15), rather than the raw CFD origin. The active table numerically restores its moments to (0,0,0); force corrections already present in the delivery are retained without applying them again. `scripts/recenter_october_database.py` reproduces this conversion from the preserved source, with independent wrench and hash checks.

## ADB Reference Handling

The delivered Nondimit export contains both nondimensional coefficient columns and dimensional force/moment columns. The FDM currently consumes the body-axis nondimensional coefficient columns:

```text
cfx, cfy, cfz, cmx, cmy, cmz
```

It dimensionalizes those coefficients at the current simulation dynamic pressure using the selected ADB reference geometry. The preserved v1.1 reference is:

```text
rho = 1.09 kg/m^3
V_ref = 25 m/s
S = 1.805 m^2
b = 3.95 m
c = 0.457 m
moment_center_body_xyz = (-0.15, 0.0, 0.15) m
```

The CSV also includes wind-axis coefficient columns (`cd`, `cs`, `cl`, `cr`, `cm`, `cn`) and dimensional body/wind columns (`Fx`, `Fy`, `Fz`, `Mx`, `My`, `Mz`, `D`, `S`, `L`, `R`, `M`, `N`) for review and regeneration checks. The runtime FDM path uses body-axis coefficients because the equations of motion expect body-axis loads.

## CG Handling

`AircraftParams.cg_body_xyz` is the target center of gravity in the same body-axis coordinate system as the ADB moment reference. The default CG is the documented nominal/reference CG:

```text
(-0.15, 0.0, 0.15) m
```

For each aerodynamic lookup, the model converts the ADB coefficients into dimensional body-axis force and moment about the database reference center, then shifts moments to `cg_body_xyz`:

```text
M_CG = M_ref + cross(reference.moment_center_body_xyz - cg_body_xyz, F_body)
```

This runtime transfer is the correct path for CG-dependent stability studies using one delivered database. Nondimit regeneration can still be used for manual reference checks or official re-exporting at a new moment center, but the runtime FDM must remain consistent with the delivered DB metadata.

If the recommended ADB file is missing, the model falls back to the older placeholder derivative model so scripts remain runnable.

## Aircraft Profiles

One FDM supports two bundled configurations in `aircraft_configs/`:

- `reference_design.yaml`: 25 kg planned MAD V8013 PRO IPE 150KV / AMPX 80A V2 /
  FLUXER PRO 26×7.8 lift hardware, with the original preliminary pusher model.
- `pt1_as_built.yaml`: 21 kg PT1.5 with four KDE5215XF-220 / KDE-UAS75HVC+ /
  21.5×7.3 lift units and a MAD V8013 PRO IPE 150KV / AMPX 80A / APC 20×12
  pusher estimate, using 12S at 44.4 V nominal. See the
  pusher estimate (local report) for sources and limits.

The Flight Deck offers independent **Aircraft**, **Initial point**, **Maneuver**
and duration choices. PT1.5 supports air-start hover or ground takeoff/hover/landing;
Spearhead supports diagnostic forward trims. New ArduPilot sessions use the
October AeroDB. The model-only **Model Lab** under **Analyze** also provides
independent aircraft/aerodynamic selections and actuator/thrust readouts.

The same profiles are available through the API and CLI:

```python
from spearhead import load_aircraft_config, SimulationConfig, run_simulation

params = load_aircraft_config("reference_design")
result = run_simulation(SimulationConfig(params=params))
```

```bash
PYTHONPATH=src/tools/flight_dynamics_model python -m spearhead.sim run \
  src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml --aircraft reference_design

PYTHONPATH=src/tools/flight_dynamics_model python -m spearhead.sim run \
  src/tools/flight_dynamics_model/scenarios/pt1_hover.yaml --aircraft pt1_as_built
```

PT1.5 mass is configured as 21 kg; CG, inertia and rotor coordinates remain preliminary. The profile
does not claim to model the wingless first-hover airframe. The original MAD lift
hardware is preserved, but its lift performance curve is not yet available.
Reference retains ADB v1.1 and its original 0.28 m tail chord. PT1.5 defaults to
the new database; its numeric updated tail chord/area are pending.
See [profile data and model limits](aircraft_profiles.md) before interpreting
PT1.5 results or replacing the AeroDB.

The base database retains 32,761 alpha/beta points. The active 11 October control
delivery contains 103 cases: 97 beta-zero cases, including four elevator/rudder
combinations, and six neutral sideslip cases preserved for audit. Sparse control
samples extend to ±15° elevator/rudder and ±30° aileron at selected alpha values;
elevator also has ±5° samples at 10 and 15 m/s. This is not rectangular coverage.
Measured increments replace provisional control terms; rate damping remains
provisional. Alpha/speed endpoint holds, deflection extrapolation and interaction
assumptions are reported explicitly. The previous 49-case revision is selectable
as `tail_update_2026_10_07_controls_2026_10_09`.
See the current review (local report)
and [reference handling](aerodb_tail_update_2026_10_07.md).

```bash
PYTHONPATH=src/tools/flight_dynamics_model python -m spearhead.sim run \
  src/tools/flight_dynamics_model/scenarios/aero_database_check.yaml \
  --aircraft reference_design --aerodynamics tail_update_2026_10_07
```

API: `load_aircraft_config("pt1_as_built", aerodynamics="adb_v1_1")` restores
old aerodynamics with KDE hardware. Scenario YAML accepts `aerodynamics`; both
simulation and stability CLIs accept `--aerodynamics`. Explicit selections win.

## Install Dependencies

From the repository root:

```bash
python -m pip install -r src/tools/flight_dynamics_model/requirements.txt
```

## Run

From the repository root:

```bash
python src/tools/flight_dynamics_model/scripts/run_open_loop.py
python src/tools/flight_dynamics_model/scripts/run_trim.py
```

The shared simulation backend is also available as a typed API:

```python
from spearhead import SimulationConfig, run_simulation

config = SimulationConfig(duration=20.0, dt=0.01, n_points=None)
result = run_simulation(config)
```

Future GUIs, notebooks, batch scripts, and report generators should call this
same backend API. They should not implement their own flight dynamics model.

Scenario YAML files can be run through the thin CLI wrapper:

```bash
PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.sim run src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml

PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.sim run src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml \
  --csv output.csv --json output.json
```

For future real-time frontends, use `RealtimeSimulation`:

```python
from spearhead import RealtimeSimulation, SimulationConfig

sim = RealtimeSimulation(SimulationConfig(duration=20.0, dt=0.01, n_points=None))
while sim.running:
    state = sim.step()
```

## Optional GUI

The Spearhead FDM includes an optional NiceGUI/Plotly frontend. The GUI is a
frontend layer over the shared backend APIs; it does not duplicate trim,
dynamics, aerodynamics, propulsion, stability, CG sweep, or export logic.

Install the core FDM dependencies first, then install the optional GUI
dependencies:

```bash
python -m pip install -r src/tools/flight_dynamics_model/requirements.txt
python -m pip install -r src/tools/flight_dynamics_model/requirements-gui.txt \
  -r src/tools/flight_dynamics_model/requirements-sitl.txt
```

Run the GUI from the repository root:

```bash
PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.gui --port 8082
```

For startup smoke checks without opening a browser:

```bash
PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.gui --no-show --port 8082
```

Open `http://127.0.0.1:8082/` in a browser. The default route opens the
ArduPilot-backed **Flight Deck**. Live firmware sessions require the pinned
SITL executable; follow the [build quickstart](ardupilot_sitl_quickstart.md).
The Python model, web server and firmware run locally; the browser renders the
instruments. The retained 3D view is disabled. The default server owns one shared session.
Docker packaging and hosted independent sessions remain future work.

### Flight Deck

Select **Aircraft**, **Initial point**, **Maneuver** and simulation duration,
then press **Start flight**. PT1.5 ground starts select **Takeoff → hover → land**;
the time field becomes **Hover time**, and climb/landing add to that dwell.
Normal calibration and arming happen before the flight sequence. **Stop flight**
terminates the local simulation processes. **Recorded flights** opens completed runs.

For the transition diagnostic select **PT1.5 → Hover · 80 m (transition) →
Hover → forward flight test → 60 s**. Choose **Transition CG** explicitly:
default forward x=+10 mm, comparison x=+25/+50 mm, or the historical aft x=−150 mm; y=0 and z=+150 mm
in body axes (x forward, z down). The source moment center stays at (0,0,0).
The default **Transition parameters** choice is the frozen 9 October proposal;
the later PT1.5 simulation overlay is available for comparison. Forward trim
is recomputed for the chosen CG without changing gains. **Technical details**
shows the actual recorded CG, control coverage and forward trim; the elevator
plot marks that forward trim. Transition results distinguish firmware handover
from sustained-flight acceptance. The recording picker also includes completed
CLI diagnostics under the local `.sitl` directory, with read-only controls.
See the CG comparison (local report)
and Flight Deck review (local report).

The attitude horizon shows the shared FDM's physical state. The 3D aircraft view is temporarily disabled: its import, construction and updates are commented out, while its implementation remains available for restoration.
The page follows **Plan → Fly → Review**. Starting conditions and frozen parameter choices are collapsible; firmware/database identifiers remain in Technical details. Startup shows observed stages, elapsed wall time and a rough wait estimate. Flight duration uses simulated time; measured simulation speed is displayed separately. UI updates at 10 Hz, with charts refreshed only while the plots tab is visible.

The five-second forward-throttle ramp begins at confirmed FBWA, including delayed confirmations; losing that mode returns throttle to neutral and requires fresh confirmation. Forward acceptance now requires fresh recorded FBWA heartbeats (maximum age 2 simulated seconds) throughout the accepted interval, plus the existing speed, lift-idle, attitude, altitude and duration checks. A passing interval does not imply that the entire maneuver stayed within those limits. The result summary separately shows peak pitch magnitude, altitude loss and the first 10 s hover peak. Plots include alpha (hidden below 1 m/s), pitch rate, pusher command and surface travel with assumed limits.

For an explicit forward-CG diagnostic, first create a local configuration from the bundled template. From the repository root, after installing/building SITL:

```sh
cp src/tools/flight_dynamics_model/sitl_configs/quad_transition.yaml \
  src/tools/flight_dynamics_model/.sitl/forward-cg.yaml
```

Set `cg_body_xyz_m: [0.01, 0.0, 0.15]` in that local copy. This selects a diagnostic CG, not a measured aircraft CG. Then run:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src/tools/flight_dynamics_model \
  python -m spearhead.sitl.quad --transition \
  --transition-parameters first_transition_2026_10_09 \
  --config src/tools/flight_dynamics_model/.sitl/forward-cg.yaml \
  --binary src/tools/flight_dynamics_model/.sitl/arduplane-4.7.1 \
  --output src/tools/flight_dynamics_model/.sitl/forward-cg-check \
  --duration 60 --wall-timeout 900 --port 19712 --mavlink-port 16470
```

Use a new output directory and free ports for each run. Source data, physical parameters and controller gains are unchanged by the ramp/acceptance corrections. See the repository preparation results (local report).

Typography uses locally installed Neue Haas Grotesk with Arial/system fallbacks; proprietary font bytes are not bundled. JetBrains Mono v2.304 is distributed with its OFL notice in `spearhead/gui/assets/fonts/`.

Interactive flights enable firmware-confirmed mode changes after initialization.
Automatic test sequences and recordings keep those commands disabled.
GPS/EKF/system status and timestamped firmware warnings remain visible.
The **Flight plots** and **Technical details** tabs expose traces and provenance.
See the [ground-control guide](flight_deck_ground_control.md).

### Model Lab

Open **Analyze → Model Lab** (`/model-lab`) for model-only realtime experiments.
Use the left setup rail to select a scenario and a predefined maneuver. Bundled
scenario YAML files are loaded through the existing scenario adapter. Maneuver
presets are represented with the existing `ControlInputCommand` model.

Available Model Lab controls:

- **Play/Pause** starts or pauses realtime stepping through `RealtimeSimulation`.
- **Reset** returns the selected scenario and maneuver to `t = 0`.
- **Step** advances one simulation step while paused.
- **Speed** selects the playback multiplier (`0.5x`, `1x`, `2x`, `5x`).
- The transport strip shows elapsed time, total scenario time, status, and a
  timeline progress bar.

The primary flight display is a PFD-style artificial horizon. The aircraft
symbol remains fixed while the sky/ground world and pitch ladder move with the
simulated attitude.

Telemetry cards show large-value readouts for:

- IAS: body-speed magnitude in `m/s`
- ALT: altitude derived from NED down position in `m`
- PITCH, ROLL, HDG: Euler angles in degrees
- ALPHA: angle of attack in degrees
- BETA: sideslip angle in degrees

Control input bars show the resolved current control command:

- ELEV: elevator deflection in degrees, centered at zero
- AIL: aileron deflection in degrees, centered at zero
- RUD: rudder deflection in degrees, centered at zero
- THR: throttle command as percent, from zero to full scale

The lower response area contains live Plotly cards for attitude, rates, control
inputs, and airspeed/altitude. These plots use a rolling display history so the
Model Lab remains responsive during playback.

### Analyze

The **Analyze** page is secondary to the Flight Deck. It can review the latest
Model Lab run and provides access to the existing batch simulation, stability,
and CG sweep views. CG sweep remains an advanced analysis tool and is not part
of the primary realtime workflow.

### GUI Limitations

- GUI dependencies are optional and are not required for backend simulations or
  stability analysis.
- Model Lab playback is timer-based inside NiceGUI. The ArduPilot runner advances
  physics in lockstep from sensor timestamps; browser timing only affects display updates.
- Joystick, gamepad, and manual control input are not implemented yet.
- Wind and disturbance editing are not implemented yet.
- The retained, disabled 3D implementation uses a schematic aircraft with
  follow/side/orbit cameras.
- Ground contact uses provisional flat-ground gear. Forward transition, real
  hardware tune, battery sag, rotor ground effect and propwash remain unvalidated.
- CG sweep remains an advanced analysis tool rather than a Flight Deck workflow.
- In restricted sandbox environments, NiceGUI startup can fail while creating a
  process pool because system semaphore queries are blocked. Run the GUI in a
  normal local shell when this occurs.

### Manual GUI Checklist

- Launch the GUI.
- Confirm Flight Deck opens by default.
- Select PT1.5, Ground, Takeoff → hover → land, 10 m and a 15 s hover dwell.
- Press Start flight; confirm ground initialization, takeoff, hover and landing.
- Confirm the completed run shows AUTO, disarmed and a validation result.
- Open Recordings and confirm its controls are read-only.
- Select an air-start Interactive flight and confirm mode controls unlock after calibration.
- Confirm the FDM horizon, external view, health/messages and flight plots update.
- Open Analyze → Model Lab for the model-only controls below.
- Select the pitch doublet maneuver.
- Press Play and confirm simulation time advances.
- Confirm Pause works.
- Confirm Reset returns the run to `t = 0`.
- Confirm Step advances while paused.
- Confirm the speed multiplier changes playback rate.
- Confirm the PFD attitude display updates.
- Confirm telemetry cards update.
- Confirm control input bars update.
- Confirm all four plot cards update.
- Select a roll/aileron maneuver.
- Select a rudder maneuver.
- Open Analyze and confirm it sees the latest Model Lab run.

## Test

```bash
PYTHONPATH=src/tools/flight_dynamics_model pytest src/tools/flight_dynamics_model/tests
```

The current stability model is preliminary. ADB v1.1 uses provisional controls; the October dataset uses partial measured control increments. Rate derivatives remain provisional for both. Failed force/moment balance is rejected even when the optimizer terminates successfully.

## Stability Analysis Reproducibility

The restored stability pipeline lives under `spearhead.stability` and consumes
the same backend as simulations:

```python
from spearhead.scenarios import load_scenario
from spearhead.stability import analyze_stability

config = load_scenario("src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml")
result = analyze_stability(config)
```

CLI examples:

```bash
PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.stability analyze \
  src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml \
  --json stability.json --csv stability_modes.csv --markdown stability.md

PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.stability cg-sweep \
  src/tools/flight_dynamics_model/scenarios/pitch_doublet.yaml \
  --cg-x -0.15 -0.10 -0.05 \
  --json cg_sweep.json --csv cg_sweep_modes.csv
```

The analysis trims through `run_simulation(config)`, then linearizes the existing
`dynamics()` function with central finite differences. Static derivatives are
estimated from the existing force/moment breakdown path. No separate aircraft
model is used.

Current limitations:

- The original ADB and October tail-update dataset are versioned independently;
  active ArduPilot sessions use the October dataset.
- CG sweep changes `AircraftParams.cg_body_xyz` at runtime and relies on the
  existing moment-shift path; it does not recreate the missing Nondimit per-CG
  CSV cases from the historical information note.
- The inertia matrix, propulsion model, rate damping, and control increments
  remain preliminary.

### PT1.5 forward-trim starts

PT1.5 now offers **Forward trim · 22 / 24 / 26 m/s** and **Custom forward trim** (22–26 m/s), alongside Ground, Hover · 20 m, Hover · 80 m and Custom hover. The Hover · 30 m shortcut is removed; old recordings remain readable as Custom hover at their original altitude. Reference-aircraft trim points are unchanged.

Forward starts compute full six-axis trim for the selected speed and diagnostic CG, retain ±30° proposed travel and the MAD pusher PWM envelope, and expose Interactive flight or Straight & level hold. Lift motors remain off and fixed-wing modes are available. Forward parameters use the selected saved baseline, with the frozen 9 October proposal as default; only the trim PWM values are derived from the computed physical trim. CG selection is explicit. The present pusher estimate does not yield a balanced PT1.5 trim at 28 m/s, so that shortcut is not offered.

Forward-start verification (local report).

### Forward-to-hover and analysis doublets

The PT1.5 FW-trim menu includes **Forward → hover transition** (10 s FBWA settlement followed by firmware-managed QLOITER braking/position hold) and **Pitch / Roll / Yaw / Throttle doublets**. They command physical elevator ±3° for 0.75 s per lobe, aileron ±5° for 0.5 s, rudder ±5° for 0.75 s, or up to ±10 throttle points for 2 s. Throttle amplitude is limited symmetrically by trim headroom and recorded. Excitation is in MANUAL, followed by FBWA recovery; frozen controller gains are retained. Run durations are at least 60 s for reverse transition and 45 s for doublets.

Interactive flight means neutral pilot sticks with user-requested mode changes, rather than a scripted maneuver. Input delivery is checked separately from recovery near trim. Body-rate, yaw/sideslip and requested/actual input plots support reviewing each experiment. Nonzero sideslip remains outside measured β=0 coverage. Maneuver verification (local report).
