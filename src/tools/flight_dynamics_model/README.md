# Flight Dynamics Model

Minimal Python 6-DOF flight dynamics and simulation tool for Project Spearhead.

The external ArduPilot JSON backend is pinned to **ArduPlane 4.7.1** and
requires the October tail-update AeroDB. See the
[SITL quickstart](docs/ardupilot_sitl_quickstart.md) for build/run commands
and the [integration plan](docs/ardupilot_sitl_plan.md) for staged acceptance.
The bootstrap smoke test verifies the real firmware connection; controlled
flight and PT1 hover are separate acceptance stages.

The model uses NED inertial axes and aircraft body axes with `x` forward, `y` right, and `z` down. Reference Spearhead retains the original ADB v1.1. PT1 KDE uses the October 2026 tail-update database with partial measured controls. Hardware and aerodynamic versions can be selected independently. The preserved original database is:

```text
docs/information-note/phase-1/Aerodynamics/0002-Preliminary-Design-Aerodynamic-Database/assets/adb_v1_1_drag_correction.csv
```

The selected ADB CSV metadata defines aerodynamic reference values. Explicit user corrections are recorded in a versioned runtime copy, with original source files preserved. For the October delivery, the user confirmed rho=1.021 kg/m³ and moment center=(0,0,0), overriding a contradictory raw CSV moment-center header.

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

- `reference_design.yaml`: planned MAD V8013 PRO IPE 150KV / AMPX 80A V2 /
  FLUXER PRO 26×7.8 lift hardware, with the original preliminary pusher model.
- `pt1_as_built.yaml`: four KDE5215XF-220 / KDE-UAS75HVC+ / 21.5×7.3 two-blade
  lift units, 12S at 44.4 V nominal, pusher absent per the KDE first-hover package.

In Flight Deck, use **Setup → Aircraft** to select **Reference Spearhead** or
**PT1 KDE**. Switching clears the previous run and selects a compatible scenario.
**Setup → Aerodynamics** selects the original or October dataset independently.
Only scenarios with a successful bounded trim are offered for trimmed flight.
The LIFT control bar and bench-estimate thrust/current/power readouts
show the separately commanded lift system.

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

PT1 mass, CG, inertia and rotor coordinates remain preliminary; the profile
does not claim to model the wingless first-hover airframe. The original MAD lift
hardware is preserved, but its lift performance curve is not yet available.
Reference retains ADB v1.1 and its original 0.28 m tail chord. PT1 defaults to
the new database; its numeric updated tail chord/area are pending.
See [profile data and model limits](docs/aircraft_profiles.md) before interpreting
PT1 results or replacing the AeroDB.

The new database includes 32,761 alpha/beta points and 17 partial deflection
cases. Measured control increments replace provisional controls; rate damping
remains provisional. Outside control coverage, endpoints are held and the GUI
reports the approximation. At the existing CG, fixed-wing trim is not achieved
within the measured ±5° elevator range. See the
[comparison and integration report](docs/aerodb_tail_update_2026_10_07.md).

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
python -m pip install -r src/tools/flight_dynamics_model/requirements-gui.txt
```

Run the GUI from the repository root:

```bash
PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.gui --port 8081
```

For startup smoke checks without opening a browser:

```bash
PYTHONPATH=src/tools/flight_dynamics_model \
python -m spearhead.gui --no-show --port 8081
```

Open `http://127.0.0.1:8081/` in a browser. The default route opens the
**Flight Deck**, which is the realtime simulation console.

### Flight Deck

Use the left setup rail to select a scenario and a predefined maneuver. Bundled
scenario YAML files are loaded through the existing scenario adapter. Maneuver
presets are represented with the existing `ControlInputCommand` model.

Available Flight Deck controls:

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
Flight Deck remains responsive during playback.

### Analyze

The **Analyze** page is secondary to the Flight Deck. It can review the latest
Flight Deck run and provides access to the existing batch simulation, stability,
and CG sweep views. CG sweep remains an advanced analysis tool and is not part
of the primary realtime workflow.

### GUI Limitations

- GUI dependencies are optional and are not required for backend simulations or
  stability analysis.
- Realtime playback is timer-based inside NiceGUI, so browser/event-loop timing
  can affect visual update cadence.
- Joystick, gamepad, and manual control input are not implemented yet.
- Wind and disturbance editing are not implemented yet.
- Full 3D aircraft rendering is not implemented yet.
- CG sweep remains an advanced analysis tool rather than a Flight Deck workflow.
- In restricted sandbox environments, NiceGUI startup can fail while creating a
  process pool because system semaphore queries are blocked. Run the GUI in a
  normal local shell when this occurs.

### Manual GUI Checklist

- Launch the GUI.
- Confirm Flight Deck opens by default.
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
- Open Analyze and confirm it sees the latest Flight Deck run.

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

- The tracked aerodynamic database is the nominal ADB only.
- CG sweep changes `AircraftParams.cg_body_xyz` at runtime and relies on the
  existing moment-shift path; it does not recreate the missing Nondimit per-CG
  CSV cases from the historical information note.
- The inertia matrix, propulsion model, rate damping, and control increments
  remain preliminary.
