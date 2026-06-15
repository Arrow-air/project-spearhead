# Flight Dynamics Model

Minimal Python 6-DOF flight dynamics and simulation tool for Project Spearhead.

The model uses NED inertial axes and aircraft body axes with `x` forward, `y` right, and `z` down. The default force/moment path uses the latest delivered drag-corrected ADB v1.1 table when it is present in the repository:

```text
docs/information-note/phase-1/Aerodynamics/0002-Preliminary-Design-Aerodynamic-Database/assets/adb_v1_1_drag_correction.csv
```

The ADB CSV metadata and the latest aerodynamics/Nondimit documentation are the source of truth for aerodynamic reference values. If those values disagree with old FDM defaults, the ADB metadata wins.

## ADB Reference Handling

The delivered Nondimit export contains both nondimensional coefficient columns and dimensional force/moment columns. The FDM currently consumes the body-axis nondimensional coefficient columns:

```text
cfx, cfy, cfz, cmx, cmy, cmz
```

It dimensionalizes those coefficients at the current simulation dynamic pressure using the ADB reference geometry from the CSV metadata:

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

## Test

```bash
PYTHONPATH=src/tools/flight_dynamics_model pytest src/tools/flight_dynamics_model/tests
```

The current trim/control derivatives are preliminary. ADB v1.1 does not include control-surface deflection dimensions, so control increments remain provisional and separate from the ADB interpolation.

## Stability Analysis Reproducibility

Preliminary stability outputs and information-note material are tracked, but the
corresponding stability-analysis source module/script is missing from the
current tree. See `docs/stability_reproducibility_gap.md`. Restoring that source
should be a follow-up change, and it should call the shared simulation backend
rather than introducing a separate model path.
