# Flight Dynamics Model

Minimal Python 6-DOF flight dynamics and simulation tool for Project Spearhead.

The model uses NED inertial axes and body axes with `x` forward, `y` right, and `z` down. The default force/moment path uses the PR #7 drag-corrected ADB v1.1 table when it is present in the repository:

```text
docs/information-note/phase-1/Aerodynamics/0002-Preliminary-Design-Aerodynamic-Database/assets/adb_v1_1_drag_correction.csv
```

If that file is missing, the model falls back to the older placeholder derivative model so scripts remain runnable.

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

## Test

```bash
PYTHONPATH=src/tools/flight_dynamics_model pytest src/tools/flight_dynamics_model/tests
```

The current trim/control derivatives are preliminary. ADB v1.1 does not include control-surface deflection dimensions, so control increments remain provisional and separate from the ADB interpolation.

