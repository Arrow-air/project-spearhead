# Spearhead Flight Deck

Local 6-DOF aircraft simulation with **ArduPilot / ArduPlane 4.7.1** running as a native Software-in-the-Loop (SITL) executable. ArduPilot supplies the motor/servo commands; the Spearhead FDM supplies physics and simulated sensors. Flight Deck provides flight tests, numeric instruments, plots and saved recordings.

## Install and run

Tested on macOS with Python 3.10. You need Git, Python, and native C/C++ build tools: Xcode Command Line Tools on macOS; `build-essential`, `python3-dev` and `python3-venv` on Ubuntu. The Linux build recipe is provided but has not been verified end to end here. Internet access is required for dependencies and the first firmware build.

```sh
git clone --branch codex/ardupilot-sitl https://github.com/Arrow-air/project-spearhead.git
cd project-spearhead/src/tools/flight_dynamics_model

python3 -m venv .sitl/venv
. .sitl/venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-gui.txt -r requirements-sitl.txt

SITL_PYTHON="$PWD/.sitl/venv/bin/python" sh scripts/build_ardupilot_sitl.sh

export PYTHONPATH="$PWD"
OPENBLAS_NUM_THREADS=1 python -m spearhead.gui --host 127.0.0.1 --port 8080 --no-show
```

Open **http://127.0.0.1:8080/**. The first build downloads the pinned ArduPilot source/submodules and can take several minutes. It creates `.sitl/arduplane-4.7.1`; local binaries, environments and flight recordings stay outside Git. It does not flash flight-controller hardware.

## First flight

Select **PT1.5 → Hover · 20 m → Hover hold → 60 s**, then **Start flight**. Calibration, navigation checks, normal arming and motor startup usually take 1–2 minutes before the flight clock starts. Simulation may run slower than real time; the page shows its progress and speed. **Stop flight** ends the run; **Recorded flights** reopens results.

Hover and forward-trim starts offer pitch/roll/yaw/throttle doublets. Transition tests support both hover → forward and forward → hover. **Interactive flight** allows mode changes with neutral pilot sticks. Advanced settings expose diagnostic CG and frozen forward-parameter choices. PT1.5 defaults to the target CG x = −40 mm, z = +65 mm from the October origin (a simulation target, not a measurement); rotor and pusher stations follow the 14 October drawing. The retained 3D view is disabled.

## Verify and troubleshoot

From the same folder and activated environment:

```sh
PYTHONPATH=. python -m spearhead.sitl validate --config sitl_configs/quad_hover.yaml
PYTHONPATH=. python -m pytest tests -q
```

If port 8080 is occupied, start the GUI with `--port 8081`. Local UDP/TCP connections must be allowed. Run in a normal local terminal; restricted sandboxes can block simulator sockets/process startup. Logs and manifests are saved under `.sitl/flight-deck/`.

The integration checks the pinned **Plane-4.7.1 / dbe79216** identity and firmware version. A newer binary requires an explicit compatibility review and fresh flight checks.

## Model scope and details

PT1.5 is a provisional 21 kg model using the October AeroDB/control revision. Measured CG, inertia (currently a build-up estimate), the pusher mount height and motor/propeller, and servo calibration remain open; ±30° travel is an aero-team proposal. Sparse aerodynamic data, vortex-lattice rate derivatives, omitted propwash and a generic battery curve (telemetry only) limit interpretation. Simulation passes do not validate the physical aircraft.

See the [validation summary](docs/validation_summary.md), [capabilities and diagrams](docs/capabilities.md), [SITL guide](docs/ardupilot_sitl_quickstart.md), and [model/API reference](docs/model_reference.md). The validation summary records the latest source-only regression and corrected-reference flight checks. Generated recordings, reports and plots remain local; repository formatting/spelling checks still require cleanup.
