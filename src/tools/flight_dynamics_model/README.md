# Spearhead Flight Deck

Local 6-DOF aircraft simulation with **ArduPilot / ArduPlane 4.7.1** running as a native Software-in-the-Loop (SITL) executable. ArduPilot supplies the motor/servo commands; the Spearhead FDM supplies physics and simulated sensors. Flight Deck provides flight tests, numeric instruments, plots and saved recordings.

## Install and run

macOS, Linux and Windows (through WSL2) are supported; see the [setup guide](docs/setup_guide.md)
for prerequisites per platform. With Git, Python 3.10+ and a C/C++ compiler installed:

```sh
git clone https://github.com/Arrow-air/project-spearhead.git
cd project-spearhead/src/tools/flight_dynamics_model
sh scripts/run_flight_deck.sh
```

The first run creates `.sitl/venv`, installs the Python packages, builds the pinned ArduPlane 4.7.1
SITL into `.sitl/arduplane-4.7.1` (several minutes) and starts the Flight Deck; later runs start in
seconds. Open the printed address, normally **http://127.0.0.1:8080/**. Local binaries,
environments and flight recordings stay under `.sitl/`, outside Git. Nothing is flashed to
flight-controller hardware.

## First flight

Select **PT1.5 → Hover · 20 m → Hover hold → 60 s**, then **Start flight**. Calibration, navigation checks, normal arming and motor startup usually take 1–2 minutes before the flight clock starts. Simulation may run slower than real time; the page shows its progress and speed. **Stop flight** ends the run; **Recorded flights** reopens results.

Hover and forward-trim starts offer pitch/roll/yaw/throttle doublets. Transition tests support both hover → forward and forward → hover. **Interactive flight** allows mode changes with neutral pilot sticks. Advanced settings choose the simulation CG (with its static margin, elevator trim and hover headroom), wind and gusts, and an injected failure. PT1.5 defaults to the target CG x = −40 mm, z = +65 mm from the October origin (a simulation target, not a measurement); rotor and pusher stations follow the 14 October drawing. **Review → Compare with aircraft** overlays a real DataFlash log and lists parameter differences. See [Flight Deck test tools](docs/flight_deck_test_tools.md).

## Verify and troubleshoot

From the same folder and activated environment:

```sh
PYTHONPATH=. python -m spearhead.sitl validate --config sitl_configs/quad_hover.yaml
PYTHONPATH=. python -m pytest tests -q
```

`sh scripts/run_flight_deck.sh --check` runs both. If port 8080 is occupied, the script takes the next free port (or set `PORT`). Local UDP/TCP connections must be allowed. Run in a normal local terminal; restricted sandboxes can block simulator sockets/process startup. Logs and manifests are saved under `.sitl/flight-deck/`.

The integration checks the pinned **Plane-4.7.1 / dbe79216** identity and firmware version. A newer binary requires an explicit compatibility review and fresh flight checks.

## Model scope and details

PT1.5 is a provisional 21 kg model using the October AeroDB/control revision. Measured CG, inertia (currently a build-up estimate), the pusher mount height and motor/propeller, and servo calibration remain open; ±30° travel is an aero-team proposal. Sparse aerodynamic data, vortex-lattice rate derivatives, omitted propwash and a generic battery curve (telemetry only) limit interpretation. Simulation passes do not validate the physical aircraft.

See the [changes relative to main](docs/update_2026_10_14.md), the [setup guide](docs/setup_guide.md), the [validation summary](docs/validation_summary.md), [capabilities and diagrams](docs/capabilities.md), [SITL guide](docs/ardupilot_sitl_quickstart.md), and [model/API reference](docs/model_reference.md). The validation summary records the latest source-only regression and corrected-reference flight checks. Generated recordings, reports and plots remain local; repository formatting/spelling checks still require cleanup.
