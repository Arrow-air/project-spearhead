# Spearhead flight dynamics model and Flight Deck

A six-degree-of-freedom model of the Spearhead quadplane, run either on its own (open-loop
simulation, trim, stability) or with the real **ArduPlane 4.7.1** firmware as software in the
loop (SITL). ArduPilot supplies the motor and servo commands; the model supplies the physics and
simulated sensors. The browser-based Flight Deck plans tests, shows the flight and keeps the
recordings. Nothing is flashed to flight-controller hardware.

## Quick start

macOS, Linux and Windows (WSL2) are supported; prerequisites per platform are in the
[setup guide](docs/setup_guide.md). With Git, Python 3.10+ and a C/C++ compiler installed:

```sh
git clone https://github.com/Arrow-air/project-spearhead.git
cd project-spearhead/src/tools/flight_dynamics_model
sh scripts/run_flight_deck.sh
```

The first run creates `.sitl/venv`, installs the Python packages, builds the pinned firmware
into `.sitl/arduplane-4.7.1` (several minutes) and starts the Flight Deck; later runs start in
seconds. Open the printed address, normally **http://127.0.0.1:8080/**, and fly **PT1.5 →
Hover · 20 m → Hover hold → 60 s**. `sh scripts/run_flight_deck.sh --check` runs the
configuration check and the test suite. Everything local stays under `.sitl/`, which Git
ignores.

## Documentation

| Document | Content |
|---|---|
| [Setup guide](docs/setup_guide.md) | Installation on macOS, Linux and Windows; options; troubleshooting |
| [Flight Deck guide](docs/flight_deck.md) | Tests, starting conditions (CG, wind, failures), live values and alerts, review, comparison with aircraft logs, Model Lab |
| [SITL integration](docs/sitl_integration.md) | Firmware pinning, protocol and timing, channels, sensors, startup, command-line runners and pass criteria, recorded files |
| [Model reference](docs/model_reference.md) | Axes, aerodynamic data and moment reference, aircraft profiles, propulsion, CG and stability analysis, Python and CLI use |
| [Validation summary](docs/validation_summary.md) | What has been checked, results, and what is not validated |
| [Changes relative to main](docs/update_2026_10_14.md) | What this work adds to the July model, with migration notes |

Engineering results for PT1 and PT1.5 are in the Flight Dynamics information notes under
[`docs/information-note/phase-1/Flight-Dynamics/`](../../../docs/information-note/phase-1/Flight-Dynamics/).

## Model scope

PT1.5 is a provisional 21 kg model on the October aerodynamic database. The simulation CG
(x −40 mm, z +65 mm) is a target, not a measurement; inertia and rate damping are estimates;
rotor and pusher stations follow the 14 October drawing; ±30° surface travel is a proposal.
Sparse control data, the low-speed aerodynamic data, missing propwash and the generic battery
curve limit interpretation. Simulation passes support test planning; they do not validate the
physical aircraft.
