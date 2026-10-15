# Setup guide: Flight Deck, FDM and ArduPilot SITL

This guide takes a fresh computer to a running Flight Deck: the Spearhead flight dynamics model
(FDM), ArduPlane 4.7.1 software-in-the-loop (SITL) and the browser interface. Nothing here flashes
flight-controller hardware.

| Platform | Status |
|---|---|
| macOS 13+ (Apple silicon or Intel) | Supported; main development platform |
| Linux (Ubuntu 22.04 / 24.04, Debian 12; other distributions with the same packages) | Supported; verified from a fresh clone on Ubuntu 24.04 |
| Windows 10/11 | Through WSL2 with Ubuntu (a Linux environment inside Windows) |
| Windows without WSL2 | Not supported (ArduPilot SITL does not build natively) |

You need about 3 GB of free disk space (firmware source and build about 1 GB, Python packages,
recordings), Python 3.10 or newer, Git, a C/C++ compiler and internet access for the first run.
The first run takes 10–20 minutes, most of it the firmware build.

## 1. Install the prerequisites

### macOS

1. Install the Xcode command line tools: `xcode-select --install`.
2. Install Python 3.10 or newer, from [python.org](https://www.python.org/downloads/macos/) or
   with Homebrew (`brew install python@3.12`). Check with `python3 --version`.

### Linux

Ubuntu or Debian:

```sh
sudo apt update
sudo apt install -y git python3 python3-venv python3-dev build-essential
```

Fedora: `sudo dnf install -y git python3 python3-devel gcc gcc-c++ make`.

If `python3 --version` reports older than 3.10 (Ubuntu 20.04), install a newer Python and pass it
as `PYTHON=python3.11` in step 3.

### Windows (WSL2)

1. In PowerShell as administrator: `wsl --install -d Ubuntu-24.04`, then restart and create the
   Ubuntu user when asked.
2. Open **Ubuntu** from the Start menu and follow the Linux steps above inside it.
3. Keep the repository in the Ubuntu home folder (`~/…`), not under `/mnt/c/…`; the build is
   much slower and file permissions break on the Windows drive.

## 2. Get the code

```sh
git clone https://github.com/Arrow-air/project-spearhead.git
cd project-spearhead/src/tools/flight_dynamics_model
```

Use the branch or tag named in the release notes if the work is not on `main` yet
(`git checkout <branch>`).

## 3. Set up and start

```sh
sh scripts/run_flight_deck.sh
```

On the first run the script:

1. creates a Python environment in `.sitl/venv` and installs `requirements-gui.txt` and
   `requirements-sitl.txt`;
2. clones ArduPilot at the pinned **Plane-4.7.1** commit `dbe79216` into `.sitl/ardupilot-4.7.1`,
   builds the SITL executable and checks its identity (`.sitl/arduplane-4.7.1`);
3. starts the Flight Deck on the first free port from 8080 and prints the address.

Open the printed address, normally **http://127.0.0.1:8080/**, in a browser. On Windows, open it
in a Windows browser; WSL2 forwards `localhost` to Ubuntu.

Later runs skip the finished steps and start in a few seconds. `Ctrl+C` stops the Flight Deck.
Options:

| Command | Effect |
|---|---|
| `sh scripts/run_flight_deck.sh --setup` | Install and build only |
| `sh scripts/run_flight_deck.sh --check` | Install and build, then run the configuration check and the test suite |
| `PORT=8090 sh scripts/run_flight_deck.sh` | Use a fixed port |
| `PYTHON=python3.12 sh scripts/run_flight_deck.sh` | Use a specific Python for the first setup |
| `SITL_JOBS=2 sh scripts/run_flight_deck.sh` | Fewer parallel build jobs (low-memory machines) |

Everything the script creates stays under `.sitl/`, which Git ignores. Deleting `.sitl/` resets
the setup.

## 4. First flight

In the Flight Deck choose **PT1.5 → Hover · 20 m → Hover hold → 60 s** and press **Start flight**.
Calibration, EKF checks, normal arming and motor start take one to two minutes of simulated time
before the flight clock starts. On a slow computer the simulation runs slower than real time; the
page shows its speed. **Recorded flights** reopens finished runs.

The [Flight Deck guide](flight_deck_ground_control.md) and the
[test tools](flight_deck_test_tools.md) describe the other tests (doublets, transitions, wind,
failures, comparison with an aircraft log).

## 5. Command-line runs (optional)

From `src/tools/flight_dynamics_model` with the environment active (`. .sitl/venv/bin/activate`)
and `export PYTHONPATH=$PWD`:

```sh
# Hover hold with the real firmware, 60 s
python -m spearhead.sitl.quad --binary .sitl/arduplane-4.7.1 --output .sitl/runs/hover-001 --duration 60

# Hover pitch doublet
python -m spearhead.sitl.quad --doublet pitch_doublet --binary .sitl/arduplane-4.7.1 \
  --output .sitl/runs/hover-pitch-001 --duration 45

# Model tests
python -m pytest tests -q
```

Each run writes a new folder with `manifest.json`, `truth.jsonl` (FDM state at 200 Hz),
`telemetry.jsonl` (MAVLink), the firmware's DataFlash log under `logs/` and `report.json`.
More commands are in the [SITL guide](ardupilot_sitl_quickstart.md).

## Troubleshooting

| Symptom | Fix |
|---|---|
| `Python 3.10+ required` | Install a newer Python; run with `PYTHON=python3.12` |
| `error: command 'gcc' failed` or `x86_64-linux-gnu-gcc` missing | Install the compiler packages of step 1 |
| `Source commit differs from pinned Plane-4.7.1; refusing build` | Delete `.sitl/ardupilot-4.7.1` and run again |
| Firmware build killed or very slow | Run again with `SITL_JOBS=1` and close other programs |
| `address already in use` | Another Flight Deck is running; stop it or set `PORT` |
| Page does not open on Windows | Use `http://localhost:<port>/`; check that WSL2 (not WSL1) is in use with `wsl -l -v` |
| A run fails with a wall-clock deadline or timeout | The computer is too busy for the run's time allowance; run one simulation at a time and close other programs |
| Restricted or sandboxed terminal | Local UDP/TCP sockets must be allowed; run in a normal terminal |

## What is and is not modelled

The simulation is a provisional PT1.5 model (October aerodynamic database, estimated inertia and
CG). It is for test planning and software checks, not flight clearance. See the
[validation summary](validation_summary.md) and
[Flight Dynamics note 0003](../../../../docs/information-note/phase-1/Flight-Dynamics/0003-PT1-5-Hover-and-Fixed-Wing-Assessment/information-note.md).
