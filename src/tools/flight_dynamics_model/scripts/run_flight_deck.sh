#!/bin/sh
# Set up (first run) and start the Spearhead Flight Deck on macOS, Linux or WSL2.
#
#   sh scripts/run_flight_deck.sh            # install if needed, build SITL if needed, start
#   sh scripts/run_flight_deck.sh --setup    # install and build only, do not start
#   sh scripts/run_flight_deck.sh --check    # install/build, then run the model tests
#
# Environment: PYTHON (default python3), PORT (default: first free port from 8080),
# SITL_JOBS (parallel firmware build jobs, default 4).
set -eu
fdm_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$fdm_root"
mode=${1:-run}
python=${PYTHON:-python3}
venv="$fdm_root/.sitl/venv"

if ! command -v "$python" >/dev/null 2>&1; then
	echo "Python 3.10 or newer is required; set PYTHON=/path/to/python3" >&2
	exit 1
fi
"$python" -c 'import sys; assert sys.version_info >= (3, 10), "Python 3.10+ required, found " + sys.version.split()[0]'

# Reinstall packages when the requirement files change.
stamp=$(cat requirements.txt requirements-gui.txt requirements-sitl.txt | cksum | cut -d' ' -f1)
if [ ! -x "$venv/bin/python" ]; then
	echo "Creating the Python environment in .sitl/venv"
	mkdir -p .sitl
	"$python" -m venv "$venv"
fi
if [ "$(cat "$venv/.requirements-stamp" 2>/dev/null || true)" != "$stamp" ]; then
	echo "Installing Python packages"
	"$venv/bin/python" -m pip install --quiet --upgrade pip
	"$venv/bin/python" -m pip install --quiet -r requirements-gui.txt -r requirements-sitl.txt
	echo "$stamp" > "$venv/.requirements-stamp"
fi

if [ ! -x .sitl/arduplane-4.7.1 ]; then
	echo "Building ArduPlane 4.7.1 SITL (first run only; several minutes)"
	SITL_PYTHON="$venv/bin/python" sh scripts/build_ardupilot_sitl.sh
fi

export PYTHONPATH="$fdm_root"
export OPENBLAS_NUM_THREADS=1
case "$mode" in
	--setup)
		echo "Setup complete. Start with: sh scripts/run_flight_deck.sh"
		exit 0 ;;
	--check)
		"$venv/bin/python" -m spearhead.sitl validate --config sitl_configs/quad_hover.yaml
		exec "$venv/bin/python" -m pytest tests -q ;;
	run) ;;
	*)
		echo "Unknown option: $mode (use --setup or --check)" >&2
		exit 2 ;;
esac

port=${PORT:-$("$venv/bin/python" -c 'import socket
for p in range(8080, 8100):
	s = socket.socket()
	try:
		s.bind(("127.0.0.1", p)); print(p); break
	except OSError:
		pass
	finally:
		s.close()')}
echo "Flight Deck: http://127.0.0.1:$port/  (Ctrl+C to stop)"
exec "$venv/bin/python" -m spearhead.gui --host 127.0.0.1 --port "$port" --no-show
