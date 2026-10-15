#!/bin/sh
# Native macOS/Linux SITL only. Does not flash a board or modify another checkout.
set -eu
fdm_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
sitl_source="$fdm_root/.sitl/ardupilot-4.7.1"
sitl_python="${SITL_PYTHON:-python3}"
expected_commit=dbe792162d06cab66c3475fd5556bf7a120f119e
mkdir -p "$fdm_root/.sitl"
if [ ! -d "$sitl_source" ]; then
	git clone --depth 1 --branch Plane-4.7.1 https://github.com/ArduPilot/ardupilot.git "$sitl_source"
fi
actual_commit=$(git -C "$sitl_source" rev-parse HEAD)
if [ "$actual_commit" != "$expected_commit" ]; then
	echo "Source commit differs from pinned Plane-4.7.1; refusing build" >&2
	exit 1
fi
if [ -n "$(git -C "$sitl_source" status --porcelain --untracked-files=no)" ]; then
	echo "Source checkout has modifications; refusing build" >&2
	exit 1
fi
git -C "$sitl_source" submodule update --init --recursive \
	modules/waf modules/mavlink modules/DroneCAN/DSDL \
	modules/DroneCAN/pydronecan modules/DroneCAN/dronecan_dsdlc \
	modules/DroneCAN/libcanard modules/littlefs modules/gtest modules/lwip
cd "$sitl_source"
"$sitl_python" ./waf configure --board sitl --no-submodule-update
"$sitl_python" ./waf plane -j "${SITL_JOBS:-4}"
PYTHONPATH="$fdm_root" "$sitl_python" -c \
	'from spearhead.sitl.backend import verify_firmware; print(verify_firmware("build/sitl/bin/arduplane"))'
cp build/sitl/bin/arduplane "$fdm_root/.sitl/arduplane-4.7.1"
