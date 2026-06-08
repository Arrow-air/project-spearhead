#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$TOOL_DIR"

if [[ -d "$SCRIPT_DIR/Nondimit.app" ]]; then
	open "$SCRIPT_DIR/Nondimit.app"
elif [[ -d "$SCRIPT_DIR/dist/Nondimit.app" ]]; then
	open "$SCRIPT_DIR/dist/Nondimit.app"
elif command -v python3 >/dev/null 2>&1; then
	python3 "nondimit.py"
else
	osascript -e 'display dialog "Python 3 is required to run Nondimit from source." buttons {"OK"} with icon stop'
	exit 1
fi
