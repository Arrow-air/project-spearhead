#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$TOOL_DIR"

APP_NAME="Nondimit"
VERSION=""
if [[ -f VERSION.txt ]]; then
	VERSION="$(tr -d '[:space:]' < VERSION.txt)"
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-macos/.venv-build}"
DIST_DIR="${DIST_DIR:-macos/dist}"
BUILD_DIR="${BUILD_DIR:-macos/build}"
TARGET_ARCH="${TARGET_ARCH:-}"

if [[ "$(uname -s)" != "Darwin" ]]; then
	echo "This build script must be run on macOS."
	exit 1
fi

"$PYTHON_BIN" -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

python - <<'PY'
import tkinter
print(f"tkinter OK: Tk {tkinter.TkVersion}")
PY

python -m pip install --upgrade pip
python -m pip install --upgrade pyinstaller

python nondimit.py --self-test

rm -rf "$DIST_DIR" "$BUILD_DIR" "$SCRIPT_DIR/${APP_NAME}.spec" "${APP_NAME}.spec"

PYINSTALLER_ARGS=(
	--noconfirm
	--clean
	--windowed
	--name "$APP_NAME"
	--osx-bundle-identifier "com.arrowair.project-spearhead.nondimit"
	--distpath "$DIST_DIR"
	--workpath "$BUILD_DIR"
	--specpath "$SCRIPT_DIR"
	--add-data "$TOOL_DIR/README.md:."
	--add-data "$TOOL_DIR/aero_frame_notation.png:."
)

if [[ -n "$TARGET_ARCH" ]]; then
	PYINSTALLER_ARGS+=(--target-arch "$TARGET_ARCH")
fi

python -m PyInstaller "${PYINSTALLER_ARGS[@]}" nondimit.py

APP_PATH="$DIST_DIR/$APP_NAME.app"
if [[ ! -d "$APP_PATH" ]]; then
	echo "Expected app bundle was not created: $APP_PATH"
	exit 1
fi

"$APP_PATH/Contents/MacOS/$APP_NAME" --self-test

ZIP_STEM="$APP_NAME-macOS"
if [[ -n "$VERSION" ]]; then
	ZIP_STEM="$ZIP_STEM-$VERSION"
fi
if [[ -n "$TARGET_ARCH" ]]; then
	ZIP_STEM="$ZIP_STEM-$TARGET_ARCH"
fi

ditto -c -k --sequesterRsrc --keepParent "$APP_PATH" "$DIST_DIR/$ZIP_STEM.zip"

echo "Build complete:"
echo "$APP_PATH"
echo "$DIST_DIR/$ZIP_STEM.zip"
