# macOS Build Notes

This document is for maintainers or agents preparing a macOS Nondimit package. Standard users should use the provided app package or the source launcher.

## Requirements

- Run this only on macOS. A macOS app bundle cannot be built from the Windows EXE or from a Windows workspace.
- Python 3 must be available as `python3` and must include Tkinter.
- Internet access is needed the first time because the script installs PyInstaller into a local virtual environment.
- Start from the repository root or from this folder. The repository-relative folder is `src/tools/Nondimit/macos`.

## Local Build

From the repository root, run:

```bash
cd src/tools/Nondimit/macos
chmod +x Build_Nondimit_macOS.sh
./Build_Nondimit_macOS.sh
```

If already inside `src/tools/Nondimit/macos`, run only:

```bash
chmod +x Build_Nondimit_macOS.sh
./Build_Nondimit_macOS.sh
```

The script handles the full build sequence:

1. Creates `macos/.venv-build`.
2. Verifies that Tkinter is available.
3. Installs PyInstaller.
4. Runs `python nondimit.py --self-test` from the tool root.
5. Builds `Nondimit.app`.
6. Runs the packaged app self-test.
7. Creates the ZIP package.

## Expected Output

After a successful build, these files should exist relative to `src/tools/Nondimit`:

```text
macos/dist/Nondimit.app
macos/dist/Nondimit-macOS-v1.1.zip
```

The script prints both paths at the end. Treat the build as failed if either file is missing or if any self-test command exits with an error.

## Optional Architecture

By default, the script builds for the native architecture of the Mac or runner. To request a specific PyInstaller target architecture, set `TARGET_ARCH`:

```bash
TARGET_ARCH=universal2 ./Build_Nondimit_macOS.sh
```

Useful values are `universal2`, `arm64`, and `x86_64`. Universal builds require a Python runtime that supports the requested architecture. If that is not available, build the native app on the target Mac architecture instead.

Use `dist/Nondimit-macOS-v1.1.zip` as the macOS package unless a signed or notarized release package is prepared separately.
