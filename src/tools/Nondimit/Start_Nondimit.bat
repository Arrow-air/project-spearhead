@echo off
setlocal
cd /d "%~dp0"

if exist "Nondimit.exe" (
    start "" "Nondimit.exe"
) else (
    python "nondimit.py"
)
