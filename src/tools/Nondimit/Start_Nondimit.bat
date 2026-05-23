@echo off
setlocal

cd /d "%~dp0"

if exist "%~dp0Nondimit.exe" (
    start "" "%~dp0Nondimit.exe"
    exit /b 0
)

python "%~dp0nondimit.py"

endlocal
