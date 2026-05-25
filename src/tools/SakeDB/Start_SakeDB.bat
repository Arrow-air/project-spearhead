@echo off
setlocal

cd /d "%~dp0"

if exist "%~dp0SakeDB.exe" (
    start "" "%~dp0SakeDB.exe" %*
    exit /b 0
)

echo SakeDB.exe was not found in %CD%
exit /b 1

endlocal
