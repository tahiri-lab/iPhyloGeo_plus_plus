@echo off
REM ============================================
REM Build script for PyInstaller using main.spec
REM ============================================

echo.
echo === Building executable with PyInstaller ===
echo.

REM Activate Poetry environment and run PyInstaller
poetry run pyinstaller main.spec

echo.
echo === Build finished! EXE is in the "dist" folder ===
pause