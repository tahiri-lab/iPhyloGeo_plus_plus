@echo off
REM ============================================
REM Clean build script for cx_Freeze
REM Deletes old build/dist
REM ============================================

echo.
echo === Cleaning old build and dist folders ===
echo.

cd "..\scripts"

IF EXIST build (
    rmdir /s /q build
)

IF EXIST dist (
    rmdir /s /q dist
)

echo.
echo === Building executable with cx_Freeze ===
echo.

poetry run python setup.py build

echo.
echo === Build finished! If it succeeded, the EXE is in the "build\exe.win-amd64-3.11" folder ===
pause
