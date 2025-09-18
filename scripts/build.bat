@echo off
REM ============================================
REM Clean build script for PyInstaller
REM Uses main.spec and deletes old build/dist
REM ============================================

echo.
echo === Cleaning old build and dist folders ===
echo.

IF EXIST build (
    rmdir /s /q build
)

IF EXIST dist (
    rmdir /s /q dist
)

echo.
echo === Building executable with PyInstaller ===
echo.

poetry run pyinstaller main.spec

echo.
echo === Build finished! EXE is in the "dist" folder ===
pause
