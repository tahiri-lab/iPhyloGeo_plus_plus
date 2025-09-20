# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files
import sys

sys.path.append('.')

a = Analysis(
    ['main.py'],
	pathex=['C:\\Users\\agaco\\Documents\Phylogeo\\iPhyloGeo_plus_plus\\.venv\\Lib\\site-packages'],
    binaries=[],
    datas=[],
    hiddenimports=["qtmodern", "qtmodern.styles", "qtmodern.windows",
	"qtpy", "qtpy.QtCore", "qtpy.QtGui", "qtpy.QtWidgets",
	"PyQt6.QtCore", "PyQt6.QtGui", "PyQt6.QtWidgets",
	"PyQt6.QtPrintSupport", "PyQt6.QtSvg", "PyQt6.QtOpenGL", "event_connector", "numpy"],
    hookspath=['hooks'],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['PySide6'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='main',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
