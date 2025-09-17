# main.spec
# Generated for a onefile build of main.py including qtmodern and PyQt5

import os

# No code encryption
block_cipher = None

from PyInstaller.utils.hooks import collect_data_files
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT

qtmodern_datas = collect_data_files("qtmodern", includes=["**/*.qss"])
pyqt5_datas = collect_data_files("PyQt5")
qt_folder = os.path.join(os.path.abspath("."), "Qt")

# Analysis: detect dependencies
a = Analysis(
    ["main.py"],
    pathex=[],
    binaries=[],
	datas=qtmodern_datas + pyqt5_datas + [(qt_folder, "Qt")],
    hiddenimports=[
        "qtmodern",
        "qtmodern.styles",
        "qtmodern.windows",
        "PyQt5",       # replace with "PySide2" if needed
        "PyQt5.QtWidgets",
        "PyQt5.QtCore",
        "PyQt5.QtGui"
    ],
    hookspath=["hooks"],  # your custom hooks folder
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)

# Build Python archive
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# Create EXE
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="main",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,  # set False if you want no console
)

# Collect everything into one folder (needed for onefile)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="main",
)
