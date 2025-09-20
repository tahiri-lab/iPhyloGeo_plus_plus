import os
from PyInstaller.utils.hooks import collect_submodules, collect_all

block_cipher = None

# --- Qt setup ---
# Ensure qtpy picks up the right API at runtime
os.environ["QT_API"] = "pyqt6"

# Collect everything from PyQt6
pyqt6_datas, pyqt6_binaries, pyqt6_hiddenimports = collect_all("PyQt6")

# Collect everything from qtmodern
qtmodern_datas, qtmodern_binaries, qtmodern_hiddenimports = collect_all("qtmodern")

# Collect everything from qtpy
qtpy_datas, qtpy_binaries, qtpy_hiddenimports = collect_all("qtpy")

a = Analysis(
    ["main.py"],
    pathex=["C:\\Users\\agaco\\Documents\Phylogeo\\iPhyloGeo_plus_plus\\.venv\\Lib\\site-packages"],
    binaries=pyqt6_binaries + qtmodern_binaries + qtpy_binaries,
    datas=pyqt6_datas + qtmodern_datas + qtpy_datas,
    hiddenimports=pyqt6_hiddenimports
        + qtmodern_hiddenimports
        + qtpy_hiddenimports
        + [
            "qtmodern",
            "qtmodern.styles",
            "qtmodern.windows",
            "qtpy",
            "qtpy.QtCore",
            "qtpy.QtGui",
            "qtpy.QtWidgets",
            "PyQt6",
            "PyQt6.QtCore",
            "PyQt6.QtGui",
            "PyQt6.QtWidgets",
			"event_connector"
        ],
    hookspath=["hooks"],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["PySide6"],  # ensure no mixed Qt bindings
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name="main",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
)
