# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['test.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=["qtpy", "qtpy.QtCore", "qtpy.QtGui", "qtpy.QtWidgets",
	"PyQt6.QtCore", "PyQt6.QtGui", "PyQt6.QtWidgets",
	"PyQt6.QtPrintSupport", "PyQt6.QtSvg", "PyQt6.QtOpenGL", "PyQt6.QtWebEngineWidgets", "qtmodern", "qtmodern.styles", "qtmodern.windows"], #add Qt too?
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='test',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='test',
)
