#hooks/hook-qtpy.py
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = collect_submodules("qtpy")