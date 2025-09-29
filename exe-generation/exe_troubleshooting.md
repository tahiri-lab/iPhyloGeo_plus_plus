# EXE Troubleshooting Notes

## Troubleshooting by Arcady

I am working on a Windows 11 computer and using the Powershell console.

## Troubleshooting steps

### September 18, 2025

1. Generated an EXE of the project with Pyinstaller
2. Got **ModuleNotFoundError** on this line in main.py: `import qtmodern.styles`
3. Created a hooks directory within scripts
4. Created hooks\hook-qtmodern.py (file contents in code snippets below)
5. Added qtmodern and PyQt5 to the Poetry venv (PyQt5 was a mistake as the project actually uses PyQt6; I fixed that later)
6. Updated main.spec to collect data files for qtmodern and PyQt5
7. Specified "qtmoderns" as an hidden import in main.spec
8. Removed and added qtmodern to the Poetry venv
10. Asked [a question](https://stackoverflow.com/questions/79767848/exe-created-with-pyinstaller-fails-to-import-qtmodern-as-specified-in-hook) on StackOverflow to get outside input
11. Got the virtual environment’s path and stored it in a variable called $venvpath in Powershell
12. Generated an EXE using main.py rather than main.spec (this makes Pyinstaller overwrite the main.spec file), with command `poetry run pyinstaller main.py --onefile --paths "$venvPath\Lib\site-packages"  --hidden-import=qtmodern  --hidden-import=qtmodern.styles --hidden-import=qtmodern.windows`
13. Got **ModuleNotFoundError** three lines further than before in main.py: `from event_connector import QtEvents, connect_decorated_methods, connect_event`
14. Since event_connector.py was located in scripts with main.py, added an empty __init__.py file to scripts
15. Edited main.spec to add event_connector to hiddenimports and "." to pathex
16. Pyinstaller **would no longer generate an EXE** file at all. Replaced "." in main.spec’s pathex with the full path on my machine
17. Ran `poetry remove PyQt5` and `poetry add PyQt6`
18. There were **traces of PyQt5 remaining**
19. Deleted the venv and recreated it with command `Invoke-Expression (poetry env activate)` (Powershell-specific syntax)
20. Now able to generate an EXE and run it, got the **ModuleNotFoundError** on the qtmodern import again
21. Removed the active folder from main.spec’s pathex after an online search showed me that Pyinstaller doesn’t handle multiple elements in the pathex list well
22. Ran `poetry remove qtmodern` and `poetry add qtmodern`
23. Pyinstaller gave this **error message** while trying to generate an EXE: `Aborting build process due to attempt to collect multiple Qt bindings packages: attempting to run hook for 'PySide6', while hook for 'PyQt6' has already been run!`
24. Added `excludes['PySide6'],` to  main.spec
25. Now able to generate an EXE and run it, got this error: `ImportError: Error importing numpy: you should not try to import numpy from its source directory; please exit the numpy source tree, and relaunch your python interpreter from there.`
26. Added `collect_data_files("numpy")` to datas in main.spec
27. Got a **ValueError** when trying to generate the EXE: `too many values to unpack (expected 2)`
28. Tried to create a function that would flatten the input and put it into two lists and pass that to datas
29. Emptied datas
30. Learned from an online search that there is a known issue when using Numpy 1.27 with 3.12 in a frozen environment
31. Downloaded, installed and added Python 3.11 to the venv
32. When running the EXE generated with Python 3.11, I get a **PyQt error**: `qtpy.QtBindingsNotFoundError: No Qt bindings could be found` from the `import qtmodern.styles` line in main.py
33. Tried to include everything that seemed relevant in hiddenimports (see code snippets below)
34. Created `hooks/hook-qtpy.py` (file contents in code snippets below)
35. Added the line `os.environ["QT_API"] = "pyqt6"` to main.py

### September 20, 2025

1. Provided ChatGPT with the troubleshooting steps from September 18 and asked it for ideas
2. Ran `poetry run python -c "import qtpy; import PyQt6; print(qtpy.API_NAME)"`, confirming PyQt6 is set up properly inside the venv
3. Replaced the contents of hook-qtpy with those suggested by ChatGPT (see code snippets below)
4. Replaced the contents of hook.qtmodern with those suggested by ChatGPT (see code snippets below)
5. Added "qtpy", "qtpy", "qtpy.QtCore", "qtpy.QtGui" and "qtpy.QtWidgets", to main.spec’s hiddenimports
6. Added the lines `from PyInstaller.utils.hooks import collect_dynamic_libs` and `binaries = collect_dynamic_libs("PyQt6")` to main.spec, and added the binaries variable to the binaries list
7. Got a **ValueError** when trying to generate the EXE: `not enough values to unpack (expected 2, got 0)`
8. Renamed main.spec to main.spec.old to replace it by a new version yet be able to consult the old version (see new version in code snippets below)
9. Added the venv path to the new main.spec and event_connector to its hiddenimports
10. An EXE does get generated; encountering the **PyQt error** again when running it
11. Using Powershell in scripts\dist, ran `pyi-archive_viewer main.exe` but **pyi-archive_viewer** wasn’t recognized as a command, despite Pyinstaller being recognized in the same folder `The term 'pyi-archive-_viewer' is not recognized as a name of a cmdlet, function, script file, or executable program.`
12. Searched my AppData folder for pyi-archive_viewer and found the path to it
13. Using Powershell in scripts\dist, ran `thepathIfound\pyi-archive_viewer.exe main.exe` and verified that PyQt6 was bundled
14. Did more research and found [a page from Qt for Python’s documentation](https://doc.qt.io/qtforpython-6.5/deployment/deployment-pyinstaller.html) indicating that "As of March 2021, Qt 6 is not supported yet. PyInstaller is unable to properly deploy Qt; the Qt plugins are not copied. With that, using `--onefile` is not possible." As PyQt6 is a set of bindings for Qt 6 so this may force us to use another tool than Pyinstaller. For now, I have decided to look into not using onefile. This means more than just the EXE will be generated, but I figure that we could create a MSI that saves all the other required files to the same folder.
15. Deleted the scripts\dist folder
16. Ran `poetry run pyinstaller main.py --clean` to generate a new .spec file and EXE
17. Edited the new spec file to add the venv’s path to pathex, all the hiddenimports from the previous version of the file and "PySide6" to excludes
18. Got the **qtpy** error again

### September 22, 2025

1. Provided ChatGPT with the troubleshooting steps from the two previous troubleshooting sessions and asked it for ideas
2. Added 3 lines of code to main.spec to collect submodules for "qtpy" and "PyQt6" as well as data files for "PyQt6" (see code snippets). I did not collect the dynamic libs as I found previously that it would lead to a ValueError (see September 20th’s step 7)
3. Added the collected submodules and data files to the Analysis in main.spec
4. Got the same **ValueError** as before when trying to generate the EXE: `not enough values to unpack (expected 2, got 0)`
5. Removed the datas variable from main.spec
6. The EXE is generated but when run it throws the **qtpy.QtBindingsNotFoundError** again
7. Researched QtPy’s Python 3.11 compatibility and [it should work](https://pypi.org/project/QtPy/) (with pip, not conda)
8. Ran `Invoke-Expression (poetry env activate)` on Powershell to use command line from the virtual environment, then ran `pip install qtpy`. It confirmed it was already installed
9. Asked [a question](https://stackoverflow.com/questions/79771969/encountering-qtpy-qtbindingsnotfounderror-when-running-an-exe-generated-by-pyins]) on StackOverflow in the hope of finding a fix
10. Tried to build a new EXE without changing anything (to my knowledge) and now there is a **new error**: `ERROR: Unable to find 'path_on_my_machine\\\AppData\\Local\\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\\LocalCache\local-packages\Python312\site-packages\\PyQt6\\__init__.py' when adding binary and data files.` I would have expected it to look into 
11. Tried the build again, same error
12. Looked up what should be the actual path to the venv, `path_on_my_machine\iPhyloGeo_plus_plus\.venv\Lib\site-packages\PyQt6\Qt6\translations` and it does have a qtwebengine_locales directory, which contains many .pak files. There is only one element in main.spec’s pathex and it is the site-packages folder from the virtual environment.
13. Modified the pyinstaller_build.bat file to set QTWEBENGINE_LOCALES_PATH to the appropriate path before running pyinstaller
14. Created a new file in scripts: test.py, which imports PyQt6 and then prints `PyQt6.__file__` (see code snippets)
15. In Powershell, ran `python test.py`. It printed `path_on_my_machine\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\LocalCache\local-packages\Python312\site-packages\PyQt6\__init__.py`
16. Ran `Invoke-Expression (poetry env activate)`, then `python test.py`. It showed the proper path, the one from the virtual environment.
17. From the scripts folder, ran `poetry remove PyQt6`. It said `The following packages were not found: PyQt6`.
18. Ran `poetry add PyQt6`
19. Went back to Powershell. Ran `poetry run python test.py' and it gave the correct path, within the .venv folder
20. Tried `poetry run pyinstaller test.py` and not only did it generate an EXE, but the EXE works, printing `path_on_my_machine\iPhyloGeo_plus_plus\scripts\dist\test\_internal\PyQt6\__init__.py`
21. Decided to go from the start again. Ran `poetry run pyinstaller main.py --clean`. An EXE was generated. It threw a ModuleNotFoundError when ran due to not finding qtmodern.
22. Added the path to the proper site-packages folder to main.spec’s pathex
23. Trying to generate the EXE failed as it was looking for the qtwebengine_locales folder in the wrong place again.
24. Added the following to main.spec’s hiddenimports list: "qtpy", "qtpy.QtCore", "qtpy.QtGui", "qtpy.QtWidgets", 	"PyQt6.QtCore", "PyQt6.QtGui", "PyQt6.QtWidgets" 
	"PyQt6.QtPrintSupport", "PyQt6.QtSvg", "PyQt6.QtOpenG
25. Read the log from Pyinstaller and there is an entry titled `425 INFO: Module search paths (PYTHONPATH)`.  The .venv\Lib\site-packages’s full path is listed.
26. Edited test.py to replace `import PyQt6` by `from PyQt6 import QtCore, QtGui, QtWidgets`, mimicking main.py, then print `QtCore.__file__`
27. Ran `poetry run python test.py` and it correctly pointed to the virtual environment
28. Edited test.py again to import much more (see code snippets
29. Got the following warning when running `poetry run python test.py`: `path_on_my_machine\iPhyloGeo_plus_plus\.venv\Lib\site-packages\Bio\Application\__init__.py:39: BiopythonDeprecationWarning: The Bio.Application modules and modules relying on it have been deprecated.` Then it printed the correct path, pointing to the virtual environment again. 
30. Copied all the imports from the beginning of main.py to test.py and tested again. The EXE did get generated but then failed with a **ModuleNotFoundError**: `ModuleNotFoundError: No module named 'PyQt6.QtWebEngineWidgets'.
31. Added "qtpy", "qtpy.QtCore", "qtpy.QtGui", "qtpy.QtWidgets", 	"PyQt6.QtCore", "PyQt6.QtGui", "PyQt6.QtWidgets" 
	"PyQt6.QtPrintSupport", "PyQt6.QtSv and", "PyQt6.QtOpenG to test.spec’s hiddenimports
32. Ran `poetry run python test.py` and it correctly pointed to the virtual environment
33. Added "PyQt6.QtWebEngineWidgets" to test.spec’s hiddenimports
34. Commented out the PyQt6.WebEngineWidgets import from test.py since it threw an error at runtime and I’m currently trying to troubleshoot EXE generation
35. Got another **ModuleNotFoundError**, this time with qtmodern
36. Added "qtmodern", "qtmodern.styles" and "qtmodern.windows" to test.spec’s hiddenimports
37. Ran `poetry remove qtmodern` and `poetry add qtmodern`
38. Commented out the qtmodern imports from test.py since it threw an error at runtime and I’m currently trying to troubleshoot EXE generation. Did the same for aphylogeo. Also commented out the event_connector and navigation imports
39. It generated test.exe without issue. When running it, it failed with a **ModuleNotFoundError** on PyQt6.WebEngineWidgets
40. In test.py, added QtWebEngineWidgets to the `from PyQt6 import` line. The EXE threw an **ImportError**:`ImportError: cannot import name 'QtWebEngineWidgets' from 'PyQt6' (path_on_my_machine\iPhyloGeo_plus_plus\scripts\dist\test\_internal\PyQt6\__init__.py)
41. In test.py, removed QtWebEngineWidgets from the `from PyQt import` line
42. Ran `poetry add PyQt6-WebEngine` and it was already there
43. Ran `poetry run python test.py` and it worked still
44. In test.spec, added "PyQt6" to the hiddenimports list
45. Went **back to** testing with **main.py**, and the **error** thrown when trying to generate an EXE: `Unable to find 'path_on_my_pc\\AppData\\Local\\Packages\\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0\\LocalCache\\local-packages\\Python312\\site-packages\\PyQt6\\Qt6\\translations\\qtwebengine_locales' when adding binary and data files.`

### September 28, 2025

Since testing with Pyinstaller wasn’t going well, the team decided to try CX-Freeze.

1. Ran `poetry run python setup.py build` in the scripts folder
2. When running the generated EXE, an **AttributeError** is thrown: `'NoneType object has no attribute 'isatty'`
3. Read the trace: it attributed the error to line 27 of `.venv\Lib\site-packages\toytree\utils\src\logger_setup.py`
4. Located the line in question, in the definition of a function named colorize() (see code snippets)
5. Edited logger_setup.py to mimic [this fix](https://github.com/chriskiehl/Gooey/issues/879#issuecomment-1586511649)
6. The new EXE throws a **TypeError**: `Cannot log to objects of type 'NoneType'`
7. Read the trace: like the previous error, it is linked to the toytree package, as one of the lines in the traceback is `import toytree` (see traceback below). More specifically, the line `set_log_level("WARNING")` in the package’s __init__.py file leads to the same logger_setup.py I edited at step 5. This time, the problem occurs when the logger.add function is called at line 57, as the logger is apparently passed a NoneType object.
8. Edited toytree’s __init__.py to comment out `set_log_level("WARNING")`
9. The new EXE throws a **FileNotFoundError**: `[Errno 2] No such file or directory: 'scripts/utils/params_default.yaml'`
10. Confirmed the existence of `scripts\utils\params_default.yaml`
11. Did some research online but didn’t find a fix

### September 29, 2025

1. Asked ChatGPT for advice
2. Looked in `\scripts\build\exe.win-amd64-3.11\lib\utils`: params_default.yaml is there
3. Examined the trace (see tracebacks)
4. Located the line in main.py that results in the error (line 262): it is within an if that checks if `scripts/utils/params.yaml` exists and if not, copies `scripts/utils/params.yaml` to it
5. At the beginning of main.py, added ChatGPT’s resource_path function
6. In main.py, replaced instances of `"./scripts/utils/params.yaml"` and `"scripts/utils.params.yaml"` by `resource_path(os.path.join("utils", "params.yaml"))` and did the equivalent for `params_default.yaml`
7. Tested running the app without the EXE (using a modified version of start.bat, start_pause.bat, which pauses so I can read the error messages): got a **FileNotFoundError** because it’s looking for a directory called `utils` directly in the project folder as opposed to inside scripts. 
8. Added this line to main.py: `current_dir = os.path.dirname(__file__)`
9. Changed each instance of `os.path.join("utils", "some_file")` to `os.path.join(current_dir, "utils", "some_file")`: start.bat now launches iPhyloGeo++ correctly
10. 

## TODO

1. Come up with a more permanent fix than September 28, step 5
2. Come up with a more permanent fix than September 28, step 8

## Code snippets and tracebacks

File created on September 18, step 4

```
#hook-qtmodern.py
from PyInstaller.utils.hooks import collect_data_files

datas = collect_data_files("qtmodern", includes=["**/*.qss"])

hiddenimports = [
    "qtmodern.styles",
    "qtmodern.windows",
]
```

main.spec’s hiddenimports on September 18, step 33

```
hiddenimports=["qtmodern", "qtmodern.styles", "qtmodern.windows", "PyQt6.QtCore", "PyQt6.QtGui", "PyQt6.QtWidgets", "PyQt6.QtPrintSupport", "PyQt6.QtSvg", "PyQt6.QtOpenGL", "event_connector", "numpy"],
```

File created on September 18, step 34

```
#hook-qtpy.py
hiddenimports = [
    "PyQt6",
    "PyQt6.QtCore",
    "PyQt6.QtGui",
    "PyQt6.QtWidgets",
    "PyQt6.QtPrintSupport",
]
```

File modified on September 20, step 3

```
#hooks/hook-qtpy.py
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = collect_submodules("qtpy")
```

File modified on September 20, step 4

```
#hooks/hook-qtmodern.py
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = collect_submodules("qtmoder
```

File modified on September 20, step 9
```
# main.spec

# PyInstaller spec file for Poetry + PyQt6 + qtmodern + qtpy
# Run with: poetry run pyinstaller main.spec

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
    pathex=[],
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
    debug=Fal
    bootloader_ignore_signals=False,
    strip=False,
    upx=[],
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
)
```

Lines added on September 22, step 2:

```
from PyInstaller.utils.hooks import collect_submodules, collect_data_files, collect_dynamic_libs

hiddenimports = collect_submodules("qtpy") + collect_submodules("PyQt6")
datas = collect_data_files("PyQt6")
```

File created on September 22, step 14:

```
#test.py
import PyQt6

print(PyQt6.__file__)
```

File modified on September 22, step 28:

```
#test.py
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QColor, QIcon
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import QGraphicsDropShadowEffect, QVBoxLayout
from Qt import main_ui
from ui_controllers import ClimatePageController, GeneticPageController, ResultPageController
from ui_helpers import create_shadow_effect, get_button_style, style_buttons
from utils import resources_rc  # noqa: F401  # Import the compiled resource module for resolving image resource path
from utils.error_dialog import show_error_dialog

print(QtCore.__file__)
```

Function found on September 28, step 4

```
def colorize():
    """colorize the logger if stderr is IPython/Jupyter or a terminal (TTY)"""
    try:
        import IPython
        tty1 = bool(IPython.get_ipython())
    except ImportError:
        tty1 = False
    tty2 = sys.stderr.isatty() # line throwing the error
    if tty1 or tty2:
        return True
    return False
```

Full traceback from September 28, step 7

```
Traceback (most recent call last):
  File "path_on_my_machine\iPhyloGeo_plus_plus\.venv\Lib\site-packages\cx_Freeze\initscripts\__startup__.py", line 133, in run
    module_init.run(f"__main__{name}")
  File "path_on_my_machine\iPhyloGeo_plus_plus\.venv\Lib\site-packages\cx_Freeze\initscripts\console.py", line 25, in run
    exec(code, main_globals)
  File "main.py", line 18, in <module>
  File "path_on_my_machine\iPhyloGeo_plus_plus\scripts\ui_controllers\__init__.py", line 2, in <module>
    from .genetics_page_controller import GeneticPageController
  File "path_on_my_machine\iPhyloGeo_plus_plus\scripts\ui_controllers\genetics_page_controller.py", line 2, in <module>
    from Genetics.genetics import Genetics
  File "path_on_my_machine\iPhyloGeo_plus_plus\scripts\Genetics\genetics.py", line 5, in <module>
    from Genetics.genetics_tree import GeneticTree
  File "path_on_my_machine\iPhyloGeo_plus_plus\scripts\Genetics\genetics_tree.py", line 4, in <module>
    from utils.tree_graph import generate_tree, init_tree
  File "path_on_my_machine\iPhyloGeo_plus_plus\scripts\utils\tree_graph.py", line 9, in <module>
    import toytree
  File "path_on_my_machine\iPhyloGeo_plus_plus\.venv\Lib\site-packages\toytree\__init__.py", line 44, in <module>
    set_log_level("WARNING")
  File "path_on_my_machine\iPhyloGeo_plus_plus\.venv\Lib\site-packages\toytree\utils\src\logger_setup.py", line 57, in set_log_level
    idx = logger.add(
          ^^^^^^^^^^^
  File "path_on_my_machine\iPhyloGeo_plus_plus\.venv\Lib\site-packages\loguru\_logger.py", line 872, in add
    raise TypeError("Cannot log to objects of type '%s'" % type(sink).__name__)
TypeError: Cannot log to objects of type 'NoneType'
```

Full traceback from September 29, step 3

```
Traceback (most recent call last):
  File "path_and_project_folder\.venv\Lib\site-packages\cx_Freeze\initscripts\__startup__.py", line 133, in run
    module_init.run(f"__main__{name}")
  File "path_and_project_folder\.venv\Lib\site-packages\cx_Freeze\initscripts\console.py", line 25, in run
    exec(code, main_globals)
  File "main.py", line 262, in <module>
  File "C:\Users\agaco\AppData\Local\Programs\Python\Python311\Lib\shutil.py", line 419, in copy
    copyfile(src, dst, follow_symlinks=follow_symlinks)
  File "C:\Users\agaco\AppData\Local\Programs\Python\Python311\Lib\shutil.py", line 256, in copyfile
    with open(src, 'rb') as fsrc:
         ^^^^^^^^^^^^^^^
FileNotFoundError: [Errno 2] No such file or directory: 'scripts/utils/params_default.yaml'
```
