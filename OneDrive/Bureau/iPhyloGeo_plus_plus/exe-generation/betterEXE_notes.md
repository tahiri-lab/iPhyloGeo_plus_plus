# Improving the EXE

## Tasks

### Going back on temporary fixes

When troubleshooting EXE generation, I employed two workarounds:

1. I edited logger_setup.py to mimic [this fix](https://github.com/chriskiehl/Gooey/issues/879#issuecomment-1586511649)
2. I edited toytree’s __init__.py to comment out `set_log_level("WARNING")`

Both workarounds are related to toytree and neither looks like a good long term solution. The second one was a workaround for the first one not fully working.

### Making the app faster

The app is somewhat slow, both in its original version and its frozen one.

## Work

### November 3, 2025

First, I had to get back to the initial issue.

1. Tried to access the main branch. It says a **lock file** is preventing it but there is no .git/index.lock
2. Cloned iPhyloGeo_plus_plus anew
3. Ran exe-generation\BUILD_frozenapp.bat: it **can’t find cx_Freeze** because I forgot to install it in the poetry environment after getting the repo
4. Ran `poetry run pip install cx_freeze'
5. Ran BUILD_frozenapp.bat: it completed (very fast) but generated an EXE which threw a **ModuleNotFoundError**: `No module named 'qtmodern'`
6. Realized I hadn’t built the poetry environment (I reimported the entire project, after all) and ran `poetry config virtualenvs.in-project true`, then `poetry install`
7. Ran BUILD_frozenapp.bat: it took a VERY long time and seemed to be **hanging** so I interrupted it
8. Ran BUILD_frozenapp.bat again with the same result
9. Ran `poetry run python --version` and it said 3.12.10 so I need to either go back to 3.11 or get into major compatibility issues
10. Ran `poetry env use 3.11` and updated pyproject.toml to reflect the change to 3.11
11. Ran `poetry lock`, then `poetry install`
12. Ran BUILD_frozenapp.bat: the EXE threw the **ModuleNotFoundError** again: `No Module named 'qtmodern'
13. Ran `poetry run pip install qtmodern`
14. Ran the build again: when running the EXE I’m back at the **AttributeError** I worked around previously. `'NoneType' object has no attribute 'isatty'`

```
---------------------------
cx_Freeze: Python error in main script
---------------------------
Traceback (most recent call last):
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\.venv\Lib\site-packages\cx_Freeze\initscripts\__startup__.py", line 133, in run
    module_init.run(f"__main__{name}")
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\.venv\Lib\site-packages\cx_Freeze\initscripts\console.py", line 25, in run
    exec(code, main_globals)
  File "main.py", line 18, in <module>
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\scripts\ui_controllers\__init__.py", line 2, in <module>
    from .genetics_page_controller import GeneticPageController
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\scripts\ui_controllers\genetics_page_controller.py", line 2, in <module>
    from Genetics.genetics import Genetics
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\scripts\Genetics\genetics.py", line 5, in <module>
    from Genetics.genetics_tree import GeneticTree
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\scripts\Genetics\genetics_tree.py", line 4, in <module>
    from utils.tree_graph import generate_tree, init_tree
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\scripts\utils\tree_graph.py", line 9, in <module>
    import toytree
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\.venv\Lib\site-packages\toytree\__init__.py", line 44, in <module>
    set_log_level("WARNING")
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\.venv\Lib\site-packages\toytree\utils\src\logger_setup.py", line 57, in set_log_level
    colorize=colorize(),
             ^^^^^^^^^^
  File "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\.venv\Lib\site-packages\toytree\utils\src\logger_setup.py", line 27, in colorize
    tty2 = sys.stderr.isatty()
           ^^^^^^^^^^^^^^^^^
AttributeError: 'NoneType' object has no attribute 'isatty'

---------------------------
OK   
---------------------------
```

This time, I am looking for a solution that won’t require editing any of the files that are part of the virtual environment. I don’t want to make any changes to cx_Freeze or toytree.

Last resort, I’ll edit the BUILD and BUILDwMSI scripts to edit the files within the venv, but first I’ll look for a more elegant solution.

1. Read (this)[https://stackoverflow.com/questions/50727131/when-is-the-attribute-isatty-not-present-in-sys-stderr] and (this)[https://stackoverflow.com/questions/50727131/when-is-the-attribute-isatty-not-present-in-sys-stderr]. It is possible our app, having a GUI, redirects the standard output, which messes with how toytree deals with error
2. Read (this)[https://forum.derivative.ca/t/fixed-stdoutcatcher-object-has-no-attribute-isatty-when-importing-python-module/145477], which suggests a fix
3. Looked at the traceback. The issue occurs in .venv\Lib\site-packages\toytree\utils\src\logger_setup.py, when it calls the colorize function. I read that function. The line that results in the the error must be `tty2 = sys.stderr.isatty()`, because sys.stderr may not have "isatty" defined. According to the conversation I read at step 2, the workaround would be the subclass stdout with that function defined.
4. Asked Claude (LLM) for advice on how to do that. It suggested creating if statements checking if the attribute "isatty" exists for sys.stderr, sys.stdout and sys.stdin, and, when it doesn’t, setting it to `lambda: False`.
5. Used my exe-generation find_string_in_scripts.py and found that only scripts\utils\tree_graph.py imports toytree
6. Added Claude’s suggested patch to tree_graph.py
7. Ran the build: the EXE still throws the **AttributeError**, now on the patch’s first line
8. Asked Claude to ask something to the patch to account for this
9. Applied Claude’s new patch
10. Ran the build: success!

This is good. Now, I’d like to move on to making the app faster. First, I’ll need to decide how to keep track of how fast the app is.

Did some reading:
* (Introduction to Windows Application Performance)[https://learn.microsoft.com/en-us/windows/apps/performance/introduction]
* (Windows app performance and fundamentals overview)[https://learn.microsoft.com/en-us/windows/apps/performance/]
* (Choosing among Visual Studio Performance Profiler, Windows Performance Toolkit, and PerfView)[https://learn.microsoft.com/en-us/windows/apps/performance/choose-between-tools]

My priorities are the following:
* Minimize the time between launching the app and getting the UI
* Reduce CPU and memory usage

### November 7, 2025

Did some more reading, this time focused on optimizing and benchmarking Python code:
* (Python Benchmarking: Unleashing the Power of Performance Analysis)[https://coderivers.org/blog/python-benchmark/]
* (4 Ways to Benchmark Python Code)[https://superfastpython.com/benchmark-python-code/]
* (Profiling in Python: How to Find Performance Bottlenecks)[https://realpython.com/python-profiling/]

1. Created start-optimized.bat, an alternative to start.bat that uses Python’s optimize mode.
2. In scripts\setup.py, set optimize to 2 instead of 1
3. Ran the frozen app build
4. Ran the EXE: it throws a **TypeError**: `unsupported operand type(s) for %: 'NoneType' and 'str'`
5. Read the traceback: the ete3 module is trying to set a docstring (optimize being set to 2 disables docstrings) for an evolutionary model
6. Asked Claude for a patch
7. Went back to the traceback to determine where the patch should be applied: aphylogeo.utils imports ete3, so the patch should be applied before importing aphylogeo.utils
8. Used exe-generation\find_string_in_scripts.py to locate the relevant imports: scripts\result.py line 2, scripts\worker.py line 4 and scripts\Climatic\climat.py line 6
9. Adapted and saved the patch suggested by Claude as scripts\utils\oo_patch
10. Updated the three files identified at step 8 to use the patch
11. Tried the live app: it throws an **AttributeError**: `'dict' object has no attribute '__build_class__'` on a line from the patch when trying to import said patch
12. Informed Claude, which provided an updated version of the patch
13. Updated the patch and tried the live app: it threw a **TypeError**: `metaclass conflict: the metaclass of a derived class must be a (non-strict) subclass of the metaclasses of all its bases`
14. Worked further with Claude before determining that the -OO option was not worth it and optimizing the app as well as aphylogeo would be a better use of my time

### November 19, 2025

Working on making aPhyloGeo faster instead.

1. In aphylogeo\tests, ran both `pytest .\test_climatic.py` and `pytest .\test_genetic.py`. Both ran. One **test failed**: `FAILED test_genetic.py::TestGenetic::test_filterResults - assert "Phylogeny(rooted=True)\n    Clade()\n        Clade(branch_length=0.0, name='OL989074')\n        Clade(bran...`
2. Ran it again: same result, so one of the tests fails on Python 3.12, I’ll be forced to use Python 3.11 until Khady fixes it.
3. Set up the poetry venv
4. Ran `poetry run pytest .\test_climatic.py` and `poetry run pytest .\test_genetic.py` and the same test failed again, which is definitely not meant to happen
5. Created a new unit test file with only the problematic test in it, test_investigation.py, and ran pytest: same result
6. Ran it from the main project directory: `poetry run pytest tests\test_investigation.py`: same result
7. Realized the environment still used Python 3.12 and changed it to Python 3.11
8. Located my install: C:\Users\agaco\AppData\Local\Programs\Python\Python311
9. Ran `poetry env use C:\Users\agaco\AppData\Local\Programs\Python\Python311\python.exe`, then `poetry install`, then ran the tests again: it still failed
10. Something is wrong, venv not working properly, `poetry install` fails, will come back to it when I have more energy