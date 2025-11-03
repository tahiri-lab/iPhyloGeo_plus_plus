# Going back on temporary fixes

When troubleshooting EXE generation, I employed two workarounds:

1. I edited logger_setup.py to mimic [this fix](https://github.com/chriskiehl/Gooey/issues/879#issuecomment-1586511649)
2. I edited toytree’s __init__.py to comment out `set_log_level("WARNING")`

Both workarounds are related to toytree and neither looks like a good long term solution.

# November 3, 2025

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



TODO

* Remove the license from setup.py since it doesn’t work and the option is not recognized
* Make a declaration for AI use on october 11th (Claude Sonnet 4.5) [https://claude.ai/share/492447c6-c91d-40f2-ba96-48dff6b21b6f]


## Original notes about the workarounds

1. Ran `poetry run python setup.py build` in the scripts folder
2. When running the generated EXE, an **AttributeError** is thrown: `'NoneType object has no attribute 'isatty'`
3. Read the trace: it attributed the error to line 27 of `.venv\Lib\site-packages\toytree\utils\src\logger_setup.py`
4. Located the line in question, in the definition of a function named colorize() (see code snippets)
5. Edited logger_setup.py to mimic [this fix](https://github.com/chriskiehl/Gooey/issues/879#issuecomment-1586511649)
6. The new EXE throws a **TypeError**: `Cannot log to objects of type 'NoneType'`
7. Read the trace: like the previous error, it is linked to the toytree package, as one of the lines in the traceback is `import toytree` (see traceback below). More specifically, the line `set_log_level("WARNING")` in the package’s `__init__.py` file leads to the same logger_setup.py I edited at step 5. This time, the problem occurs when the logger.add function is called at line 57, as the logger is apparently passed a NoneType object.
8. Edited toytree’s `__init__.py` to comment out `set_log_level("WARNING")`