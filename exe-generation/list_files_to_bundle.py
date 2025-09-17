# This script creates a list of the files that must explicitely be bundled when generating an EXE with Pyinstaller

import os

def list_non_python_files(folder_path):
    """
    Print local path and filename for all files in folder_path
    that do not end with .py or .pyc.
    """
    for root, dirs, files in os.walk(folder_path):
        #we must exclude certain folders which are not part of iPhyloGeo++ itself
        if ".git" in dirs:
            dirs.remove(".git")
        if ".github" in dirs:
            dirs.remove(".github")
        if ".venv" in dirs:
            dirs.remove(".venv")
        if ".vscode" in dirs:
            dirs.remove(".vscode")
        for file in files:
            excluded = False
            #we must exclude .py scripts since they will be automatically bundled by pyinstaller
            #we must exclude .pyc caches since they are not necessary to run iPhyloGeo++
            if (file.endswith(".py") or file.endswith(".pyc")):
                excluded = True
            #we must exclude .DS_Store and .gitignore
            if (file == ".DS_Store" or file == ".gitignore"):
                excluded = True
            if not (excluded):
                full_path = os.path.join(root, file)
                print(full_path[2:])

if __name__ == "__main__":
    folder = ".."
    list_non_python_files(folder)

