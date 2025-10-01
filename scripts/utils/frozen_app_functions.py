import os
import sys

def find_utils(filename):
    if getattr(sys, 'frozen', False):
        # The application is frozen
        datadir = os.path.join(os.path.dirname(sys.executable), "lib", "utils")
    else:
        # The application is not frozen
        datadir = os.path.dirname(__file__)
    return os.path.join(datadir, filename)