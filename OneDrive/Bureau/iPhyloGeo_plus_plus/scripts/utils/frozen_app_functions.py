# This file contains functions which were created in order for iPhyloGeo++ to work not only as a live app, but also as a frozen one.

import os
import sys

def find_utils(filename):
    """
    Method that provides the correct path to a file within the utils directory.
    args:
        filename (str): the name of the file, including the extension
    """
    if getattr(sys, 'frozen', False):
        # The application is frozen
        datadir = os.path.join(os.path.dirname(sys.executable), "lib", "utils")
    else:
        # The application is not frozen
        datadir = os.path.dirname(__file__)
    return os.path.join(datadir, filename)