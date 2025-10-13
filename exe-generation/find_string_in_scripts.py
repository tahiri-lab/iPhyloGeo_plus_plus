# This script asks the user for a string and prints every line of code in the scripts folder containing that string.
# It makes troubleshooting easier as EXE errors may be unclear.
# Example:
# - User inputs "frozen_app_functions"
# - The program prints a line for each occurence of "frozen_app_functions", like so:
#   "Found in ..\scripts\main.py (Line 22): from utils.frozen_app_functions import find_utils"

import os

def find_string_in_scripts(search_string):
    for root, _, files in os.walk("..\\scripts"):
        for file in files:
            if file.endswith(".py"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    try:
                        for line_number, line in enumerate(f, start=1):
                            if search_string in line:
                                print(f"Found in {file_path} (Line {line_number}): {line.strip()}")
                    except UnicodeDecodeError:
                        print(f"Could not read {file_path} due to encoding issues.")

print("Input the string to search for, then press the ENTER key.")
find_string_in_scripts(input())