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

find_string_in_scripts("Ui_MainWindow")
