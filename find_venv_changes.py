import os
import time
from datetime import datetime

# The script will look for changes after this timestamp
datetime_str = "2025-12-19 12:08:00"

datetime_obj = datetime.strptime(datetime_str, "%Y-%m-%d %H:%M:%S")
check_ti = datetime_obj.timestamp()

def print_if_later(path, check_ti):
    # Both the variables would contain time
    # elapsed since EPOCH in float
    ti_c = os.path.getctime(path)
    ti_m = os.path.getmtime(path)

    # Converting the time in seconds to a timestamp
    c_ti = time.ctime(ti_c)
    m_ti = time.ctime(ti_m)

    if (ti_m >= check_ti):
        print(f"{path} \
        was created at {c_ti} and was "
              f"last modified at {m_ti}")


for root, _, files in os.walk(".venv"):
    for file in files:
        file_path = os.path.join(root, file)
        print_if_later(file_path, check_ti)

#print_if_later(".venv/.gitignore", check_ti)