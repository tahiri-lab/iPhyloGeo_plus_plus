import threading
import sys
import time
import io

class InterceptOutput(io.StringIO):
    def __init__(self, g_func, original_stdout):
        super().__init__()
        self.g_func = g_func
        self.original_stdout = original_stdout

    def write(self, s):
        if s.strip():  # Ignore empty lines
            self.original_stdout.write(s)  # Ensure the original print still happens
            self.original_stdout.flush()  # Make sure it gets printed immediately
            self.g_func(s.strip())

def f():
    for i in range(5):
        print(f"Output {i}")
        time.sleep(1)  # Simulating a long-running process

def g(message):
    # Simulate sending the message to another microservice
    print(f"g received: {message}", file=sys.__stdout__)  # Print to the original stdout

def run_f_and_intercept_output():
    original_stdout = sys.stdout
    try:
        sys.stdout = InterceptOutput(g, original_stdout)
        f()
    finally:
        sys.stdout = original_stdout

if __name__ == "__main__":
    thread = threading.Thread(target=run_f_and_intercept_output)
    thread.start()
    thread.join()
