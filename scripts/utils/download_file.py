import os
import shutil

import matplotlib.pyplot as plt
import plotly.io as pio
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QFileDialog
from utils.error_dialog import show_error_dialog


def download_file_temporary_PLT(name, fig=None):
    img_path = format_string_create_directory(name)
    plt.savefig(img_path)
    plt.close(fig)
    return QPixmap(img_path)


def download_file_temporary_PIO(name, fig):
    img_path = format_string_create_directory(name)
    pio.write_image(fig, img_path, format="png")
    return QPixmap(img_path)


def download_file_local(name, parent):
    img_path = format_string_create_directory(name)

    if not os.path.exists(img_path):
        show_error_dialog("No plot found to download.")
        return

    # Prompt the user to select a location to save the plot
    file_path, _ = QFileDialog.getSaveFileName(
        parent,
        "Save Plot As",
        os.path.basename(img_path),
        "PNG Files (*.png);;All Files (*)",
    )
    if file_path:
        if not file_path.lower().endswith(".png"):
            file_path += ".png"
        shutil.copy(img_path, file_path)


results_dir = "results"


def format_string_create_directory(name):
    img_path = os.path.join(results_dir, f"{name.lower().replace(' ', '_')}.png")
    os.makedirs(results_dir, exist_ok=True)
    return img_path
