from PyQt6 import QtWidgets


def show_error_dialog(message, title="error"):
    """
    Display a professional error dialog with the given title and message.

    Args:
        title (str): The title of the error dialog.
        message (str): The error message to display.
    """
    msgBox = QtWidgets.QMessageBox()
    msgBox.setIcon(QtWidgets.QMessageBox.Icon.Critical)
    msgBox.setWindowTitle(title)
    msgBox.setText(message)
    msgBox.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
    msgBox.setDefaultButton(QtWidgets.QMessageBox.StandardButton.Ok)
    msgBox.exec()
