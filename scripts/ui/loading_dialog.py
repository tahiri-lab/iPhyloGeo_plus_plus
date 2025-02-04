from PyQt6 import QtWidgets

from scripts.Qt.loading_ui import Ui_LoadingDialog


class LoadingDialog(Ui_LoadingDialog, QtWidgets.QDialog):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
