from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(400, 130)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.waitingLabel = QtWidgets.QLabel(Dialog)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.waitingLabel.setFont(font)
        self.waitingLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.waitingLabel.setObjectName("waitingLabel")
        self.verticalLayout.addWidget(self.waitingLabel)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.progressBar = QtWidgets.QProgressBar(Dialog)
        self.progressBar.setMinimumSize(QtCore.QSize(0, 20))
        self.progressBar.setStyleSheet("QProgressBar {\n"
"    border: 2px solid #5A5A5A;\n"
"    border-radius: 5px;\n"
"    background-color: #E0E0E0;\n"
"}\n"
"QProgressBar::chunk {\n"
"    background-color: #3498DB;\n"
"    width: 20px;\n"
"}")
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(0)
        self.progressBar.setAlignment(QtCore.Qt.AlignCenter)
        self.progressBar.setTextVisible(False)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout.addWidget(self.progressBar)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Please Wait"))
        self.waitingLabel.setText(_translate("Dialog", "Please wait..."))
