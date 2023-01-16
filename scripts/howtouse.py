

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_sec(object):
    def setupUi(self, sec):
        sec.setObjectName("sec")
        sec.resize(1258, 1134)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(sec.sizePolicy().hasHeightForWidth())
        sec.setSizePolicy(sizePolicy)
        self.centralwidget = QtWidgets.QWidget(sec)
        self.centralwidget.setObjectName("centralwidget")
        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(20, 20, 1211, 1031))
        self.textBrowser.setObjectName("textBrowser")
        sec.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(sec)
        self.statusbar.setObjectName("statusbar")
        sec.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(sec)
        self.toolBar.setObjectName("toolBar")
        sec.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionprint = QtWidgets.QAction(sec)
        self.actionprint.setObjectName("actionprint")
        self.actionprint_2 = QtWidgets.QAction(sec)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("print.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionprint_2.setIcon(icon)
        self.actionprint_2.setObjectName("actionprint_2")
        self.toolBar.addAction(self.actionprint_2)

        self.retranslateUi(sec)
        QtCore.QMetaObject.connectSlotsByName(sec)

    def retranslateUi(self, sec):
        _translate = QtCore.QCoreApplication.translate
        sec.setWindowTitle(_translate("sec", "How To Use"))
        self.textBrowser.setHtml(_translate("sec", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:12pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">How to use the application</p></body></html>"))
        self.toolBar.setWindowTitle(_translate("sec", "toolBar"))
        self.actionprint.setText(_translate("sec", "print"))
        self.actionprint_2.setText(_translate("sec", "print"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    sec = QtWidgets.QMainWindow()
    ui = Ui_sec()
    ui.setupUi(sec)
    sec.show()
    sys.exit(app.exec_())
