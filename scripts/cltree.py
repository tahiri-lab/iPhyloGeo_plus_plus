
import os
import PyPDF2
from PyQt5 import QtCore, QtGui, QtWidgets
import tree
import toytree
import random
import toyplot.pdf
import sys

class Ui_ct(object):
    def setupUi(self, ct):
        ct.setObjectName("ct")
        ct.resize(1258, 1100)
        self.centralwidget = QtWidgets.QWidget(ct)
        self.centralwidget.setObjectName("centralwidget")
        self.tb = QtWidgets.QTextBrowser(self.centralwidget)#, triggered = lambda: self.openWindow())
        self.tb.setGeometry(QtCore.QRect(20, 20, 1200, 1031))
        self.tb.setObjectName("tb")
        ct.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(ct)
        self.statusbar.setObjectName("statusbar")
        ct.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(ct)
        self.toolBar.setObjectName("toolBar")
        ct.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionPrint = QtWidgets.QAction(ct)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("print.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPrint.setIcon(icon)
        self.actionPrint.setObjectName("actionPrint")
        self.toolBar.addAction(self.actionPrint)

    
        self.retranslateUi(ct)
        QtCore.QMetaObject.connectSlotsByName(ct)

    def retranslateUi(self, ct):
        _translate = QtCore.QCoreApplication.translate
        ct.setWindowTitle(_translate("ct", "Climate Tree"))
        
        self.tb.setHtml(_translate("tb", '../test/climatic_trees.pdf'))
        self.toolBar.setWindowTitle(_translate("ct", "toolBar"))
        self.actionPrint.setText(_translate("ct", "Print"))


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    ct = QtWidgets.QMainWindow()
    ui = Ui_ct()
    ui.setupUi(ct)
    ct.show()
    sys.exit(app.exec_())



    
