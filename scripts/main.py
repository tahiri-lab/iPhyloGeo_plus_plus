import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import pathlib
from howtouse import Ui_sec
#from cltree import Ui_ct
import toytree
import random
import toyplot.pdf



class Ui_Main_Window(object):
    #[

    #def openWindow1(self):
     #   self.window = QtWidgets.QMainWindow()
      #  self.ui = Ui_sec()
       # self.ui.setupUi(QtWidgets.QMainWindow())
        #self.window.show()
        #self.howtouse.setText(self.ui)


    def useWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_sec()
        self.ui.setupUi(self.window)
        self.window.show()



# open cltree window
    def openWindow(self):
        from cltree import Ui_ct
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_ct()
        self.ui.setupUi(self.window)
        self.window.show()        
        
# for clim_tree
    #def openWindow1(self):
        #def openWindow(self):
     #   self.window1 = QtWidgets.QMainWindow()
        #toyplot.pdf.render(canvas,'../proj/climactic_trees.pdf')
        #import tree1
        #self.ui.setupUi(self.window1)
      #  self.window1.show()
     

    #]

    def setupUi(self, Main_Window):
        Main_Window.setObjectName("Main_Window")
        Main_Window.resize(1400, 1000)
        font = QtGui.QFont()
        font.setFamily("aakar")
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        Main_Window.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("p1zbTlLQ_400x400.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon.addPixmap(QtGui.QPixmap("icons8-sequence-64.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Main_Window.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(Main_Window)
        self.centralwidget.setObjectName("centralwidget")
        
#[
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(5, 0, 1350, 950))
        #self.tabWidget.setGeometry(QtCore.QRect(10, 30, 691, 461))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
#]
        self.seq_alig = QtWidgets.QPushButton(self.tab)
        #self.seq_alig = QtWidgets.QPushButton(self.centralwidget)
        self.seq_alig.setGeometry(QtCore.QRect(15, 200, 141, 51))
        #self.seq_alig.setGeometry(QtCore.QRect(10, 150, 141, 51))        
        self.seq_alig.setIcon(icon)
        self.seq_alig.setObjectName("seq_alig")
        self.browse = QtWidgets.QPushButton(self.tab, clicked = lambda: self.press_it())
        #self.browse = QtWidgets.QPushButton(self.centralwidget, clicked = lambda: self.press_it())
        #self.browse.setGeometry(QtCore.QRect(10, 10, 141, 51))
        self.browse.setGeometry(QtCore.QRect(15, 60, 141, 51))
        self.browse.setObjectName("browse")
        self.clear = QtWidgets.QPushButton(self.tab, clicked = lambda: self.clear_it())
        #self.clear = QtWidgets.QPushButton(self.centralwidget, clicked = lambda: self.clear_it())
        #self.clear.setGeometry(QtCore.QRect(10, 80, 141, 51))
        self.clear.setGeometry(QtCore.QRect(15, 130, 141, 51))
        self.clear.setObjectName("clear")

        # gens data
        self.te = QtWidgets.QTextEdit(self.tab)
        #self.te = QtWidgets.QTextEdit(self.centralwidget)
        self.te.setGeometry(QtCore.QRect(210, 60, 1000, 351))
    # Font type and size[
        font = QtGui.QFont()
        font.setFamily("Tibetan Machine Uni")
        font.setPointSize(10)
        self.te.setFont(font)
     #]
        self.te.setObjectName("te")
#[      seq_alig gens data
        self.te1 = QtWidgets.QTextEdit(self.tab)
        #self.te = QtWidgets.QTextEdit(self.centralwidget)
        self.te1.setGeometry(QtCore.QRect(210, 450, 1000, 351))
        self.te1.setObjectName("te1")
  

  #now { # climatic data text box
        self.te2 = QtWidgets.QTextEdit(self.tab_2)
        #self.te2 = QtWidgets.QTextEdit(self.centralwidget)
        self.te2.setGeometry(QtCore.QRect(210, 60, 1000, 351))
        self.te2.setObjectName("te2")

 #[now
        # button to draw tree
        #self.to_draw = QtWidgets.QPushButton(self.tab_2, clicked =  lambda: self.to_draw())
        self.to_draw = QtWidgets.QPushButton(self.tab_2, clicked = lambda: self.openWindow())
        #self.clear.setGeometry(QtCore.QRect(10, 80, 141, 51))
        self.to_draw.setGeometry(QtCore.QRect(15, 130, 141, 51))
        self.to_draw.setObjectName("to_draw")
      #]
        #}
        self.browse_c = QtWidgets.QPushButton(self.tab_2, clicked = lambda: self.pressit())
        #self.browse = QtWidgets.QPushButton(self.centralwidget, clicked = lambda: self.press_it())
        #self.browse.setGeometry(QtCore.QRect(10, 10, 141, 51))
        self.browse_c.setGeometry(QtCore.QRect(15, 60, 141, 51))
        self.browse_c.setObjectName("browse")
#]
        Main_Window.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(Main_Window)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1800, 25))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuEdit = QtWidgets.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        self.menuAbout = QtWidgets.QMenu(self.menubar)
        self.menuAbout.setObjectName("menuAbout")
        Main_Window.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(Main_Window)
        self.statusbar.setObjectName("statusbar")
        Main_Window.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(Main_Window)
        self.toolBar.setObjectName("toolBar")
        Main_Window.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionOpen = QtWidgets.QAction(Main_Window, triggered = lambda: self.press_it())
        self.actionOpen.setObjectName("actionOpen")
        self.actionCopy = QtWidgets.QAction(Main_Window)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("icons/blue-document-copy.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionCopy.setIcon(icon1)
        self.actionCopy.setObjectName("actionCopy")
        self.actionCut = QtWidgets.QAction(Main_Window)
        self.actionCut.setObjectName("actionCut")
        self.actionRead_Me = QtWidgets.QAction(Main_Window)
        self.actionRead_Me.setObjectName("actionRead_Me")
        self.actionHelp = QtWidgets.QAction(Main_Window, triggered = lambda: self.useWindow())#.........
        self.actionHelp.setObjectName("actionHelp")
        self.actionopen = QtWidgets.QAction(Main_Window)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("icons/folder-open.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionopen.setIcon(icon2)
        self.actionopen.setObjectName("actionopen")
        self.menuFile.addAction(self.actionOpen)
        self.menuEdit.addAction(self.actionCopy)
        self.menuEdit.addAction(self.actionCut)
        self.menuAbout.addAction(self.actionRead_Me)
        self.menuAbout.addAction(self.actionHelp)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())
        self.toolBar.addAction(self.actionopen)
        self.toolBar.addAction(self.actionCopy)

        self.retranslateUi(Main_Window)

#[
        self.tabWidget.setCurrentIndex(0)
#]
        QtCore.QMetaObject.connectSlotsByName(Main_Window)

    
    # press the button to get gens data
    def press_it(self):
        filename = QFileDialog.getOpenFileName()
        path = pathlib.Path(filename[0])
        okPressed = True
        with open(path) as f:
            #print(f.read())
            #print(path)
            self.te.setText(f.read())

  
    # press the button to get climatic data
    def pressit(self):
        filename = QFileDialog.getOpenFileName()
        path = pathlib.Path(filename[0])
        okPressed = True
        with open(path) as c:
            #print(f.read())
            #print(path)
            self.te2.setText(c.read())

#[now
    # press the button to draw climatic tree
   # def to_draw(self):
     #   import tree1
    #    self.te2.clear()
        #filename = QFileDialog.getOpenFileName()
        #path = pathlib.Path(filename[0])
        #okPressed = True
        #with open(path) as ff:
            #print(f.read())
            #print(path)
         #   self.te2.setText(ff.read())
    
#]

      # press the button to delet data
    def clear_it(self):
        self.te.clear()


    def retranslateUi(self, Main_Window):
        _translate = QtCore.QCoreApplication.translate
        Main_Window.setWindowTitle(_translate("Main_Window", "aPhylogeo"))
#[
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "GENS DATA"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "CLIMATIC DATA"))
#]


        self.seq_alig.setToolTip(_translate("Main_Window", "<html><head/><body><p>sequence alignment</p></body></html>"))
        self.seq_alig.setText(_translate("Main_Window", "  Seq_Alig"))
        self.browse.setText(_translate("Main_Window", "Browse File"))
        self.browse_c.setText(_translate("Main_Window", "Browse \n Climatic Data"))
        self.to_draw.setText(_translate("Main_Window", "Draw \n Climatic Tree "))
        self.clear.setText(_translate("Main_Window", "Clear"))
        self.menuFile.setTitle(_translate("Main_Window", "File"))
        self.menuEdit.setTitle(_translate("Main_Window", "Edit"))
        self.menuAbout.setTitle(_translate("Main_Window", "Help"))
        self.toolBar.setWindowTitle(_translate("Main_Window", "toolBar"))
        self.actionOpen.setText(_translate("Main_Window", "Open"))
        self.actionCopy.setText(_translate("Main_Window", "Copy"))
        self.actionCut.setText(_translate("Main_Window", "Cut"))
        self.actionRead_Me.setText(_translate("Main_Window", "Read Me"))
        self.actionHelp.setText(_translate("Main_Window", "How to use"))
        self.actionopen.setText(_translate("Main_Window", "open"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Main_Window = QtWidgets.QMainWindow()
    ui = Ui_Main_Window()
    ui.setupUi(Main_Window)
    Main_Window.show()
    sys.exit(app.exec_())
