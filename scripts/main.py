import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QApplication, QMainWindow, QTableWidget, QTableWidgetItem
import pathlib
from howtouse import Ui_sec
from cltree import Ui_ct 
import toytree
import random
import toyplot.pdf
import pandas as pd
from aPhyloGeo.aPhyloGeo import create_and_save_tree


class Ui_Main_Window(object):
   

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
        

    def setupUi(self, Main_Window):
        Main_Window.setObjectName("Main_Window")
        Main_Window.resize(1400, 1000)
        font = QtGui.QFont()
        font.setFamily("aakar")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(50)
        Main_Window.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("p1zbTlLQ_400x400.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon.addPixmap(QtGui.QPixmap("icons8-sequence-64.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Main_Window.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(Main_Window)
        self.centralwidget.setObjectName("centralwidget")
           
             
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(5, 0, 1350, 950))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")

        self.comboBox = QtWidgets.QComboBox(self.tab_2)
        self.comboBox.setGeometry(QtCore.QRect(420, 480, 201, 30))
        
        #self.comboBox.setSizePolicy(sizePolicy)
        self.comboBox.setFrame(True)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.frame_3 = QtWidgets.QFrame(self.tab_2)
        self.frame_3.setGeometry(QtCore.QRect(730, 460, 381, 61))
        self.frame_3.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.comboBox_2 = QtWidgets.QComboBox(self.frame_3)
        self.comboBox_2.setGeometry(QtCore.QRect(170, 20, 201, 30))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.comboBox_2.sizePolicy().hasHeightForWidth())
        self.comboBox_2.setSizePolicy(sizePolicy)
        self.comboBox_2.setFrame(True)
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.label_16 = QtWidgets.QLabel(self.frame_3)
        self.label_16.setGeometry(QtCore.QRect(80, 20, 121, 31))
        self.label_16.setObjectName("label_16")
        self.label_13 = QtWidgets.QLabel(self.tab_2)
        self.label_13.setGeometry(QtCore.QRect(340, 480, 201, 30))
        self.label_13.setObjectName("label_13")
        
        self.seq_alig = QtWidgets.QPushButton(self.tab)
        self.seq_alig.setGeometry(QtCore.QRect(15, 200, 141, 51))
        self.seq_alig.setIcon(icon)
        self.seq_alig.setObjectName("seq_alig")
        self.browse = QtWidgets.QPushButton(self.tab, clicked = lambda: self.press_it())
        self.browse.setGeometry(QtCore.QRect(15, 60, 141, 51))
        self.browse.setObjectName("browse")
        self.clear = QtWidgets.QPushButton(self.tab, clicked = lambda: self.clear_it())
        self.clear.setGeometry(QtCore.QRect(15, 130, 141, 51))
        self.clear.setObjectName("clear")


        self.clear1 = QtWidgets.QPushButton(self.tab_2, clicked = lambda: self.clear_cl())
        self.clear1.setGeometry(QtCore.QRect(15, 130, 141, 51))
        self.clear1.setObjectName("clear1")
        
        # Genetic data
        self.te = QtWidgets.QTextEdit(self.tab)
        self.te.setGeometry(QtCore.QRect(210, 60, 1000, 351))
        self.te.setReadOnly(True)

        # Font type and size[
        font = QtGui.QFont()
        font.setFamily("Tibetan Machine Uni")
        font.setPointSize(10)
        self.te.setFont(font)
     
        self.te.setObjectName("te")
#       seq_alig Genetic data
        self.te1 = QtWidgets.QTextEdit(self.tab)
        self.te1.setGeometry(QtCore.QRect(210, 450, 1000, 351))
        self.te1.setObjectName("te1")
  

        # climatic data text box
        self.te2 = QtWidgets.QTextEdit(self.tab_2)
        self.te2.setGeometry(QtCore.QRect(210, 60, 1000, 351))
        self.te2.setObjectName("te2")
        self.te2.setReadOnly(True)

 
        # button to draw tree
        self.to_draw = QtWidgets.QPushButton(self.tab_2, clicked = lambda: self.openWindow())
        self.to_draw.setGeometry(QtCore.QRect(15, 200, 141, 51))
        self.to_draw.setObjectName("to_draw")
        self.browse_c = QtWidgets.QPushButton(self.tab_2, clicked = lambda: self.pressit())
        self.browse_c.setGeometry(QtCore.QRect(15, 60, 141, 51))
        self.browse_c.setObjectName("browse")

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
        self.tabWidget.setCurrentIndex(0)
        Main_Window.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
        QtCore.QMetaObject.connectSlotsByName(Main_Window)

        self.comboBox.currentIndexChanged.connect(self.show_frame)
        self.frame_3.setHidden(True)
    def show_frame(self, value):
        if value is not None:
            self.frame_3.setHidden(False)
        else:
            self.frame_3.setHidden(True)
    
    def press_it(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fileName, _ = QFileDialog.getOpenFileName(None,"QFileDialog.getOpenFileName()", "","All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            with open(fileName, "r") as f:
                content = f.read()
                self.te.setText(content)
                global sequence
                sequence = ''
                for char in content:
                    sequence += char.strip()
                self.child_window = QtWidgets.QMainWindow()
                self.ui = Ui_sec()
                self.ui.setupUi(self.child_window)
                self.child_window.setWindowModality(QtCore.Qt.NonModal)

        def color_background(letter):
            if letter == 'A':
                return 'background-color: yellow'
            elif letter == 'C':
                return 'background-color: blue'
            elif letter == 'G':
                return 'background-color: red'
            elif letter == 'T':
                return 'background-color: orange'
            else:
                return ''

        if 'sequence' in globals():
            formatted_sequence = ''.join(f'<span style="{color_background(l)}">{l}</span>' for l in sequence)
            self.te.setText(formatted_sequence)
            
  
    
    
    # press the button to get climatic data
    def pressit(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fileName, _ = QFileDialog.getOpenFileName(None,"QFileDialog.getOpenFileName()", "","All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            with open(fileName, "r") as c:
                lines = c.readlines()
                num_rows = len(lines)
                first_line = lines[1].split(",")
                num_columns = len(first_line)
                self.te2.clear()   #modified by Yannick
                cursor = QtGui.QTextCursor(self.te2.textCursor())
                cursor.insertTable(num_rows, num_columns)
                for line in lines:
                    line_split = line.split(",")
                    for value in line_split:
                        cursor.insertText(value)
                        cursor.movePosition(QtGui.QTextCursor.NextCell)
                self.child_window = QtWidgets.QMainWindow()
                self.ui = Ui_sec()
                self.ui.setupUi(self.child_window)
                self.child_window.setWindowModality(QtCore.Qt.NonModal)
                #self.child_window.show()
            


      # press the button to delet data
    def clear_it(self):
        self.te.clear()
    def clear_cl(self):
        self.te2.clear()


    def retranslateUi(self, Main_Window):
        _translate = QtCore.QCoreApplication.translate
        Main_Window.setWindowTitle(_translate("Main_Window", "aPhyloGeo_Plus_plus"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "GENETIC DATA"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "CLIMATIC DATA"))
        self.seq_alig.setToolTip(_translate("Main_Window", "<html><head/><body><p>sequence alignment</p></body></html>"))
        self.seq_alig.setText(_translate("Main_Window", "  Seq_Alig"))
        self.browse.setText(_translate("Main_Window", "Browse File"))
        self.browse_c.setText(_translate("Main_Window", "Browse \n Climatic Data"))
        self.to_draw.setText(_translate("Main_Window", "Draw \n Climatic Tree "))
        self.clear.setText(_translate("Main_Window", "Clear "))
        self.clear1.setText(_translate("Main_Window", "Clear "))
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

        self.comboBox.setItemText(0, _translate("MainWindow", "None"))
        self.comboBox.setItemText(1, _translate("MainWindow", "Bar Chart"))
        self.comboBox.setItemText(2, _translate("MainWindow", "Line Chart"))
        self.comboBox.setItemText(3, _translate("MainWindow", "Pie Chart"))
        self.comboBox.setItemText(4, _translate("MainWindow", "Area Chart"))
        self.comboBox.setItemText(5, _translate("MainWindow", "Scatter Chart"))
        self.comboBox_2.setItemText(0, _translate("MainWindow", "All"))
        self.comboBox_2.setItemText(1, _translate("MainWindow", "Temperature"))
        self.comboBox_2.setItemText(2, _translate("MainWindow", "Wind"))
        self.comboBox_2.setItemText(3, _translate("MainWindow", "Humidity"))
        self.comboBox_2.setItemText(4, _translate("MainWindow", "Altitude"))
        self.label_16.setText(_translate("MainWindow", "Climate \n"
"condition"))
        self.label_13.setText(_translate("MainWindow", "Statistics"))


if __name__ == "__main__":
    import sys
    #create_and_save_tree() #modified by Yannick  #removed by Yannick, added to cltree instead
    app = QtWidgets.QApplication(sys.argv)
    Main_Window = QtWidgets.QMainWindow()
    ui = Ui_Main_Window()
    ui.setupUi(Main_Window)
    Main_Window.show()
    Main_Window.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
    sys.exit(app.exec_())
