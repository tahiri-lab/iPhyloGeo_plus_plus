import io
import re
from decimal import Decimal
import resources_rc
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QFileDialog
import sys
import aPhyloGeo.Alignement
import aPhyloGeo.aPhyloGeo
import folium
import matplotlib.pyplot as plt
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from help import UiHowToUse
from parameters import UiDialog
from PyQt5 import QtWidgets, uic
import qtmodern.styles
import qtmodern.windows
class UiMainWindow(QtWidgets.QMainWindow):

    # added code from her[

    def useWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = UiHowToUse()
        self.ui.setupUi(self.window)
        self.window.show()

    def paramWin(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = UiDialog()
        self.ui.setupUi(self.window)
        self.window.show()

    # open cl tree window
    def openWindow(self):
        from cltree import Ui_ct
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_ct()
        self.ui.setupUi(self.window)
        self.window.show()

    def enableButton(self):
        if self.textEdit_4.toPlainText():
            self.pushButton_4.setEnabled(True)

    # to her]
    def __init__(self):
        super(UiMainWindow, self).__init__()
        uic.loadUi("Qt/main.ui", self)
        self.setupUi()

    def setupUi(self):
        self.setObjectName("MainWindow")
        self.pushButton_6.clicked.connect(self.press_it)
        self.pushButton_12.clicked.connect(self.showSeqAlinFrame19)
        self.pushButton_7.clicked.connect(self.clearIt)
        self.pushButton_13.clicked.connect(self.showGenStatFrame4)
        self.pushButton_16.clicked.connect(self.showGenTreeFrame6)
        self.pushButton_2.clicked.connect(self.showPage)
        self.pushButton_2.clicked.connect(self.enableFrame)
        self.pushButton_3.clicked.connect(self.showPage4)
        self.pushButton_11.clicked.connect(self.useWindow)
        self.pushButton_4.clicked.connect(self.changeIconAndShowPage3)
        self.pushButton_15.clicked.connect(self.showGenStatFrame4)
        self.pushButton_15.clicked.connect(self.enableFrame4)
        self.pushButton_8.clicked.connect(self.press_it)
        self.pushButton_8.clicked.connect(self.changeIconAndShowPage)
        self.pushButton_8.clicked.connect(self.enableFrame)
        self.pushButton_9.clicked.connect(self.clearGenStat)
        self.pushButton_9.clicked.connect(self.resetCom2)
        self.back.clicked.connect(self.showPage)
        self.back.clicked.connect(self.enableFrame)
        self.pushButton_10.clicked.connect(self.press_it)
        self.pushButton_10.clicked.connect(self.changeIconAndShowPage)
        self.pushButton_10.clicked.connect(self.enableFrame)
        self.pushButton_19.clicked.connect(self.clearGenTree)
        self.back_2.clicked.connect(self.showPage)
        self.pushButton_21.clicked.connect(self.clearCl)
        self.pushButton_22.clicked.connect(self.openWindow)
        self.pushButton_22.clicked.connect(self.showClimTreeFrame13)
        self.pushButton_22.clicked.connect(self.enableFrame13)
        self.pushButton_23.clicked.connect(self.pressIt)
        self.pushButton_24.clicked.connect(self.showClimStatBar)
        self.pushButton_24.clicked.connect(self.showClimStatFrame10)
        self.pushButton_24.clicked.connect(self.enableFrame10)
        self.textBrowser_3.textChanged.connect(self.onTextChanged)
        self.textBrowser_3.textChanged.connect(self.onTextChangeClim)
        self.pushButton_25.clicked.connect(self.resetCom)
        self.pushButton_25.clicked.connect(self.clearClimStat)
        self.pushButton_27.clicked.connect(self.pressIt)
        self.pushButton_27.clicked.connect(self.showPage4)
        self.back_3.clicked.connect(self.showPage4)
        self.back_3.clicked.connect(self.enableFrame)
        self.pushButton_29.clicked.connect(self.clearClimTree)
        self.pushButton_31.clicked.connect(self.pressIt)
        self.pushButton_31.clicked.connect(self.showPage4)
        self.back_4.clicked.connect(self.showPage4)
        self.back_4.clicked.connect(self.enableFrame)
        self.pushButton.clicked.connect(self.paramWin)
        self.pushButton_33.clicked.connect(self.showFilteredResults)
        self.pushButton_35.clicked.connect(self.showResultStatFrame16)
        self.pushButton_35.clicked.connect(self.enableFrame17)
        self.pushButton_43.clicked.connect(self.clearResult)
        self.pushButton_5.clicked.connect(self.paramWin)
        self.pushButton_5.clicked.connect(self.changeIconAndShowPage3)
        self.pushButton_44.clicked.connect(self.resetCom_4_5)
        self.pushButton_44.clicked.connect(self.clearResultStat)
        self.back_5.clicked.connect(self.showPage7)
        self.pushButton_40.clicked.connect(self.press_it)
        self.pushButton_40.clicked.connect(self.changeIconAndShowPage)
        self.pushButton_40.clicked.connect(self.enableFrame)
        self.pushButton_41.clicked.connect(self.clearSeq)
        self.back_6.clicked.connect(self.showPage)
        self.back_6.clicked.connect(self.enableFrame)

        self.translateUi()
        self.stackedWidget.setCurrentIndex(0)
        self.pushButton_2.clicked.connect(self.showPage)
        self.pushButton_3.clicked.connect(self.showPage4)
        self.pushButton_4.clicked.connect(self.showPage7)
        self.pushButton_12.clicked.connect(self.showSeqAlinFrame19)
        self.pushButton_14.clicked.connect(self.showSeqAlinFrame19)
        self.pushButton_17.clicked.connect(self.showSeqAlinFrame19)
        self.pushButton_13.clicked.connect(self.showGenStatFrame4)
        self.pushButton_18.clicked.connect(self.showGenStatFrame4)
        self.pushButton_39.clicked.connect(self.showGenStatFrame4)
        self.pushButton_16.clicked.connect(self.showGenTreeFrame6)
        self.pushButton_34.clicked.connect(self.showGenTreeFrame6)
        self.pushButton_42.clicked.connect(self.showGenTreeFrame6)
        self.pushButton_26.clicked.connect(self.showClimTreeFrame13)
        self.pushButton_24.clicked.connect(self.showClimStatFrame10)
        self.pushButton_32.clicked.connect(self.showClimStatFrame10)
        self.pushButton_35.clicked.connect(self.showResultStatFrame16)

        buttons = [self.pushButton, self.pushButton_2, self.pushButton_3, self.pushButton_4, self.pushButton_5,
                   self.pushButton_6, self.pushButton_7, self.pushButton_8,
                   self.pushButton_9, self.pushButton_10, self.pushButton_11, self.pushButton_12, self.pushButton_13,
                   self.pushButton_14, self.pushButton_15, self.pushButton_16,
                   self.pushButton_17, self.pushButton_18, self.pushButton_19, self.pushButton_20, self.pushButton_21,
                   self.pushButton_22, self.pushButton_23, self.pushButton_24,
                   self.pushButton_25, self.pushButton_26, self.pushButton_27, self.pushButton_28, self.pushButton_29,
                   self.pushButton_30, self.pushButton_31, self.pushButton_32,
                   self.pushButton_33, self.pushButton_34, self.pushButton_35, self.pushButton_35, self.pushButton_36,
                   self.pushButton_37, self.pushButton_38,
                   self.pushButton_39, self.pushButton_40, self.pushButton_41, self.pushButton_42, self.pushButton_43,
                   self.pushButton_44, self.back,

                   ]

        # Définir le curseur et la feuille de style pour tous les boutons
        for button in buttons:
            button.setCursor(Qt.PointingHandCursor)
            button.setStyleSheet("""
            QPushButton:hover {
            background-color: grey;
            color: white;
            border: none;
            padding: 20px;
            border-radius: 10px;
            font-size: 13px;
            }
            """)

        QtCore.QMetaObject.connectSlotsByName(self)

        self.comboBox.currentIndexChanged.connect(self.showFrame)
        self.frame_11.setHidden(True)

    def press_it(self):
        '''
        Retrieve data from genetic file and show it in color
        '''
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fileName, _ = QFileDialog.getOpenFileName(None, "Select FASTA file", "", " (*.fasta);; (*.fasta)",
                                                  options=options)
        if fileName:
            aPhyloGeo.Alignement.userData_align.set_referenceGeneFile(fileName)
            with open(fileName, "r") as f:
                self.clearIt()
                content = f.read()
                self.textEdit_4.setText(content)
                sequence = ""
                for line in content.splitlines():
                    if not line.startswith('>'):
                        new_line = ''
                        for char in line:
                            if char == 'A':
                                new_line += f'<span style="background-color: yellow">{char}</span>'
                            elif char == 'C':
                                new_line += f'<span style="background-color: blue">{char}</span>'
                            elif char == 'G':
                                new_line += f'<span style="background-color: red">{char}</span>'
                            elif char == 'T':
                                new_line += f'<span style="background-color: orange">{char}</span>'
                        line = new_line
                    sequence += line + "<br>"
                self.textEdit_4.setText(sequence)

    def callSeqAlign(self):
        '''
        Initiate sequence alignment
        Return: genetic dictionary used for the final filter
        '''
        if self.textEd_4.toPlainText() == "":
            align_obj = aPhyloGeo.Alignement.AlignSequences()
            seq_al = align_obj.aligned
            obj = str(seq_al)
            self.textEd_4.setText(obj)
            self.genTree = aPhyloGeo.aPhyloGeo.createGenTree(align_obj)
        return self.genTree

    def retrieveDataNames(self, list):
        '''
        Retrieve data from a list, except for first element
        Args: 
         list (from which we get data)
        Return: 
         names_to_retrieve (Retrieved elements)
        '''
        names_to_retrieve = []
        for data in list:
            if data != list[0]:
                names_to_retrieve.append(data)
        return names_to_retrieve

    def populateMap(self, lat, long):
        '''
        Create folium map
        Args:
         lat (latitude), long (longitude)
        '''
        mean_lat = 0
        mean_long = 0
        for y in lat:
            mean_lat = mean_lat + Decimal(y)
        mean_lat = mean_lat / len(lat)
        for x in long:
            mean_long = mean_long + Decimal(x)
        mean_long = mean_long / len(long)
        m = folium.Map(location=[mean_lat, mean_long],
                       zoom_start=14,
                       tiles="OpenStreetMap")
        i = 0
        while i < len(lat):
            folium.Marker([Decimal(lat[i]), Decimal(long[i])]).add_to(m)
            i = i + 1
        # m.save('m.html')                  #folium map does not appear correctly on MAC with webbrowser.open(),
        # web browser.open('m.html')         #but it does not appear correctly on Linux with webview.show()
        # will have to be fixed
        data = io.BytesIO()
        m.save(data, close_file=False)
        self.webview = QWebEngineView()
        self.webview.setWindowTitle("Climate Map")
        self.webview.setHtml(data.getvalue().decode())
        self.webview.show()

    def pressIt(self):
        '''
        Retrieve data from climatic file and show it in a table
        If the last 2 columns of the data are 'LAT' and 'LONG',
        generate a folium map with these columns
        '''
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fileName, _ = QFileDialog.getOpenFileName(None, "Select CSV file", "", "Comma Separated Values (*.csv)",
                                                  options=options)

        if fileName:
            aPhyloGeo.aPhyloGeo.userData.set_fileName(fileName)
            with open(fileName, "r") as c:
                lines = c.readlines()
                num_rows = len(lines)
                first_line = lines[0].split(",")
                lat = []
                long = []
                self.species = []
                self.factors = [[], [], [], [], []]
                loc = False
                num_columns = len(first_line)
                if first_line[len(first_line) - 2] == 'LAT':
                    first_line_without_loc = first_line
                    first_line_without_loc.pop(len(first_line_without_loc) - 1)
                    first_line_without_loc.pop(len(first_line_without_loc) - 1)
                    clim_data_names = self.retrieveDataNames(first_line_without_loc)
                    aPhyloGeo.aPhyloGeo.userData.set_names(first_line_without_loc)
                    aPhyloGeo.Alignement.userData_align.set_names(first_line_without_loc)
                    loc = True
                else:
                    clim_data_names = self.retrieveDataNames(first_line)
                    aPhyloGeo.aPhyloGeo.userData.set_names(first_line)
                    aPhyloGeo.Alignement.userData_align.set_names(first_line)
                aPhyloGeo.aPhyloGeo.userData.set_dataNames(clim_data_names)
                aPhyloGeo.Alignement.userData_align.set_dataNames(clim_data_names)
                self.textBrowser_3.clear()
                cursor = QtGui.QTextCursor(self.textBrowser_3.textCursor())
                clim_data_table = cursor.insertTable(num_rows, num_columns)
                fmt = clim_data_table.format()
                fmt.setWidth(QtGui.QTextLength(QtGui.QTextLength.PercentageLength, 100))
                clim_data_table.setFormat(fmt)
                format = QtGui.QTextCharFormat()
                format.setForeground(QtGui.QColor('#006400'))
                i = 0
                for line in lines:
                    line_split = line.split(",")
                    if line != lines[0] and loc == True:
                        lat.append(line_split[len(line_split) - 2])
                        long.append(line_split[len(line_split) - 1])
                        self.species.append(line_split[0])
                        for j in [1, 2, 3, 4, 5]:
                            self.factors[i].append(line_split[j])
                        i += 1
                    if line != lines[0] and loc == False:
                        self.species.append(line_split[0])
                        for j in [1, 2, 3, 4, 5]:
                            self.factors[i].append(line_split[j])
                        i += 1
                    for value in line_split:
                        if line == lines[0]:
                            cursor.setCharFormat(format)
                        if re.search("^[0-9\\-]*\\.[0-9]*", value) is not None:
                            cursor.insertText(str(round(Decimal(value), 3)))
                            cursor.movePosition(QtGui.QTextCursor.NextCell)
                        else:
                            cursor.insertText(value)
                            cursor.movePosition(QtGui.QTextCursor.NextCell)
                if loc == True:
                    self.populateMap(lat, long)
                self.child_window = QtWidgets.QMainWindow()
                self.ui = UiHowToUse()
                self.ui.setupUi(self.child_window)
                self.child_window.setWindowModality(QtCore.Qt.NonModal)

    def showClimStatBar(self):
        # Condition should be added to show plot
        # according to user choice
        self.showClimStatBarAllFact()

    def showClimStatBarAllFact(self):
        '''
        Generate a bar graph that includes every factor for every species
        '''
        self.factors[0] = [float(v) for v in self.factors[0]]
        self.factors[1] = [float(v) for v in self.factors[1]]
        self.factors[2] = [float(v) for v in self.factors[2]]
        self.factors[3] = [float(v) for v in self.factors[3]]
        self.factors[4] = [float(v) for v in self.factors[4]]

        barWidth = 0.15
        fig = plt.subplots(figsize=(15, 10))
        br1 = np.arange(len(self.factors[0]))
        br2 = [x + barWidth for x in br1]
        br3 = [x + barWidth for x in br2]
        br4 = [x + barWidth for x in br3]
        br5 = [x + barWidth for x in br4]
        plt.bar(br1, self.factors[0], color='black', width=barWidth,
                edgecolor='grey', label=self.species[0])
        plt.bar(br2, self.factors[1], color='red', width=barWidth,
                edgecolor='grey', label=self.species[1])
        plt.bar(br3, self.factors[2], color='green', width=barWidth,
                edgecolor='grey', label=self.species[2])
        plt.bar(br4, self.factors[3], color='blue', width=barWidth,
                edgecolor='grey', label=self.species[3])
        plt.bar(br5, self.factors[4], color='cyan', width=barWidth,
                edgecolor='grey', label=self.species[4])
        plt.xlabel('COVID-19 Variant', fontweight='bold')
        plt.ylabel("Climatic data", fontweight='bold')
        plt.xticks([r + barWidth for r in range(len(self.factors[0]))],
                   aPhyloGeo.aPhyloGeo.userData.get_dataNames())
        plt.yticks(np.arange(0, 30))
        plt.title("Distribution of climatic variables for each COVID Variant")
        plt.legend()
        plt.show()

    def showFrame(self, value):
        if value is not None:
            self.frame_11.setHidden(False)
        else:
            self.frame_11.setHidden(True)

    def showPage(self):
        self.stackedWidget.setCurrentIndex(0)

    def showGenStatFrame4(self):
        self.stackedWidget.setCurrentIndex(1)

    def showGenTreeFrame6(self):
        self.stackedWidget.setCurrentIndex(2)

    def showPage4(self):
        self.stackedWidget.setCurrentIndex(3)

    def showClimStatFrame10(self):
        self.stackedWidget.setCurrentIndex(4)

    def showClimTreeFrame13(self):
        self.stackedWidget.setCurrentIndex(5)

    def showPage7(self):
        self.stackedWidget.setCurrentIndex(6)

    def showResultStatFrame16(self):
        self.stackedWidget.setCurrentIndex(7)

    def showSeqAlinFrame19(self):
        self.geneticTreeDict = self.callSeqAlign()
        self.stackedWidget.setCurrentIndex(8)

    def showFilteredResults(self):
        '''
        Show the results filtered with a metric threshold provided by user
        '''
        aPhyloGeo.aPhyloGeo.filterResults(
            aPhyloGeo.aPhyloGeo.climaticPipeline(aPhyloGeo.aPhyloGeo.userData.get_fileName(),
                                                 aPhyloGeo.aPhyloGeo.userData.get_names()),
            self.geneticTreeDict)
        with open("output.csv", "r") as f:
            content = f.read()
        self.textBrowser_7.setText(str(content))

    # Enable_button():
    def onTextChanged(self):
        if self.textEdit_4.toPlainText() and self.textBrowser_3.toPlainText():
            self.pushButton_4.setEnabled(True)
            self.pushButton_4.setIcon(QIcon(":inactive/result.svg"))
        else:
            self.pushButton_4.setEnabled(False)

    def onTextChangeGen(self):
        if self.textEdit_4.toPlainText():
            self.pushButton_12.setEnabled(True)
            self.pushButton_13.setEnabled(True)
            self.pushButton_16.setEnabled(True)
            self.pushButton_12.setIcon(QIcon(":inactive/sequence.svg"))
            self.pushButton_13.setIcon(QIcon(":inactive/statistics.svg"))
            self.pushButton_16.setIcon(QIcon(":inactive/tree.svg"))

        else:
            self.pushButton_12.setEnabled(False)
            self.pushButton_13.setEnabled(False)
            self.pushButton_16.setEnabled(False)

    def onTextChangeClim(self):
        if self.textBrowser_3.toPlainText():
            self.pushButton_24.setEnabled(True)
            self.pushButton_22.setEnabled(True)
            self.pushButton_24.setIcon(QIcon(":inactive/statistics.svg"))
            self.pushButton_22.setIcon(QIcon(":inactive/tree.svg"))

        else:
            self.pushButton_24.setEnabled(False)
            self.pushButton_22.setEnabled(False)

    # enable the frame when the push button is clicked
    def enableFrame(self):
        self.frame_2.setEnabled(True)
        self.frame.setEnabled(True)

    def enableFrame8(self):
        self.frame_8.setEnabled(False)
        self.frame_7.setEnabled(False)

    def enableFrame4(self):
        self.frame_4.setEnabled(True)
        self.showGenStatFrame4()

    def enableFrame10(self):
        self.frame_10.setEnabled(True)
        self.showGenStatFrame4()

    def enableFrame13(self):
        self.frame_13.setEnabled(True)
        # self.showGenStatFrame4()

    def enableFrame17(self):
        self.frame_17.setEnabled(True)

    def changeIconAndShowPage2(self):
        if self.pushButton_3.icon().isNull():
            self.pushButton_3.setIcon(QIcon("icon2.png"))
        else:
            self.pushButton_3.setIcon(QIcon(":active/climatic.svg"))
            self.pushButton_2.setIcon(QIcon(":inactive/genetic.svg"))
            self.pushButton_4.setIcon(QIcon(":inactive/result.svg"))
        self.showPage4()

    def changeIconAndShowPage(self):
        if self.pushButton_2.icon().isNull():
            self.pushButton_2.setIcon(QIcon("icon1.png"))
        else:
            self.pushButton_3.setIcon(QIcon(":inactive/climatic.svg"))
            self.pushButton_2.setIcon(QIcon(":active/genetic.svg"))
            self.pushButton_4.setIcon(QIcon(":inactive/result.svg"))
        self.showPage()

    def changeIconAndShowPage3(self):
        if self.pushButton_4.icon().isNull():
            self.pushButton_4.setIcon(QIcon("icon3.png"))
        else:
            self.pushButton_3.setIcon(QIcon(":inactive/climatic.svg"))
            self.pushButton_2.setIcon(QIcon(":inactive/genetic.svg"))
            self.pushButton_4.setIcon(QIcon(":active/result.svg"))
        self.showPage7()

    # press the button to delete data
    def clearIt(self):
        self.textEd_4.clear()
        self.textEdit_4.clear()

    def clearSeq(self):
        self.textEd_4.clear()

    def clearGenStat(self):
        self.textBrowser.clear()

    def clearGenTree(self):
        self.textBrowser_2.clear()

    def clearCl(self):
        self.textBrowser_3.clear()

    def clearClimStat(self):
        self.textBrowser_4.clear()

    def clearClimTree(self):
        self.textBrowser_5.clear()

    def clearResult(self):
        self.textBrowser_7.clear()
        # tableView

    def clearResultStat(self):
        self.textBrowser_6.clear()
        # graphicsView_4

    # Set the combo box to its default value
    def resetCom2(self):
        self.comboBox_2.setCurrentIndex(0)
        # reset_button.clicked.connect(self.resetCom2resetCom2)

    def resetCom(self):
        self.comboBox.setCurrentIndex(0)
        self.comboBox_3.setCurrentIndex(0)

    def resetCom_4_5(self):
        self.comboBox_4.setCurrentIndex(0)
        self.comboBox_5.setCurrentIndex(0)

    def translateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.pushButton_11.setToolTip(_translate("MainWindow", "How to use the application"))
        self.pushButton_3.setToolTip(_translate("MainWindow", "Climatic Data"))
        self.pushButton_4.setToolTip(_translate("MainWindow", "Results"))
        self.pushButton_2.setToolTip(_translate("MainWindow", "Genetic Data"))
        self.pushButton_2.setShortcut(_translate("MainWindow", "Down"))
        self.label_2.setText(_translate("MainWindow", "File Browser"))
        self.label_5.setText(_translate("MainWindow", " Sequence \n"
                                                      "Alignment"))
        self.label_3.setText(_translate("MainWindow", "Clear"))
        self.label_7.setText(_translate("MainWindow", "Statistics"))
        self.label_11.setText(_translate("MainWindow", "Genetic Tree"))
        self.label_4.setText(_translate("MainWindow", "File Browser"))
        self.label_8.setText(_translate("MainWindow", " Sequence \n"
                                                      "Alignment"))
        self.label_9.setText(_translate("MainWindow", "Clear"))
        self.label_10.setText(_translate("MainWindow", "Statistics"))
        self.label_37.setText(_translate("MainWindow", "Genetic Tree"))
        self.back.setText(_translate("MainWindow", "back"))
        self.label.setText(_translate("MainWindow", "Species name"))
        self.comboBox_2.setItemText(0, _translate("MainWindow", "All"))
        self.comboBox_2.setItemText(1, _translate("MainWindow", "species1"))
        self.comboBox_2.setItemText(2, _translate("MainWindow", "species2"))
        self.comboBox_2.setItemText(3, _translate("MainWindow", "species3"))
        self.comboBox_2.setItemText(4, _translate("MainWindow", "species4"))
        self.label_38.setText(_translate("MainWindow", "Statistics"))
        self.label_12.setText(_translate("MainWindow", "File Browser"))
        self.label_13.setText(_translate("MainWindow", " Sequence \n"
                                                       "Alignment"))
        self.label_14.setText(_translate("MainWindow", "Clear"))
        self.label_15.setText(_translate("MainWindow", "Statistics"))
        self.label_16.setText(_translate("MainWindow", "Genetic Tree"))
        self.label_17.setText(_translate("MainWindow", "Genetic Tree"))
        self.back_2.setText(_translate("MainWindow", "back"))
        self.label_18.setText(_translate("MainWindow", "File Browser"))
        self.label_19.setText(_translate("MainWindow", "Climatic Tree"))
        self.label_20.setText(_translate("MainWindow", "Clear"))
        self.label_21.setText(_translate("MainWindow", "Statistics"))
        self.label_22.setText(_translate("MainWindow", "File Browser"))
        self.label_23.setText(_translate("MainWindow", "Climatic Tree"))
        self.label_24.setText(_translate("MainWindow", "Clear"))
        self.label_25.setText(_translate("MainWindow", "Statistics"))
        self.comboBox_3.setItemText(0, _translate("MainWindow", "All"))
        self.comboBox_3.setItemText(1, _translate("MainWindow", "Temperature"))
        self.comboBox_3.setItemText(2, _translate("MainWindow", "Wind"))
        self.comboBox_3.setItemText(3, _translate("MainWindow", "Humidity"))
        self.comboBox_3.setItemText(4, _translate("MainWindow", "Altitude"))
        self.label_26.setText(_translate("MainWindow", "  Climate condition"))
        self.comboBox.setItemText(0, _translate("MainWindow", "None"))
        self.comboBox.setItemText(1, _translate("MainWindow", "Bar Chart"))
        self.comboBox.setItemText(2, _translate("MainWindow", "Line Chart"))
        self.comboBox.setItemText(3, _translate("MainWindow", "Pie Chart"))
        self.comboBox.setItemText(4, _translate("MainWindow", "Area Chart"))
        self.comboBox.setItemText(5, _translate("MainWindow", "Scatter Chart"))
        self.label_27.setText(_translate("MainWindow", "Chart Type "))
        self.back_3.setText(_translate("MainWindow", "back"))
        self.label_39.setText(_translate("MainWindow", "Statistics"))
        self.label_28.setText(_translate("MainWindow", "File Browser"))
        self.label_29.setText(_translate("MainWindow", "Climatic Tree"))
        self.label_30.setText(_translate("MainWindow", "Clear"))
        self.label_31.setText(_translate("MainWindow", "Statistics"))
        self.label_32.setText(_translate("MainWindow", "Climatic Tree"))
        self.back_4.setText(_translate("MainWindow", "back"))
        self.label_33.setText(_translate("MainWindow", "Settings"))
        self.label_34.setText(_translate("MainWindow", "Submit"))
        self.label_40.setText(_translate("MainWindow", "Statistics"))
        self.label_52.setText(_translate("MainWindow", "Clear"))
        self.label_36.setText(_translate("MainWindow", "Result  "))
        self.label_41.setText(_translate("MainWindow", "Settings"))
        self.label_42.setText(_translate("MainWindow", "Submit"))
        self.label_43.setText(_translate("MainWindow", "Statistics"))
        self.label_53.setText(_translate("MainWindow", "Clear"))
        self.back_5.setText(_translate("MainWindow", "back"))
        self.comboBox_4.setItemText(0, _translate("MainWindow", "None"))
        self.comboBox_4.setItemText(1, _translate("MainWindow", "Bar Chart"))
        self.comboBox_4.setItemText(2, _translate("MainWindow", "Line Chart"))
        self.comboBox_4.setItemText(3, _translate("MainWindow", "Pie Chart"))
        self.comboBox_4.setItemText(4, _translate("MainWindow", "Area Chart"))
        self.comboBox_4.setItemText(5, _translate("MainWindow", "Scatter Chart"))
        self.label_44.setText(_translate("MainWindow", "   condition"))
        self.label_45.setText(_translate("MainWindow", "Chart Type "))
        self.comboBox_5.setItemText(0, _translate("MainWindow", "All"))
        self.comboBox_5.setItemText(1, _translate("MainWindow", "Temperature"))
        self.comboBox_5.setItemText(2, _translate("MainWindow", "Wind"))
        self.comboBox_5.setItemText(3, _translate("MainWindow", "Humidity"))
        self.comboBox_5.setItemText(4, _translate("MainWindow", "Altitude"))
        self.label_51.setText(_translate("MainWindow", "Statistics"))
        self.label_35.setText(_translate("MainWindow", "File Browser"))
        self.label_46.setText(_translate("MainWindow", " Sequence \n"
                                                       "Alignment"))
        self.label_47.setText(_translate("MainWindow", "Clear"))
        self.label_48.setText(_translate("MainWindow", "Statistics"))
        self.label_49.setText(_translate("MainWindow", "Genetic Tree"))
        self.label_50.setText(_translate("MainWindow", "Sequence Alignment"))
        self.back_6.setText(_translate("MainWindow", "back"))


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    window = UiMainWindow()
    mw = qtmodern.windows.ModernWindow(window)
    mw.show()
    window.show()
    sys.exit(app.exec_())
