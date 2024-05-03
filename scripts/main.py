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
        if self.textEditPage1.toPlainText():
            self.pushButton_4.setEnabled(True)

    # to her]
    def __init__(self):
        super(UiMainWindow, self).__init__()
        uic.loadUi("Qt/main.ui", self)
        self.setupUi()

    def setupUi(self):
        self.setObjectName("MainWindow")
        self.fileBrowserButtonPage1.clicked.connect(self.press_it)
        self.sequenceAlignmentButtonPage1.clicked.connect(self.showSeqAlinFrame19)
        self.clearButtonPage1.clicked.connect(self.clearIt)
        self.statisticsButtonPage1.clicked.connect(self.showGenStatFrame4)
        self.geneticTreeButtonPage1.clicked.connect(self.showGenTreeFrame6)
        self.geneticDataButton.clicked.connect(self.showPage)
        self.geneticDataButton.clicked.connect(self.enableFrame)
        self.climaticDataButton.clicked.connect(self.showPage4)
        self.helpButton.clicked.connect(self.useWindow)
        self.resultsButton.clicked.connect(self.changeIconAndShowPage3)
        self.statisticsButtonPage2.clicked.connect(self.showGenStatFrame4)
        self.statisticsButtonPage2.clicked.connect(self.enableFrame4)
        self.fileBrowserButtonPage2.clicked.connect(self.press_it)
        self.fileBrowserButtonPage2.clicked.connect(self.changeIconAndShowPage)
        self.fileBrowserButtonPage2.clicked.connect(self.enableFrame)
        self.clearButtonPage2.clicked.connect(self.clearGenStat)
        self.clearButtonPage2.clicked.connect(self.resetCom2)
        self.backButtonPage2.clicked.connect(self.showPage)
        self.backButtonPage2.clicked.connect(self.enableFrame)
        self.fileBrowserButtonPage3.clicked.connect(self.press_it)
        self.fileBrowserButtonPage3.clicked.connect(self.changeIconAndShowPage)
        self.fileBrowserButtonPage3.clicked.connect(self.enableFrame)
        self.clearButtonPage3.clicked.connect(self.clearGenTree)
        self.backButtonPage3.clicked.connect(self.showPage)
        self.clearButtonPage4.clicked.connect(self.clearCl)
        self.climaticTreeButtonPage4.clicked.connect(self.openWindow)
        self.climaticTreeButtonPage4.clicked.connect(self.showClimTreeFrame13)
        self.climaticTreeButtonPage4.clicked.connect(self.enableFrame13)
        self.fileBrowserButtonPage4.clicked.connect(self.pressIt)
        self.statisticsButtonPage4.clicked.connect(self.showClimStatBar)
        self.statisticsButtonPage4.clicked.connect(self.showClimStatFrame10)
        self.statisticsButtonPage4.clicked.connect(self.enableFrame10)
        self.textEditPage4.textChanged.connect(self.onTextChangeClim)
        self.clearButtonPage5.clicked.connect(self.resetCom)
        self.clearButtonPage5.clicked.connect(self.clearClimStat)
        self.fileBrowserButtonPage5.clicked.connect(self.pressIt)
        self.fileBrowserButtonPage5.clicked.connect(self.showPage4)
        self.backButtonPage5.clicked.connect(self.showPage4)
        self.backButtonPage5.clicked.connect(self.enableFrame)
        self.clearButtonPage6.clicked.connect(self.clearClimTree)
        self.fileBrowserButtonPage6.clicked.connect(self.pressIt)
        self.fileBrowserButtonPage6.clicked.connect(self.showPage4)
        self.backButtonPage6.clicked.connect(self.showPage4)
        self.backButtonPage6.clicked.connect(self.enableFrame)
        self.settingsButtonPage7.clicked.connect(self.paramWin)
        self.submitButtonPage7.clicked.connect(self.showFilteredResults)
        self.statisticsButtonPage7.clicked.connect(self.showResultStatFrame16)
        self.statisticsButtonPage7.clicked.connect(self.enableFrame17)
        self.clearButtonPage7.clicked.connect(self.clearResult)
        self.settingsButtonPage8.clicked.connect(self.paramWin)
        self.settingsButtonPage8.clicked.connect(self.changeIconAndShowPage3)
        self.clearButtonPage8.clicked.connect(self.resetCom_4_5)
        self.clearButtonPage8.clicked.connect(self.clearResultStat)
        self.backButtonPage8.clicked.connect(self.showPage7)
        self.fileBrowserButtonPage9.clicked.connect(self.press_it)
        self.fileBrowserButtonPage9.clicked.connect(self.changeIconAndShowPage)
        self.fileBrowserButtonPage9.clicked.connect(self.enableFrame)
        self.clearButtonPage9.clicked.connect(self.clearSeq)
        self.backButtonPage9.clicked.connect(self.showPage)
        self.backButtonPage9.clicked.connect(self.enableFrame)

        self.translateUi()
        self.stackedWidget.setCurrentIndex(0)
        self.geneticDataButton.clicked.connect(self.showPage)
        self.climaticDataButton.clicked.connect(self.showPage4)
        self.resultsButton.clicked.connect(self.showPage7)
        self.sequenceAlignmentButtonPage1.clicked.connect(self.showSeqAlinFrame19)
        self.sequenceAlignmentButtonPage2.clicked.connect(self.showSeqAlinFrame19)
        self.sequenceAlignmentButtonPage3.clicked.connect(self.showSeqAlinFrame19)
        self.statisticsButtonPage1.clicked.connect(self.showGenStatFrame4)
        self.clearButtonPage3.clicked.connect(self.showGenStatFrame4)
        self.statisticsButtonPage9.clicked.connect(self.showGenStatFrame4)
        self.geneticTreeButtonPage1.clicked.connect(self.showGenTreeFrame6)
        self.geneticTreeButtonPage2.clicked.connect(self.showGenTreeFrame6)
        self.geneticTreeButtonPage9.clicked.connect(self.showGenTreeFrame6)
        self.climaticTreeButtonPage5.clicked.connect(self.showClimTreeFrame13)
        self.statisticsButtonPage4.clicked.connect(self.showClimStatFrame10)
        self.statisticButtonPage6.clicked.connect(self.showClimStatFrame10)
        self.statisticsButtonPage7.clicked.connect(self.showResultStatFrame16)

        buttons = [self.settingsButtonPage7, self.geneticDataButton, self.climaticDataButton, self.resultsButton,
                   self.settingsButtonPage8,
                   self.fileBrowserButtonPage1, self.clearButtonPage1, self.fileBrowserButtonPage2,
                   self.clearButtonPage2, self.fileBrowserButtonPage3, self.helpButton, self.sequenceAlignmentButtonPage1,
                   self.statisticsButtonPage1,
                   self.sequenceAlignmentButtonPage2, self.statisticsButtonPage2, self.geneticTreeButtonPage1,
                   self.sequenceAlignmentButtonPage3, self.clearButtonPage3, self.clearButtonPage3, self.geneticTreeButtonPage3, self.clearButtonPage4,
                   self.climaticTreeButtonPage4, self.fileBrowserButtonPage4, self.statisticsButtonPage4,
                   self.clearButtonPage5, self.climaticTreeButtonPage5, self.fileBrowserButtonPage5, self.statisticsButtonPage5, self.clearButtonPage6,
                   self.climaticTreeButtonPage6, self.fileBrowserButtonPage6, self.statisticButtonPage6,
                   self.submitButtonPage7, self.geneticTreeButtonPage2, self.statisticsButtonPage7, self.statisticsButtonPage7, self.submitButtonPage8,
                   self.statisticsButtonPage8, self.sequenceAlignmentButtonPage9,
                   self.statisticsButtonPage9, self.fileBrowserButtonPage9, self.clearButtonPage9, self.geneticTreeButtonPage9, self.clearButtonPage7,
                   self.clearButtonPage8, self.backButtonPage2,

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

        self.chartTypeComboBoxPage5.currentIndexChanged.connect(self.showFrame)
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
                self.textEditPage1.setText(content)
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
                self.textEditPage1.setText(sequence)

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
                self.textEditPage4.clear()
                cursor = QtGui.QTextCursor(self.textEditPage4.textCursor())
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
        self.textEditPage7.setText(str(content))

    # Enable_button():
    def onTextChanged(self):
        if self.textEditPage1.toPlainText() and self.textEditPage4.toPlainText():
            self.resultsButton.setEnabled(True)
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
        else:
            self.resultsButton.setEnabled(False)

    def onTextChangeGen(self):
        if self.textEditPage1.toPlainText():
            self.sequenceAlignmentButtonPage1.setEnabled(True)
            self.statisticsButtonPage1.setEnabled(True)
            self.geneticTreeButtonPage1.setEnabled(True)
            self.sequenceAlignmentButtonPage1.setIcon(QIcon(":inactive/sequence.svg"))
            self.statisticsButtonPage1.setIcon(QIcon(":inactive/statistics.svg"))
            self.geneticTreeButtonPage1.setIcon(QIcon(":inactive/tree.svg"))

        else:
            self.sequenceAlignmentButtonPage1.setEnabled(False)
            self.statisticsButtonPage1.setEnabled(False)
            self.geneticTreeButtonPage1.setEnabled(False)

    def onTextChangeClim(self):
        if self.textEditPage4.toPlainText():
            self.statisticsButtonPage4.setEnabled(True)
            self.climaticTreeButtonPage4.setEnabled(True)
            self.statisticsButtonPage4.setIcon(QIcon(":inactive/statistics.svg"))
            self.climaticTreeButtonPage4.setIcon(QIcon(":inactive/tree.svg"))

        else:
            self.statisticsButtonPage4.setEnabled(False)
            self.climaticTreeButtonPage4.setEnabled(False)

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
        if self.climaticDataButton.icon().isNull():
            self.climaticDataButton.setIcon(QIcon("icon2.png"))
        else:
            self.climaticDataButton.setIcon(QIcon(":active/climatic.svg"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
        self.showPage4()

    def changeIconAndShowPage(self):
        if self.geneticDataButton.icon().isNull():
            self.geneticDataButton.setIcon(QIcon("icon1.png"))
        else:
            self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
            self.geneticDataButton.setIcon(QIcon(":active/genetic.svg"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
        self.showPage()

    def changeIconAndShowPage3(self):
        if self.resultsButton.icon().isNull():
            self.resultsButton.setIcon(QIcon("icon3.png"))
        else:
            self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.resultsButton.setIcon(QIcon(":active/result.svg"))
        self.showPage7()

    # press the button to delete data
    def clearIt(self):
        self.textEd_4.clear()
        self.textEditPage1.clear()

    def clearSeq(self):
        self.textEd_4.clear()

    def clearGenStat(self):
        self.textEditPage2.clear()

    def clearGenTree(self):
        self.textEditpage3.clear()

    def clearCl(self):
        self.textEditPage4.clear()

    def clearClimStat(self):
        self.textBrowser_4.clear()

    def clearClimTree(self):
        self.textEditPage6.clear()

    def clearResult(self):
        self.textEditPage7.clear()
        # tableView

    def clearResultStat(self):
        self.textEditPage8.clear()
        # graphicsView_4

    # Set the combo box to its default value
    def resetCom2(self):
        self.speciesNamesList.setCurrentIndex(0)
        # reset_button.clicked.connect(self.resetCom2resetCom2)

    def resetCom(self):
        self.chartTypeComboBoxPage5.setCurrentIndex(0)
        self.climateConditionComboBoxPage5.setCurrentIndex(0)

    def resetCom_4_5(self):
        self.chartTypeComboBoxPage8.setCurrentIndex(0)
        self.conditionComboBoxPage8.setCurrentIndex(0)

    def translateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.helpButton.setToolTip(_translate("MainWindow", "How to use the application"))
        self.climaticDataButton.setToolTip(_translate("MainWindow", "Climatic Data"))
        self.resultsButton.setToolTip(_translate("MainWindow", "Results"))
        self.geneticDataButton.setToolTip(_translate("MainWindow", "Genetic Data"))
        self.geneticDataButton.setShortcut(_translate("MainWindow", "Down"))
        self.fileBrowserLabelPage1.setText(_translate("MainWindow", "File Browser"))
        self.sequenceAlignmentLabelPage1.setText(_translate("MainWindow", " Sequence \n"
                                                                          "Alignment"))
        self.clearLabelPage1.setText(_translate("MainWindow", "Clear"))
        self.statisticsLabelPage1.setText(_translate("MainWindow", "Statistics"))
        self.geneticTreeLabelPage1.setText(_translate("MainWindow", "Genetic Tree"))
        self.fileBrowserLabelPage2.setText(_translate("MainWindow", "File Browser"))
        self.sequenceAlignmentLabelPage2.setText(_translate("MainWindow", " Sequence \n"
                                                      "Alignment"))
        self.clearLabelPage2.setText(_translate("MainWindow", "Clear"))
        self.statisticsLabelPage2.setText(_translate("MainWindow", "Statistics"))
        self.geneticTreeLabelPage2.setText(_translate("MainWindow", "Genetic Tree"))
        self.backButtonPage2.setText(_translate("MainWindow", "back"))
        self.speciesNamesLabel.setText(_translate("MainWindow", "Species name"))
        self.speciesNamesList.setItemText(0, _translate("MainWindow", "All"))
        self.speciesNamesList.setItemText(1, _translate("MainWindow", "species1"))
        self.speciesNamesList.setItemText(2, _translate("MainWindow", "species2"))
        self.speciesNamesList.setItemText(3, _translate("MainWindow", "species3"))
        self.speciesNamesList.setItemText(4, _translate("MainWindow", "species4"))
        self.speciesTitleLabel.setText(_translate("MainWindow", "Statistics"))
        self.fileBrowserLabelPage3.setText(_translate("MainWindow", "File Browser"))
        self.sequenceAlignmentLabelPage3.setText(_translate("MainWindow", " Sequence \n"
                                                       "Alignment"))
        self.clearLabelPage3.setText(_translate("MainWindow", "Clear"))
        self.statisticsLabelPage3.setText(_translate("MainWindow", "Statistics"))
        self.geneticTreeLabelPage3.setText(_translate("MainWindow", "Genetic Tree"))
        self.geneticTreeLabel.setText(_translate("MainWindow", "Genetic Tree"))
        self.backButtonPage3.setText(_translate("MainWindow", "back"))
        self.fileBrowserLabelPage4.setText(_translate("MainWindow", "File Browser"))
        self.climaticTreeLabelPage4.setText(_translate("MainWindow", "Climatic Tree"))
        self.clearLabelPage4.setText(_translate("MainWindow", "Clear"))
        self.statisticsLabelPage4.setText(_translate("MainWindow", "Statistics"))
        self.fileBrowserLabelPage5.setText(_translate("MainWindow", "File Browser"))
        self.fileBrowserLabelPage5.setText(_translate("MainWindow", "Climatic Tree"))
        self.clearLabelPage5.setText(_translate("MainWindow", "Clear"))
        self.statisticsTitleLabelPage5.setText(_translate("MainWindow", "Statistics"))
        self.climateConditionComboBoxPage5.setItemText(0, _translate("MainWindow", "All"))
        self.climateConditionComboBoxPage5.setItemText(1, _translate("MainWindow", "Temperature"))
        self.climateConditionComboBoxPage5.setItemText(2, _translate("MainWindow", "Wind"))
        self.climateConditionComboBoxPage5.setItemText(3, _translate("MainWindow", "Humidity"))
        self.climateConditionComboBoxPage5.setItemText(4, _translate("MainWindow", "Altitude"))
        self.climateConditionLabelPage5.setText(_translate("MainWindow", "  Climate condition"))
        self.chartTypeComboBoxPage5.setItemText(0, _translate("MainWindow", "None"))
        self.chartTypeComboBoxPage5.setItemText(1, _translate("MainWindow", "Bar Chart"))
        self.chartTypeComboBoxPage5.setItemText(2, _translate("MainWindow", "Line Chart"))
        self.chartTypeComboBoxPage5.setItemText(3, _translate("MainWindow", "Pie Chart"))
        self.chartTypeComboBoxPage5.setItemText(4, _translate("MainWindow", "Area Chart"))
        self.chartTypeComboBoxPage5.setItemText(5, _translate("MainWindow", "Scatter Chart"))
        self.chartTypeLabelPage5.setText(_translate("MainWindow", "Chart Type "))
        self.backButtonPage5.setText(_translate("MainWindow", "back"))
        self.statisticsTitleLabelPage5.setText(_translate("MainWindow", "Statistics"))
        self.fileBrowserLabelPage6.setText(_translate("MainWindow", "File Browser"))
        self.climaticTreeLabelPage6.setText(_translate("MainWindow", "Climatic Tree"))
        self.clearLabelPage6.setText(_translate("MainWindow", "Clear"))
        self.statisticsLabelPage6.setText(_translate("MainWindow", "Statistics"))
        self.climaticTreeTitlePage6.setText(_translate("MainWindow", "Climatic Tree"))
        self.backButtonPage6.setText(_translate("MainWindow", "back"))
        self.settingsLabelPage7.setText(_translate("MainWindow", "Settings"))
        self.submitLabelPage7.setText(_translate("MainWindow", "Submit"))
        self.statisticsLabelPage7.setText(_translate("MainWindow", "Statistics"))
        self.clearLabelPage7.setText(_translate("MainWindow", "Clear"))
        self.resultTitlePage7.setText(_translate("MainWindow", "Result  "))
        self.settingsLabelPage8.setText(_translate("MainWindow", "Settings"))
        self.submitLabelPage8.setText(_translate("MainWindow", "Submit"))
        self.statisticsLabelPage8.setText(_translate("MainWindow", "Statistics"))
        self.clearLabelPage8.setText(_translate("MainWindow", "Clear"))
        self.backButtonPage8.setText(_translate("MainWindow", "back"))
        self.chartTypeComboBoxPage8.setItemText(0, _translate("MainWindow", "None"))
        self.chartTypeComboBoxPage8.setItemText(1, _translate("MainWindow", "Bar Chart"))
        self.chartTypeComboBoxPage8.setItemText(2, _translate("MainWindow", "Line Chart"))
        self.chartTypeComboBoxPage8.setItemText(3, _translate("MainWindow", "Pie Chart"))
        self.chartTypeComboBoxPage8.setItemText(4, _translate("MainWindow", "Area Chart"))
        self.chartTypeComboBoxPage8.setItemText(5, _translate("MainWindow", "Scatter Chart"))
        self.conditionLabelPage8.setText(_translate("MainWindow", "   condition"))
        self.chartTypeLabelPage8.setText(_translate("MainWindow", "Chart Type "))
        self.conditionComboBoxPage8.setItemText(0, _translate("MainWindow", "All"))
        self.conditionComboBoxPage8.setItemText(1, _translate("MainWindow", "Temperature"))
        self.conditionComboBoxPage8.setItemText(2, _translate("MainWindow", "Wind"))
        self.conditionComboBoxPage8.setItemText(3, _translate("MainWindow", "Humidity"))
        self.conditionComboBoxPage8.setItemText(4, _translate("MainWindow", "Altitude"))
        self.statisticsTitleLabelPage8.setText(_translate("MainWindow", "Statistics"))
        self.fileBrowserLabelPage9.setText(_translate("MainWindow", "File Browser"))
        self.sequenceAlignmentLabelPage9.setText(_translate("MainWindow", " Sequence \n"
                                                       "Alignment"))
        self.clearLabelPage9.setText(_translate("MainWindow", "Clear"))
        self.statisticsLabelPage9.setText(_translate("MainWindow", "Statistics"))
        self.geneticTreeLabelPage9.setText(_translate("MainWindow", "Genetic Tree"))
        self.sequenceAlignmentTitlePage9.setText(_translate("MainWindow", "Sequence Alignment"))
        self.backButtonPage9.setText(_translate("MainWindow", "back"))


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    window = UiMainWindow()
    mw = qtmodern.windows.ModernWindow(window)
    mw.show()
    window.show()
    sys.exit(app.exec_())
