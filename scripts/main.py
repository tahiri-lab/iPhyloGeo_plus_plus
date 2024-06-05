import io
import os
import re
import sys
from decimal import Decimal
import folium
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qtmodern.styles
import qtmodern.styles
import qtmodern.windows
import qtmodern.windows
import yaml
from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon, QColor
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QFileDialog, QGraphicsDropShadowEffect
from aphylogeo import utils
from aphylogeo.alignement import AlignSequences
from aphylogeo.genetic_trees import GeneticTrees
from aphylogeo.params import Params
import resources_rc
from help import UiHowToUse
from parameters import UiDialog


class MyDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(MyDumper, self).increase_indent(flow, False)

    def represent_list(self, data):
        return self.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.add_representer(list, MyDumper.represent_list, Dumper=MyDumper)


def update_yaml_param(params, file_path, property_name, new_value):
    """
    Updates a specified property within a YAML file with a new value.

    Args:
        file_path (str): The path to the YAML file.
        property_name (str): The name of the property to modify (e.g., 'file_name').
        new_value: The new value to set for the property (can be any valid YAML type).
    """
    if isinstance(new_value, list):
        new_value = [element.strip() for element in new_value]
    params.update_from_dict({property_name: new_value})

    # 1. Load existing YAML data
    with open(file_path, "r") as yaml_file:
        data = yaml.safe_load(yaml_file)  # Use safe_load for security

    # 2. Update the specified property
    if property_name in data:
        data[property_name] = new_value
    else:
        print(f"Warning: Property '{property_name}' not found in '{file_path}'.")

    # 3. Write the updated data back to the file
    with open(file_path, "w") as yaml_file:
        yaml.dump(data, yaml_file, default_flow_style=None, Dumper=MyDumper, sort_keys=False)


Params.load_from_file("params.yaml")


class UiMainWindow(QtWidgets.QMainWindow):

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
    def openClimTree(self):
        from cltree import Ui_ct
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_ct()
        self.ui.setupUi(self.window)
        self.window.show()

        self.stackedWidget.setCurrentIndex(2)
        self.tabWidget2.setCurrentIndex(3)

    # to her]
    def __init__(self):
        super(UiMainWindow, self).__init__()
        uic.loadUi("Qt/main.ui", self)
        self.setupUi()

    def toggleDarkMode(self):
        self.isDarkMode = not self.isDarkMode
        if self.isDarkMode:
            qtmodern.styles.dark(app)
            self.darkModeButton.setIcon(QIcon(":other/light.png"))  # Set the 'light' icon for dark mode
        else:
            qtmodern.styles.light(app)
            self.darkModeButton.setIcon(QIcon(":other/dark.png"))  # Set the 'dark' icon

    def setupUi(self):
        self.setObjectName("MainWindow")

        self.homeButton.clicked.connect(self.showHomePage)
        self.geneticDataButton.clicked.connect(self.showGenDatPage)
        self.climaticDataButton.clicked.connect(self.showClimDatPage)
        self.helpButton.clicked.connect(self.useWindow)
        self.darkModeButton.clicked.connect(self.toggleDarkMode)
        self.darkModeButton.setCursor(Qt.PointingHandCursor)
        self.isDarkMode = False  # Keep track of the state

        self.fileBrowserButtonPage1.clicked.connect(self.pressItFasta)
        self.sequenceAlignmentButtonPage1.clicked.connect(self.SeqAlign)
        self.clearButtonPage1.clicked.connect(self.clearGen)
        self.statisticsButtonPage1.clicked.connect(self.showClimStatBarAllFact)

        self.clearButtonPage2.clicked.connect(self.clearClim)
        self.climaticTreeButtonPage2.clicked.connect(self.openClimTree)
        self.fileBrowserButtonPage2.clicked.connect(self.pressItCSV)
        self.statisticsButtonPage2.clicked.connect(self.showClimStatBarAllFact)
        self.resultsButtonPage2.clicked.connect(self.showResultsPage)

        self.settingsButtonPage4.clicked.connect(self.paramWin)
        self.submitButtonPage3.clicked.connect(self.showFilteredResults)
        self.clearButtonPage4.clicked.connect(self.clearResult)
        self.statisticsButtonPage4.clicked.connect(self.showResultsStatsPage)

        self.clearButtonPage4.clicked.connect(self.clearResultStat)

        self.stackedWidget.setCurrentIndex(0)

        buttons = [self.geneticDataButton, self.climaticDataButton, self.helpButton, self.homeButton]

        buttons_Vertical = [self.fileBrowserButtonPage1, self.sequenceAlignmentButtonPage1,
                            self.clearButtonPage1, self.statisticsButtonPage1, self.geneticTreeButtonPage1,
                            self.fileBrowserButtonPage2, self.clearButtonPage2, self.resultsButtonPage2,
                            self.climaticTreeButtonPage2, self.statisticsButtonPage2,
                            self.settingsButtonPage3,
                            self.settingsButtonPage4,
                            self.submitButtonPage3,
                            self.statisticsButtonPage3,
                            self.submitButtonPage4,
                            self.statisticsButtonPage4,
                            self.clearButtonPage3,
                            self.clearButtonPage4]

        # Définir le curseur et la feuille de style pour tous les boutons
        for button in buttons:
            button.setCursor(Qt.PointingHandCursor)
            button.setStyleSheet("""
                QPushButton { 
                    padding: 10px 20px;
                    font-weight: bold;
                    background-color: #DEDDDA;
                    border-radius: 14px;
                    transition: background-color 0.3s ease; /* Add transition */
                }
                QPushButton:hover {
                    background-color: #B7B7B6; 
                }
                QPushButton:pressed {
                    background-color: #DEDDDA;
                }
            """)
            # Créer l'effet d'ombre
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)  # Ajuster le flou de l'ombre
            shadow_effect.setColor(QColor(0, 0, 0, 140))  # Couleur de l'ombre (noir avec transparence)
            shadow_effect.setOffset(3, 3)  # Décalage de l'ombre

            # Appliquer l'effet d'ombre au bouton
            button.setGraphicsEffect(shadow_effect)

        # Définir le curseur et la feuille de style pour tous les boutons
        for buttonV in buttons_Vertical:
            buttonV.setCursor(Qt.PointingHandCursor)
            buttonV.setStyleSheet("""
                QPushButton {
                    border-radius: 14px;
                    background-color: #EEEEEE;
                    padding: 10px 20px;
                    font-weight: bold;
                    transition: background-color 0.3s ease; /* Add transition */
                }
                QPushButton:hover {
                    background-color: #D7D7D7; 
                }
                QPushButton:pressed {
                    background-color: #EEEEEE;
                }
            """)
            # Créer l'effet d'ombre
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)  # Ajuster le flou de l'ombre
            shadow_effect.setColor(QColor(0, 0, 0, 110))  # Couleur de l'ombre (noir avec transparence)
            shadow_effect.setOffset(3, 3)  # Décalage de l'ombre

            # Appliquer l'effet d'ombre au bouton
            buttonV.setGraphicsEffect(shadow_effect)

        self.darkModeButton.setCursor(Qt.PointingHandCursor)
        self.darkModeButton.setStyleSheet("""
                QPushButton { 
                    padding: 10px 20px;
                    font-weight: bold;
                    background-color: #DEDDDA;
                    border-radius: 14px;
                    transition: background-color 0.3s ease; /* Add transition */
                }
                QPushButton:hover {
                    background-color: #B7B7B6; 
                }
                QPushButton:pressed {
                    background-color: #DEDDDA;
                }
            """)
        # Créer l'effet d'ombre
        shadow_effect = QGraphicsDropShadowEffect()
        shadow_effect.setBlurRadius(10)  # Ajuster le flou de l'ombre
        shadow_effect.setColor(QColor(0, 0, 0, 140))  # Couleur de l'ombre (noir avec transparence)
        shadow_effect.setOffset(3, 3)  # Décalage de l'ombre

        # Appliquer l'effet d'ombre au bouton
        self.darkModeButton.setGraphicsEffect(shadow_effect)
        QtCore.QMetaObject.connectSlotsByName(self)

    def pressItFasta(self):
        '''
        Retrieve data from genetic file and show it in color
        '''
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fullFileName, _ = QFileDialog.getOpenFileName(None, "Select FASTA file", "../datasets",
                                                      " (*.fasta);; (*.fasta)",
                                                      options=options)
        if fullFileName:
            update_yaml_param(Params, "params.yaml", "reference_gene_file", os.path.basename(fullFileName))
            update_yaml_param(Params, "params.yaml", "reference_gene_dir", os.path.dirname(fullFileName))
            with open(fullFileName, "r") as f:
                self.clearGen()
                content = f.read()
                self.textEditFasta.setText(content)
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
                self.textEditFasta.setText(sequence)
                self.tabWidget.setCurrentIndex(1)
                self.sequenceAlignmentButtonPage1.setEnabled(True)
                self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
                self.sequenceAlignmentButtonPage1.setIcon(QIcon(":inactive/sequence.svg"))

    def SeqAlign(self):
        self.geneticTreeDict = self.callSeqAlign()

    def callSeqAlign(self):
        import time
        from PyQt5 import QtCore
        def update_progress(loading_screen, step):
            if step < loading_screen.checkListWidget.count():
                item = loading_screen.checkListWidget.item(step)
                item.setCheckState(QtCore.Qt.Checked)
                progress_value = int((step + 1) * (100 / loading_screen.checkListWidget.count()))
                loading_screen.progressBar.setValue(progress_value)
            else:
                loading_screen.progressBar.setValue(100)
            time.sleep(0.8)

        '''
        Initiate sequence alignment
        Return: genetic dictionary used for the final filter
        '''
        loading_screen = uic.loadUi("Qt/loading.ui")
        loading_screen.setWindowModality(QtCore.Qt.ApplicationModal)

        loading_screen.show()

        QtWidgets.QApplication.processEvents()
        update_progress(loading_screen, 0)
        QtWidgets.QApplication.processEvents()
        try:
            # Step 1: Load sequences
            sequenceFile = utils.loadSequenceFile(Params.reference_gene_filepath)
            update_progress(loading_screen, 1)
            QtWidgets.QApplication.processEvents()

            # Step 2: Align sequences
            align_sequence = AlignSequences(sequenceFile)
            alignements = align_sequence.align()
            update_progress(loading_screen, 2)
            QtWidgets.QApplication.processEvents()

            geneticTrees = utils.geneticPipeline(alignements.msa)
            trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
            update_progress(loading_screen, 3)
            QtWidgets.QApplication.processEvents()

            # Step 3: Preparing results
            obj = str(alignements.to_dict())
            self.textEditSeqAlign.setText(obj)
            self.tabWidget.setCurrentIndex(2)
            self.statisticsButtonPage1.setEnabled(True)
            self.statisticsButtonPage1.setIcon(QIcon(":inactive/statistics.svg"))
            self.resultsButtonPage2.setEnabled(True)
            update_progress(loading_screen, 4)
            QtWidgets.QApplication.processEvents()

            # Step 4: Save results
            alignements.save_to_json(f"./results/aligned_{Params.reference_gene_file}.json")
            trees.save_trees_to_json("./results/geneticTrees.json")
            update_progress(loading_screen, 5)
            QtWidgets.QApplication.processEvents()
            time.sleep(0.8)
        finally:
            loading_screen.close()

        return geneticTrees

    def load_sequences(self):
        # Simulate loading sequences
        self.sequenceFile = utils.loadSequenceFile(Params.reference_gene_filepath)
        print("Sequences loaded")

    def align_sequences(self):
        # Simulate aligning sequences
        align_sequence = AlignSequences(self.sequenceFile)
        self.alignements = align_sequence.align()
        print("Sequences aligned")

    def construct_genetic_trees(self):
        # Simulate constructing genetic trees
        geneticTrees = utils.geneticPipeline(self.alignements.msa)
        self.trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
        self.geneticTrees = geneticTrees
        print("Genetic trees constructed")

    def prepare_results(self):
        # Simulate preparing results
        obj = str(self.alignements.to_dict())
        self.textEd_4.setText(obj)
        self.resultsButton.setEnabled(True)
        print("Results prepared")

    def save_results(self):
        # Simulate saving results
        self.alignements.save_to_json(f"./results/aligned_{Params.reference_gene_file}.json")
        self.trees.save_trees_to_json("./results/geneticTrees.json")
        print("Results saved")

    def finalize_alignment(self):
        return self.geneticTrees

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

    def pressItCSV(self):
        '''
        Retrieve data from climatic file and show it in a table
        If the last 2 columns of the data are 'LAT' and 'LONG',
        generate a folium map with these columns
        '''
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fullFilePath, _ = QFileDialog.getOpenFileName(None, "Select CSV file", "../datasets",
                                                      "Comma Separated Values (*.csv)",
                                                      options=options)

        if fullFilePath:
            update_yaml_param(Params, "params.yaml", "file_name", fullFilePath)

            with open(fullFilePath, "r") as c:
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
                    update_yaml_param(Params, "params.yaml", "names", first_line_without_loc)
                    loc = True
                else:
                    clim_data_names = self.retrieveDataNames(first_line)
                    update_yaml_param(Params, "params.yaml", "names", first_line)
                update_yaml_param(Params, "params.yaml", "data_names", clim_data_names)
                self.textEditClimData.clear()
                cursor = QtGui.QTextCursor(self.textEditClimData.textCursor())
                clim_data_table = cursor.insertTable(num_rows, num_columns)
                fmt = clim_data_table.format()
                fmt.setWidth(QtGui.QTextLength(QtGui.QTextLength.PercentageLength, 98))
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
        self.tabWidget2.setCurrentIndex(1)
        self.resultsButtonPage2.setEnabled(True)
        self.resultsButtonPage2.setIcon(QIcon(":inactive/result.svg"))
        self.climaticTreeButtonPage2.setEnabled(True)
        self.climaticTreeButtonPage2.setIcon(QIcon(":inactive/climatic.svg"))
        self.statisticsButtonPage2.setEnabled(True)
        self.statisticsButtonPage2.setIcon(QIcon(":inactive/statistics.svg"))

    def showClimStatBarAllFact(self):
        '''
        Generate a bar graph that includes every factor for every species
        '''
        try:
            print(self.factors)
        except AttributeError:
            self.textEditGenStats.setText("Valeurs non existantes, veuillez choisir un fichier !")
            self.tabWidget.setCurrentIndex(3)

        else:
            self.tabWidget.setCurrentIndex(3)
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
                       Params.data_names)
            plt.yticks(np.arange(0, 30))
            plt.title("Distribution of climatic variables for each COVID Variant")
            plt.legend()
            plt.show()

    def showHomePage(self):
        self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
        self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
        self.stackedWidget.setCurrentIndex(0)

    def showGenDatPage(self):
        self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
        self.geneticDataButton.setIcon(QIcon(":active/genetic.svg"))
        self.stackedWidget.setCurrentIndex(1)
        self.tabWidget.setCurrentIndex(0)

    def showClimDatPage(self):
        self.climaticDataButton.setIcon(QIcon(":active/climatic.svg"))
        self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
        self.stackedWidget.setCurrentIndex(2)
        self.tabWidget2.setCurrentIndex(0)

    def showResultsPage(self):
        self.stackedWidget.setCurrentIndex(3)

    def showResultsStatsPage(self):
        self.stackedWidget.setCurrentIndex(4)

    def showFilteredResults(self):
        '''
        Show the results filtered with a metric threshold provided by user
        '''
        try:
            df = pd.read_csv(Params.file_name)
            climaticTrees = utils.climaticPipeline(df)
            utils.filterResults(climaticTrees, self.geneticTreeDict, df)
            df_results = pd.read_csv("./results/output.csv")

            # Convert to HTML table with basic styling
            html_table = df_results.to_html(index=False, border=1, classes="dataframe")  # Basic styling

            # Display in QTextEdit (assuming self.textEditPage7 exists)
            self.textEditResults.setHtml(html_table)
        except AttributeError:
            self.textEditClimTree.setText("Please do the sequence alignment before attempting to generate the tree !")
            self.stackedWidget.setCurrentIndex(2)
            self.tabWidget2.setCurrentIndex(2)

    # Enable_button():
    def onTextChanged(self):
        if self.textEditPage1.toPlainText() and self.textEditPage4.toPlainText():
            self.resultsButton.setEnabled(True)
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
        else:
            df = pd.read_csv(Params.file_name)
            climaticTrees = utils.climaticPipeline(df)
            utils.filterResults(climaticTrees, self.geneticTreeDict, df)
            with open("./results/output.csv", "r") as f:
                content = f.read()
                print(content)

    def toggleDarkMode(self):
        self.isDarkMode = not self.isDarkMode
        buttons_Vertical = [self.fileBrowserButtonPage1, self.sequenceAlignmentButtonPage1,
                            self.clearButtonPage1, self.statisticsButtonPage1, self.geneticTreeButtonPage1,
                            self.fileBrowserButtonPage2, self.clearButtonPage2, self.resultsButtonPage2,
                            self.climaticTreeButtonPage2, self.statisticsButtonPage2,
                            self.settingsButtonPage3,
                            self.settingsButtonPage4,
                            self.submitButtonPage3,
                            self.statisticsButtonPage3,
                            self.submitButtonPage4,
                            self.statisticsButtonPage4,
                            self.clearButtonPage3,
                            self.clearButtonPage4]

        if self.isDarkMode:
            qtmodern.styles.dark(app)
            self.darkModeButton.setIcon(QIcon(":other/light.png"))  # Set the 'light' icon for dark mode
            self.darkModeButton.setCursor(Qt.PointingHandCursor)
            self.darkModeButton.setStyleSheet("""
                QPushButton { 
                    padding: 10px 20px;
                    font-weight: bold;
                    background-color: #464645;
                    border-radius: 14px;
                    transition: background-color 0.3s ease; /* Add transition */
                }
                QPushButton:hover {
                    background-color: #B7B7B6; 
                }
                QPushButton:pressed {
                    background-color: #DEDDDA;
                }
            """)
            # Créer l'effet d'ombre
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)  # Ajuster le flou de l'ombre
            shadow_effect.setColor(QColor(0, 0, 0, 140))  # Couleur de l'ombre (noir avec transparence)
            shadow_effect.setOffset(3, 3)  # Décalage de l'ombre

            # Appliquer l'effet d'ombre au bouton
            self.darkModeButton.setGraphicsEffect(shadow_effect)
            #######################################################
            for buttonV in buttons_Vertical:
                buttonV.setCursor(Qt.PointingHandCursor)
                buttonV.setStyleSheet("""
                    QPushButton {
                        color: #EFEFEF;
                        border-radius: 14px;
                        background-color: #464645;
                        padding: 10px 20px;
                        font-weight: bold;
                        transition: background-color 0.3s ease; /* Add transition */
                    }
                    QPushButton:hover {
                        background-color: #9F9F9F; 
                    }
                    QPushButton:pressed {
                        background-color: #EEEEEE;
                    }
                """)
                # Créer l'effet d'ombre
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)  # Ajuster le flou de l'ombre
                shadow_effect.setColor(QColor(0, 0, 0, 110))  # Couleur de l'ombre (noir avec transparence)
                shadow_effect.setOffset(3, 3)  # Décalage de l'ombre

                # Appliquer l'effet d'ombre au bouton
                buttonV.setGraphicsEffect(shadow_effect)

        else:
            qtmodern.styles.light(app)
            self.darkModeButton.setIcon(QIcon(":other/dark.png"))  # Set the 'dark' icon
            self.darkModeButton.setCursor(Qt.PointingHandCursor)
            self.darkModeButton.setStyleSheet("""
                    QPushButton { 
                        padding: 10px 20px;
                        font-weight: bold;
                        background-color: #DEDDDA;
                        border-radius: 14px;
                        transition: background-color 0.3s ease; /* Add transition */
                    }
                    QPushButton:hover {
                        background-color: #B7B7B6; 
                    }
                    QPushButton:pressed {
                        background-color: #DEDDDA;
                    }
                """)
            # Créer l'effet d'ombre
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)  # Ajuster le flou de l'ombre
            shadow_effect.setColor(QColor(0, 0, 0, 140))  # Couleur de l'ombre (noir avec transparence)
            shadow_effect.setOffset(3, 3)  # Décalage de l'ombre

            # Appliquer l'effet d'ombre au bouton
            self.darkModeButton.setGraphicsEffect(shadow_effect)
            ##########################################################
            for buttonV in buttons_Vertical:
                buttonV.setCursor(Qt.PointingHandCursor)
                buttonV.setStyleSheet("""
                    QPushButton {
                        border-radius: 14px;
                        background-color: #EEEEEE;
                        padding: 10px 20px;
                        font-weight: bold;
                        transition: background-color 0.3s ease; /* Add transition */
                    }
                    QPushButton:hover {
                        background-color: #D7D7D7; 
                    }
                    QPushButton:pressed {
                        background-color: #EEEEEE;
                    }
                """)
                # Créer l'effet d'ombre
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)  # Ajuster le flou de l'ombre
                shadow_effect.setColor(QColor(0, 0, 0, 110))  # Couleur de l'ombre (noir avec transparence)
                shadow_effect.setOffset(3, 3)  # Décalage de l'ombre

                # Appliquer l'effet d'ombre au bouton
                buttonV.setGraphicsEffect(shadow_effect)

    # press the button to delete data
    def clearGen(self):
        self.textEditFasta.clear()
        self.textEditSeqAlign.clear()
        self.textEditGenTree.clear()
        self.GenStatsList.setCurrentIndex(0)

    def clearClim(self):
        self.textEditClimData.clear()
        self.textEditClimStats.clear()
        self.textEditClimTree.clear()
        self.graphicsViewClimData.clear()
        self.ClimStatsListCondition.setCurrentIndex(0)
        self.ClimStatsListChart.setCurrentIndex(0)

    def clearResult(self):
        self.textEditResults.clear()

    def clearResultStat(self):
        self.ResultsStatsListCondition.setCurrentIndex(0)
        self.ResultsStatsListChart.setCurrentIndex(0)


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    window = UiMainWindow()
    mw = qtmodern.windows.ModernWindow(window)

    # Get screen geometry
    screen_geometry = app.primaryScreen().availableGeometry()

    # Calculate center position
    center_point = screen_geometry.center()
    x = center_point.x() - mw.width() // 2
    y = center_point.y() - mw.height() // 2

    # Move window to center
    mw.move(x, y)

    mw.show()
    sys.exit(app.exec_())
