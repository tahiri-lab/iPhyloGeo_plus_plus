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
    """
     Custom YAML Dumper to modify the default indentation and list representation behavior.

     Methods:
         increase_indent(flow=False, indentless=False):
             Increase the indentation level in the YAML output.

         represent_list(data):
             Represent Python lists in a flow style in the YAML output.
     """

    def increase_indent(self, flow=False, indentless=False):
        """
        Increase the indentation level in the YAML output.

        Args:
            flow (bool): Indicates whether the current context is a flow style. Defaults to False.
            indentless (bool): Indicates whether to use an indentless format. This argument is ignored. Defaults to False.

        Returns:
            The result from the superclass's increase_indent method with modified behavior.
        """
        return super(MyDumper, self).increase_indent(flow, False)

    def represent_list(self, data):
        """
        Represent Python lists in a flow style in the YAML output.

        Args:
            data (list): The list to represent in the YAML output.

        Returns:
            The YAML representation of the list in a flow style.
        """
        return self.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.add_representer(list, MyDumper.represent_list, Dumper=MyDumper)
"""
Add a representer for the list type to the YAML dumper.

This ensures that lists are represented in a flow style using the MyDumper class.

Args:
    list (type): The Python list type to represent.
    MyDumper.represent_list (method): The method that defines how to represent lists.
    Dumper (yaml.Dumper): The custom dumper class to use, in this case, MyDumper.
"""


def update_yaml_param(params, file_path, property_name, new_value):
    """
    Updates a specified property within a YAML file with a new value.

    Args:
        params: An object with an update_from_dict method, typically used for updating parameters.
        file_path (str): The path to the YAML file.
        property_name (str): The name of the property to modify (e.g., 'file_name').
        new_value: The new value to set for the property (can be any valid YAML type).

    Raises:
        FileNotFoundError: If the specified YAML file does not exist.
        KeyError: If the specified property name is not found in the YAML file.
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
        """
        Initialize and display the 'How to Use' window.

        This method creates a new QMainWindow instance, sets up its UI using the UiHowToUse class, and displays the window.
        """
        self.window = QtWidgets.QMainWindow()
        self.ui = UiHowToUse()
        self.ui.initUI()
        self.ui.show()

    def paramWin(self):
        """
        Initialize and display the parameters window.
        This method creates a new QMainWindow instance, sets up its UI using the UiDialog class, and displays the window.
        """
        self.window = QtWidgets.QMainWindow()
        self.ui = UiDialog()
        self.ui.setupUi(self.window)
        self.window.show()

    def openClimTree(self):
        """
        Initialize and display the climatic tree window.

        This method imports the Ui_ct class, creates a new QMainWindow instance, sets up its UI using the Ui_ct class, and displays the window. It also sets the current index for stackedWidget and tabWidget2.
        """
        from cltree import Ui_ct
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_ct()
        self.ui.setupUi(self.window)
        self.window.show()

        self.stackedWidget.setCurrentIndex(2)
        self.tabWidget2.setCurrentIndex(3)

    def __init__(self):
        super(UiMainWindow, self).__init__()
        uic.loadUi("Qt/main.ui", self)
        self.setupUi()

    def setupUi(self):
        """
         Setup the UI components and initialize the main window.

         This method connects various UI buttons to their corresponding event handlers, sets up styles and effects for UI elements, and initializes the state of the application.
         """
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
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)
            shadow_effect.setColor(QColor(0, 0, 0, 140))
            shadow_effect.setOffset(3, 3)
            button.setGraphicsEffect(shadow_effect)

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
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)
            shadow_effect.setColor(QColor(0, 0, 0, 110))
            shadow_effect.setOffset(3, 3)

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
        shadow_effect.setBlurRadius(10)
        shadow_effect.setColor(QColor(0, 0, 0, 140))
        shadow_effect.setOffset(3, 3)
        self.darkModeButton.setGraphicsEffect(shadow_effect)
        QtCore.QMetaObject.connectSlotsByName(self)

    def pressItFasta(self):
        """
        Open a dialog to select a FASTA file, update parameters, and display the content with color-coded sequences.

        This method allows the user to select a FASTA file from the file system. It updates the relevant YAML parameters
        with the file's path and name, reads the file content, and displays the sequences with color-coded nucleotides
        in a text edit widget.

        Actions:
            - Opens a file dialog to select a FASTA file.
            - Updates 'reference_gene_file' and 'reference_gene_dir' parameters in the YAML file.
            - Reads and displays the content of the selected FASTA file.
            - Color-codes nucleotides (A, C, G, T) in the displayed sequence.
            - Enables the sequence alignment button and updates icons.
        """
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
                sequence = ""
                for line in content.splitlines():
                    if line.startswith('>'):
                        line = f'<span style="color: green; font-weight: bold; font-size: 20px;">{line}</span>'
                        sequence += "<br>" + line + "<br>"
                    else:
                        nucleotide_colors = {
                            'A': 'yellow',
                            'C': 'blue',
                            'G': 'red',
                            'T': 'orange'
                        }
                        colored_line = ''
                        for char in line:
                            color = nucleotide_colors.get(char, '')
                            if color:
                                colored_line += f'<span style="color: {color}; font-weight: bold; font-size: 20px;">{char}</span>'
                            else:
                                colored_line += char
                        line = colored_line
                        sequence += line

                self.textEditFasta.setHtml(
                    f"<div style='background-color: #000000; color: #ffffff; padding: 10px; white-space: pre-wrap; word-wrap: break-word;'>{sequence}</div>"
                )
                self.tabWidget.setCurrentIndex(1)
                self.sequenceAlignmentButtonPage1.setEnabled(True)
                self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
                self.sequenceAlignmentButtonPage1.setIcon(QIcon(":inactive/sequence.svg"))

    def SeqAlign(self):
        """
        Perform sequence alignment and store the resulting genetic tree dictionary.

        This method calls the callSeqAlign method to perform sequence alignment and stores the resulting genetic tree dictionary in the geneticTreeDict attribute.
        """
        self.geneticTreeDict = self.callSeqAlign()

    def callSeqAlign(self):
        """
        Execute the sequence alignment pipeline and display progress.

        This method performs the following steps:
        1. Loads sequences from the reference gene file.
        2. Aligns the loaded sequences.
        3. Generates genetic trees based on the aligned sequences.
        4. Prepares and displays the results in the UI.
        5. Saves the alignment and genetic tree results to JSON files.

        Returns:
            dict: A dictionary containing the genetic trees.

        Raises:
            Exception: Any exception that occurs during the alignment process will be handled, and the loading screen will close.
        """
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

    def retrieveDataNames(self, list):
        """
        Retrieve data from a list, excluding the first element.

        Args:
            list (list): The list to retrieve data from.

        Returns:
            list: A list of data excluding the first element.
        """
        names_to_retrieve = []
        for data in list:
            if data != list[0]:
                names_to_retrieve.append(data)
        return names_to_retrieve

    from decimal import Decimal
    from PyQt5.QtWebEngineWidgets import QWebEngineView
    import io
    import folium
    from PyQt5 import QtWidgets, QtGui
    from PyQt5.QtWidgets import QFileDialog

    def pressItCSV(self):
        """
        Retrieve data from a climatic file and display it in a table.

        This method allows the user to select a CSV file from the file system. It updates the relevant YAML parameters with the file's path, reads the file content, and displays the data in a table widget. It also processes location data (latitude and longitude) if available and populates a map.

        Actions:
            - Opens a file dialog to select a CSV file.
            - Updates 'file_name' parameter in the YAML file.
            - Reads and displays the content of the selected CSV file.
            - Processes and stores species and factor data.
            - Populates a map with location data if available.
            - Updates the UI to reflect the loaded data.
        """

        def is_valid_decimal(value):
            try:
                Decimal(value)
                return True
            except:
                return False

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
                    first_line_without_loc = first_line[:-2]
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
                self.tabWidget2.setCurrentIndex(1)
                for i, line in enumerate(lines):
                    line_split = line.split(",")
                    for value in line_split:
                        if i == 0:
                            cursor.setCharFormat(format)
                        cursor.insertText(value.strip())
                        cursor.movePosition(QtGui.QTextCursor.NextCell)
                    if i > 0 and loc:
                        lat_value = line_split[-2].strip()
                        long_value = line_split[-1].strip()
                        if is_valid_decimal(lat_value) and is_valid_decimal(long_value):
                            lat.append(Decimal(lat_value))
                            long.append(Decimal(long_value))
                        else:
                            print(f"Invalid coordinate data at row {i}: {lat_value}, {long_value}")

                if loc and lat and long:
                    self.populateMap(lat, long)


    def populateMap(self, lat, long):
        """
        Create and display a folium map with given latitude and longitude.

        This method generates a map centered on the mean latitude and longitude of the provided coordinates.
        It places markers on the map for each coordinate pair and displays the map in a QWebEngineView.

        Args:
            lat (list): List of latitudes.
            long (list): List of longitudes.
        """
        mean_lat = sum(lat) / len(lat)
        mean_long = sum(long) / len(long)

        m = folium.Map(location=[mean_lat, mean_long],
                       zoom_start=14,
                       tiles="OpenStreetMap")
        for latitude, longitude in zip(lat, long):
            folium.Marker([latitude, longitude]).add_to(m)

        data = io.BytesIO()
        m.save(data, close_file=False)

        web_view = QWebEngineView(self.graphicsViewClimData)  # Embed the map inside graphicsViewClimData
        web_view.setHtml(data.getvalue().decode())
        layout = QtWidgets.QVBoxLayout(self.graphicsViewClimData)
        layout.addWidget(web_view)
        self.graphicsViewClimData.setLayout(layout)

    def showClimStatBarAllFact(self):
        """
        Generate a bar graph that includes every factor for every species.

        This method processes the climatic factors for each species and generates a bar graph to visualize the distribution of climatic variables. It handles cases where the data is not available and updates the UI accordingly.

        Raises:
            AttributeError: If the factors attribute does not exist, an error message is displayed, and the UI is updated.
        """
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
        """
        Display the home page of the application.

        This method sets the icons for the climatic data and genetic data buttons to their inactive states and displays the home page by setting the stacked widget's current index to 0.
        """
        self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
        self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
        self.homeButton.setIcon(QIcon(":active/home.png"))

        self.stackedWidget.setCurrentIndex(0)

    def showGenDatPage(self):
        """
        Display the genetic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the genetic data page by setting the stacked widget's current index to 1, and sets the tab widget's current index to 0.
        """
        self.climaticDataButton.setIcon(QIcon(":inactive/climatic.svg"))
        self.geneticDataButton.setIcon(QIcon(":active/genetic.svg"))
        self.homeButton.setIcon(QIcon(":other/home.svg"))
        self.stackedWidget.setCurrentIndex(1)
        self.tabWidget.setCurrentIndex(0)

    def showClimDatPage(self):
        """
        Display the climatic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the climatic data page by setting the stacked widget's current index to 2, and sets the tab widget's current index to 0.
        """
        self.climaticDataButton.setIcon(QIcon(":active/climatic.svg"))
        self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
        self.homeButton.setIcon(QIcon(":other/home.svg"))
        self.stackedWidget.setCurrentIndex(2)
        self.tabWidget2.setCurrentIndex(0)

    def showResultsPage(self):
        """
        Display the results page of the application.

        This method sets the stacked widget's current index to 3 to display the results page.
        """
        self.stackedWidget.setCurrentIndex(3)

    def showResultsStatsPage(self):
        """
        Display the results statistics page of the application.

        This method sets the stacked widget's current index to 4 to display the results statistics page.
        """
        self.stackedWidget.setCurrentIndex(4)

    def showFilteredResults(self):
        """
        Show the results filtered with a metric threshold provided by the user.

        This method reads the data from a CSV file, processes it through the climatic pipeline, filters the results, and displays the filtered results in an HTML table format within a QTextEdit widget. It handles exceptions related to missing sequence alignment.

        Raises:
            AttributeError: If the sequence alignment has not been performed before attempting to generate the tree.
        """
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

    def onTextChanged(self):
        """
        Handle text change events in the text edit widgets.

        This method enables or disables the results button based on whether both textEditPage1 and textEditPage4 have content. If the text edit widgets have content, the results button is enabled and its icon is set. If not, it processes the CSV file, generates climatic trees, filters results, and prints the content of the output CSV file.

        Actions:
            - Enable the results button if both text edit widgets have content.
            - Process the CSV file, generate climatic trees, filter results, and print output if the text edit widgets do not have content.
        """
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
        """
        Toggle the application's dark mode setting.

        This method switches the application's theme between dark mode and light mode. It updates the isDarkMode attribute, applies the corresponding style, and changes the icon of the darkModeButton.

        Attributes:
            isDarkMode (bool): A flag indicating whether dark mode is currently enabled.

        Actions:
            If dark mode is enabled, apply the dark style and set the darkModeButton icon to the 'light' icon.
            If dark mode is disabled, apply the light style and set the darkModeButton icon to the 'dark' icon.
        """

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
            self.top_frame.setStyleSheet("background-color: #646464;")
            self.homeButton.setStyleSheet("background-color: #838383;")
            self.geneticDataButton.setStyleSheet("background-color: #838383;")
            self.climaticDataButton.setStyleSheet("background-color: #838383;")
            self.helpButton.setStyleSheet("background-color: #838383;")

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
            self.top_frame.setStyleSheet("background-color: rgb(222, 221, 218);")
            self.homeButton.setStyleSheet("background-color: rgb(222, 221, 218);")
            self.geneticDataButton.setStyleSheet("background-color: rgb(222, 221, 218);")
            self.climaticDataButton.setStyleSheet("background-color: rgb(222, 221, 218);")
            self.helpButton.setStyleSheet("background-color: rgb(222, 221, 218);")
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
        """
        Clear the genetic data fields.

        This method clears the content of the text edit widgets related to FASTA sequences, sequence alignments, and genetic trees. It also resets the current index of the genetic statistics list to 0.
        """
        self.textEditFasta.clear()
        self.textEditSeqAlign.clear()
        self.textEditGenTree.clear()
        self.GenStatsList.setCurrentIndex(0)

    def clearClim(self):
        """
            Clear the text fields related to climatic data.
            """
        self.textEditClimData.clear()
        self.textEditClimStats.clear()
        self.textEditClimTree.clear()
        self.graphicsViewClimData.clear()
        self.ClimStatsListCondition.setCurrentIndex(0)
        self.ClimStatsListChart.setCurrentIndex(0)

    def clearResult(self):
        """
        Clear the text fields related to climatic data.

        This method clears the content of the text edit widgets for climatic data, climatic statistics, and climatic trees. It also clears the graphics view for climatic data and resets the current index of the climatic statistics and chart lists to 0.
        """
        self.textEditResults.clear()

    def clearResultStat(self):
        """
        Clear the statistics result lists.

        This method resets the current index of the results statistics condition list and the results statistics chart list to 0.
        """
        self.ResultsStatsListCondition.setCurrentIndex(0)
        self.ResultsStatsListChart.setCurrentIndex(0)


if __name__ == "__main__":
    # Create the application instance
    app = QtWidgets.QApplication([])

    # Create the main window instance
    window = UiMainWindow()

    # Wrap the main window with the ModernWindow style
    mw = qtmodern.windows.ModernWindow(window)

    # Get screen geometry to determine the available screen space
    screen_geometry = app.primaryScreen().availableGeometry()

    # Calculate the center position of the screen
    center_point = screen_geometry.center()
    x = center_point.x() - mw.width() // 2
    y = center_point.y() - mw.height() // 2

    # Move the main window to the center of the screen
    mw.move(x, y)

    # Show the main window
    mw.show()

    # Execute the application's event loop
    sys.exit(app.exec_())

