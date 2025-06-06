import pandas as pd
from aphylogeo import utils
from aphylogeo.params import Params
from event_connector import blocked_signals
from utils.custom_table import create_sleek_table
from utils.download_file import download_file_local
from utils.error_dialog import show_error_dialog
from PyQt6 import QtGui, QtWidgets
from utils.result_settings_dialog import ResultSettingsDialog
from utils.tree_graph import generate_tree_with_bar, init_tree
from utils.tree_map import generate_tree_map

RESULT_TAB = 0
STATS_TAB = 1	
MAP_TAB = 2

class Result:
    def __init__(self, main):
        self.main = main
        self.table = None

    def open_result_settings_window(self):
        """
        Initialize and display the parameters window.
        This method creates a new QMainWindow instance, sets up its UI using the UiDialog class, and displays the window.
        """

        dialog = ResultSettingsDialog()
        dialog.exec()

    def show_filtered_results(self):
        """
        Show the results filtered with a metric threshold provided by the user.

        This method reads the data from a CSV file, processes it through the climatic pipeline, filters the results,
        and displays the filtered results in an HTML table format within a QTextBrowser widget.
        It handles exceptions related to missing sequence alignment.

        Raises:
            AttributeError: If the sequence alignment has not been performed before attempting to generate the tree.
        """
        self.main.tabWidgetResult.setCurrentIndex(0)

        df = self.main.climatePage.climat.data.copy()
        try:
            utils.filterResults(self.main.climatePage.climat.climaticTree.climaticTrees, self.main.geneticsPage.genetics.geneticAlignment.geneticTrees, df)
        except Exception:
            show_error_dialog("The data given is not correct, make sure that you loaded the correct files in the previous steps.")
            return
        df_results = pd.read_csv("./results/output.csv")
        df_results["Name of species"] = df_results["Name of species"].str.replace("_", " ")
        # Replace the first column values with Params.file_name just before visualization
        df_results.iloc[:, 0] = Params.reference_gene_file

        if self.table is not None:
            self.main.resultTableLayout.removeWidget(self.table)
            self.table.deleteLater()
            self.table = None

        self.table = create_sleek_table(df_results)
        self.main.resultTableLayout.addWidget(self.table)

    def clear_result(self):
        """
        Clear the text fields related to climatic data.

        This method clears the content of the text edit widgets for climatic data, climatic statistics, and climatic trees.
        It also clears the graphics view for climatic data and resets the current index of the climatic statistics and chart lists to 0.
        """
        try:
            self.main.textEditResults.clear()
            self.main.geneticsPage.genetics.clear_genetic_data()
            self.main.climatePage.climat.clear_climmatic_data()
            self.main.stackedWidget.setCurrentIndex(0)
            self.main.resultTableLayout.removeWidget(self.table)
            self.table.deleteLater()
            self.table = None
            self.main.resultsButton.setEnabled(False)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}", "Error")

    def on_tab_changed(self, index):
        if index == STATS_TAB:
            self.display_phylogeographic_trees()
        if index == MAP_TAB:
            self.show_map()

    def show_map(self):
        self.main.tabWidgetResult.setCurrentIndex(2)
        ### LOADING ###
        def update_progress_bar(value):
                self.main.mapImageBar.setValue(value)
                QtWidgets.QApplication.processEvents()
        if self.main.mapImage.pixmap() is None:

            self.main.mapImageBar.setVisible(True)
            QtWidgets.QApplication.processEvents()
        ### LOADING ###

        try:
            # Get Genetic data
            trees = self.main.geneticsPage.genetics.geneticStatistics.geneticAlignment.geneticTrees
            firstTreeKey, _ = next(iter(trees.items()))
            currentTree = self.tree_keys[self.main.phyloTreescomboBox.currentIndex()] if self.main.phyloTreescomboBox.currentIndex() != -1 else firstTreeKey
            treeData = trees.get(currentTree) or trees.get(firstTreeKey)
            # Get Climatic data
            gpsFile = self.main.climatePage.climat.data.copy()
        except Exception:
            show_error_dialog("The data given is not correct, make sure that you loaded the correct files in the previous steps.")
            return
        
        generate_tree_map(treeData, gpsFile, progress_callback=update_progress_bar)
        self.main.mapImage.setPixmap(QtGui.QPixmap("results/figure_tree2map.png"))
        self.main.mapImageBar.setVisible(False)

    def display_phylogeographic_trees(self):
        self.jsonData, self.tree_keys, formatted_tree_keys = init_tree()

        with blocked_signals(self.main.phyloTreescomboBox, self.main.criteriaComboBox, self.main.tabWidgetResult):
            self.main.tabWidgetResult.setCurrentIndex(1)

            self.main.phyloTreescomboBox.clear()
            self.main.phyloTreescomboBox.addItems(formatted_tree_keys)

            self.climatic_data = self.main.climatePage.climat.data.copy()
            self.main.criteriaComboBox.clear()
            self.main.criteriaComboBox.addItems(self.climatic_data.columns[1:])

        self.render_tree()

    def render_tree(self):
        self.key = self.tree_keys[self.main.phyloTreescomboBox.currentIndex()]

        pixmap = generate_tree_with_bar(self.key, self.jsonData, self.main.criteriaComboBox.currentText(), self.climatic_data, self.main.isDarkMode)

        self.main.PhyloTreeLabel.clear()
        self.main.PhyloTreeLabel.setPixmap(pixmap)
        self.main.PhyloTreeLabel.adjustSize()

    def save_tree_graph(self):
        download_file_local(self.key)
