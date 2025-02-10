import pandas as pd
from aphylogeo import utils
from aphylogeo.params import Params
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import QVBoxLayout

from utils.error_dialog import show_error_dialog
from utils.result_settings_dialog import ResultSettingsDialog
from utils.download_file import download_file_local
from utils.tree_graph import init_tree, generate_tree_with_bar
from event_connector import blocked_signals


class Result:
    def __init__(self, main):
        self.main = main
        
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
        self.main.stackedWidget.setCurrentIndex(3)

        df = pd.read_csv(Params.file_name)
        try:
            utils.filterResults(self.main.climat.climaticTree.climaticTrees, self.main.genetics.geneticTreeDict, df)
        except Exception as e:
            show_error_dialog(str(e), "Aphylogeo Utils Error")
            return
        df_results = pd.read_csv("./results/output.csv")
        df_results["Name of species"] = df_results["Name of species"].str.replace("_", " ")
        # Replace the first column values with Params.file_name just before visualization
        df_results.iloc[:, 0] = Params.reference_gene_file

        # Convert to HTML table with sleek, modern styling
        html_table = df_results.to_html(index=False, border=0, classes="dataframe", table_id="styled-table")

        # Add CSS styles for a professional look with smooth hover animations
        html_style = """
        <style>
            #styled-table {
                font-family: 'Trebuchet MS', Arial, Helvetica, sans-serif;
                border-collapse: collapse;
                width: 100%;
                border-radius: 5px;
                box-shadow: 0 2px 5px rgba(0, 0, 0, 0.15);
                overflow: hidden;
            }
            #styled-table td, #styled-table th {
                border: 1px solid #ddd;
                padding: 12px;
            }
            #styled-table tr:nth-child(even) {
                background-color: #f9f9f9;
            }
            #styled-table tr:hover {
                background-color: #f1f1f1;
                transform: scale(1.01);
            }
            #styled-table th {
                padding-top: 12px;
                padding-bottom: 12px;
                text-align: left;
                background-color: #4CAF50;
                color: white;
            }
            #styled-table td {
                padding-left: 12px;
                padding-right: 12px;
            }
        </style>
        """

        # Combine the HTML style and table
        html_content = html_style + html_table

        # Set up the QWebEngineView
        web_engine_view = QWebEngineView()
        web_engine_view.setHtml(html_content)

        # Clear the existing content and layout of the QTextBrowser (if any)
        layout = QVBoxLayout(self.main.textEditResults)
        for i in reversed(range(layout.count())):
            widget_to_remove = layout.itemAt(i).widget()
            layout.removeWidget(widget_to_remove)
            widget_to_remove.setParent(None)

        # Add the QWebEngineView to the QTextBrowser
        layout.addWidget(web_engine_view)

    def clear_result(self):
        """
        Clear the text fields related to climatic data.

        This method clears the content of the text edit widgets for climatic data, climatic statistics, and climatic trees.
        It also clears the graphics view for climatic data and resets the current index of the climatic statistics and chart lists to 0.
        """
        try:
            self.main.textEditResults.clear()
            self.main.genetics.clear_genetic_data()
            self.main.climat.clear_climmatic_data()
            self.main.stackedWidget.setCurrentIndex(0)
            self.main.resultsButton.setEnabled(False)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}", "Error")

    def on_tab_changed(self, index):
        if index == 1:
            self.display_phylogeographic_trees()     
            
    def display_phylogeographic_trees(self):
        self.jsonData, self.tree_keys, formatted_tree_keys = init_tree()
        
        with blocked_signals(self.main.phyloTreescomboBox, self.main.criteriaComboBox, self.main.tabWidgetResult):
            self.main.tabWidgetResult.setCurrentIndex(1)
            
            self.main.phyloTreescomboBox.clear()
            self.main.phyloTreescomboBox.addItems(formatted_tree_keys)

            self.climatic_data = pd.read_csv(Params.file_name)
            self.main.criteriaComboBox.clear()
            self.main.criteriaComboBox.addItems(self.climatic_data.columns[1:])        

        self.render_tree()

    def render_tree(self):
        self.key = self.tree_keys[self.main.phyloTreescomboBox.currentIndex()]
    
        pixmap = generate_tree_with_bar(self.key, self.jsonData, self.main.criteriaComboBox.currentText(), self.climatic_data)

        self.main.PhyloTreeLabel.clear()
        self.main.PhyloTreeLabel.setPixmap(pixmap)
        self.main.PhyloTreeLabel.adjustSize()
        


    def save_tree_graph(self):
        download_file_local(self.key, self.main)

