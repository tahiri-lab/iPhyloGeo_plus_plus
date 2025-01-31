import json
import os
import toytree
import toyplot.png
import numpy as np
from numpy.typing import ArrayLike
from typing import cast

from utils.result_settings_dialog import ResultSettingsDialog
import pandas as pd
from aphylogeo import utils
from aphylogeo.params import Params
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QFileDialog, QVBoxLayout
from utils.error_dialog import show_error_dialog

class Result:
    
    def __init__(self, main):
        self.main = main
        
    # Initializes the result page tabs and connects the tab change event to the on_tab_changed method
    def initialize_result_page(self):
        self.main.tabWidgetResult.currentChanged.connect(self.on_tab_changed)
        self.main.tabWidgetResult.setCurrentIndex(0)
        
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
        utils.filterResults(self.main.climat.climaticTrees, self.main.genetics.geneticTreeDict, df)
        df_results = pd.read_csv("./scripts/results/output.csv")
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
                transition: all 0.3s ease-in-out;
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
            # self.graphicsViewClimData.clear()
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}", "Error")
            
    def on_tab_changed(self, index):
        if index == 1:
            self.handle_tab_statistic_refresh()     

    def handle_tab_statistic_refresh(self):
        self.display_phylogeographic_trees()
            
    def display_phylogeographic_trees(self):
        self.main.tabWidgetResult.setCurrentIndex(1)
        self.results_dir = "scripts/results"
        file_path = os.path.join(self.results_dir, "geneticTrees.json")

        with open(file_path, "r") as file:
            self.phylo_json = json.load(file)

        self.tree_keys = list(self.phylo_json.keys())
        self.total_trees = len(self.tree_keys)

        self.current_index1 = 0
        self.main.phyloTreescomboBox.clear()

        formatted_tree_keys = [self.rename_tree_key(key) for key in self.tree_keys]
        self.main.phyloTreescomboBox.addItems(formatted_tree_keys)

        self.climatic_data = pd.read_csv(Params.file_name)
        self.main.criteriaComboBox.clear()
        self.main.criteriaComboBox.addItems(self.climatic_data.columns[1:])

        self.render_tree(self.current_index1)

    def rename_tree_key(self, tree_key):
        parts = tree_key.split("_")
        if len(parts) == 2:
            return f"{parts[0]} nt {parts[1]} nt"
        return tree_key

    def render_tree(self, index):
        if 0 <= index < self.total_trees:
            self.current_index1 = index
            key = self.tree_keys[index]
            newick_str = self.phylo_json[key]

            tree = toytree.tree(newick_str)
            tip_labels = tree.get_tip_labels()

            if not hasattr(self, "climatic_data"):
                return

            selected_column = self.main.criteriaComboBox.currentText()

            # Match the original labels from the CSV
            ordered_climatic_data = self.climatic_data.set_index("id").reindex(tip_labels).reset_index()
            bar_values = ordered_climatic_data[selected_column].values

            # Ensure bar_values have no negative values
            bar_values = np.clip(cast(ArrayLike, bar_values), a_min=0, a_max=None)

            # Create a canvas for the plot
            canvas = toyplot.Canvas(width=921, height=450)

            # Add tree to canvas
            ax0 = canvas.cartesian(bounds=(50, 300, 50, 400), ymin=0, ymax=tree.ntips, padding=15)

            # Replace underscores with spaces for display purposes only
            tree.draw(axes=ax0, tip_labels=[label.replace("_", " ") for label in tip_labels], tip_labels_colors="black")
            ax0.show = False

            # Add bar plot to canvas
            ax1 = canvas.cartesian(bounds=(325, 900, 50, 400), ymin=0, ymax=tree.ntips, padding=15)
            ax1.bars(
                np.arange(tree.ntips),
                bar_values,
                along="y",
            )

            # Adjust the appearance for better visibility
            ax1.show = True
            ax1.y.show = False
            ax1.x.ticks.show = True
            ax1.x.label.text = selected_column

            # Save and display the plot
            self.temp_img_path = os.path.join(self.results_dir, f"{key}.png")
            toyplot.png.render(canvas, self.temp_img_path)

            pixmap = QPixmap(self.temp_img_path)
            self.main.PhyloTreeLabel.clear()
            self.main.PhyloTreeLabel.setPixmap(pixmap)
            self.main.PhyloTreeLabel.adjustSize()

    def save_tree_graph(self):
        current_key = self.tree_keys[self.current_index1]
        default_file_name = f"{current_key}.png"

        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Graph As",
            default_file_name,
            "PNG Files (*.png);;All Files (*)",
            options=options,
        )

        if file_path:
            if not file_path.lower().endswith(".png"):
                file_path += ".png"

            with open(self.temp_img_path, "rb") as temp_file:
                with open(file_path, "wb") as file:
                    file.write(temp_file.read())