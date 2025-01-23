import io
import json
import os
import re
import shutil
import sys
from decimal import Decimal

import folium
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import qtmodern.styles
import qtmodern.windows
import seaborn as sns
import toyplot.png
import toytree
from aphylogeo import utils
from aphylogeo.params import Params
from Bio import Phylo
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor, QIcon, QPixmap
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QDialog, QFileDialog, QGraphicsDropShadowEffect, QTableWidget, QTableWidgetItem, QVBoxLayout
from ui_helpers import create_shadow_effect, get_button_style, style_buttons
from utils import resources_rc  # noqa: F401  # Import the compiled resource module for resolving image resource path
from utils.genetic_params_dialog import ParamDialog
from utils.help import HelpDialog
from utils.MyDumper import update_yaml_param
from utils.PreferencesDialog import PreferencesDialog  # Import PreferencesDialog
from utils.resultSettingsDialog import ResultSettingsDialog
from utils.settings import Params2
from genetics import Genetics
from event_connector import QtEvents, connect_event, connect_decorated_methods

try:
    Params.load_from_file("./scripts/utils/params.yaml")
except FileNotFoundError:
    Params.validate_and_set_params(Params2.PARAMETER_KEYS)

window_size = 50
starting_position = 1


class UiMainWindow(QtWidgets.QMainWindow):
    def open_help_window(self):
        """
        Initialize and display the 'How to Use' window.

        This method creates a new QMainWindow instance, sets up its UI using the UiHowToUse class, and displays the window.
        """

        dialog = HelpDialog()
        dialog.exec_()

    def open_result_settings_window(self):
        """
        Initialize and display the parameters window.
        This method creates a new QMainWindow instance, sets up its UI using the UiDialog class, and displays the window.
        """

        dialog = ResultSettingsDialog()
        dialog.exec_()

    def show_error_dialog(self, message, title="error"):
        """
        Display a professional error dialog with the given title and message.

        Args:
            title (str): The title of the error dialog.
            message (str): The error message to display.
        """
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setWindowTitle(title)
        msgBox.setText(message)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.Ok)
        msgBox.exec_()

    def __init__(self):
        super(UiMainWindow, self).__init__()
        self.genetics = Genetics(self)
        uic.loadUi("scripts/Qt/main.ui", self)
        connect_decorated_methods(self)
        self.setup_ui()

    def setup_ui(self):
        """
        Setup the UI components and initialize the main window.

        This method connects various UI buttons to their corresponding event handlers,
        sets up styles and effects for UI elements, and initializes the state of the application.
        """
        try:
            self.preferences = {
                "label_color": "black",
                "edge_color": "blue",
                "reticulation_color": "red",
                "layout": "horizontal",
                "proportional_edge_lengths": False,
                "label_internal_vertices": False,
                "use_leaf_names": True,
                "show_branch_length": False,
                "view_type": "network",
            }
            self.tree_keys = []
            self.total_trees = 0
            self.current_index = 0
            self.current_index1 = 0
            self.setObjectName("MainWindow")
            self.window_size_spinbox_2.setRange(1, 1000)
            self.starting_position_spinbox_2.setRange(1, 1000)
            self.starting_position_spinbox_2.valueChanged.connect(self.genetics.update_plot)
            self.window_size_spinbox_2.valueChanged.connect(self.genetics.update_plot)
            self.homeButton.clicked.connect(self.show_home_section)
            self.geneticDataButton.clicked.connect(self.show_genetic_section)
            self.clearButtonPage3.clicked.connect(self.clear_results)
            self.climaticDataButton.clicked.connect(self.show_climate_section)
            self.helpButton.clicked.connect(self.open_help_window)
            self.darkModeButton.clicked.connect(self.toggle_dark_mode)
            self.climaticTreeButtonPage2.clicked.connect(self.display_climatic_trees)
            self.climaticTreescomboBox.currentIndexChanged.connect(self.show_selected_climatic_tree)
            self.downloadGraphButton2.clicked.connect(self.download_climatic_tree_graph)
            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
            self.isDarkMode = False  # Keep track of the state
            self.downloadResultsPlotButton.clicked.connect(self.save_tree_graph)
            self.fileBrowserButtonPage1.clicked.connect(self.genetics.select_fasta_file)
            self.geneticTreeButtonPage1.clicked.connect(self.genetics.display_newick_trees)
            self.statisticsButtonPage3.clicked.connect(self.display_phylogeographic_trees)
            self.sequenceAlignmentButtonPage1.clicked.connect(self.genetics.show_sequence_alignment_page)
            self.clearButtonPage1.clicked.connect(self.genetics.clear_genetic_data)
            self.downloadSimilarityButton.clicked.connect(self.genetics.download_similarity_plot_chart)
            self.statisticsButtonPage1.clicked.connect(self.genetics.initialize_species_list)
            self.clearButtonPage2.clicked.connect(self.clear_climmatic_data)
            self.fileBrowserButtonPage2.clicked.connect(self.load_csv_climate_file)
            self.resultsButton.clicked.connect(self.show_results_section)
            self.statisticsButtonPage2.clicked.connect(self.load_climate_statistics)
            self.ClimaticChartSettingsAxisX.currentIndexChanged.connect(self.generate_climate_graph)
            self.ClimaticChartSettingsAxisY.currentIndexChanged.connect(self.generate_climate_graph)
            self.PlotTypesCombobox.currentIndexChanged.connect(self.generate_climate_graph)
            self.climatePlotDownloadButton.clicked.connect(self.download_climate_plot)
            self.geneticTreescomboBox.currentIndexChanged.connect(self.genetics.show_tree)
            self.criteriaComboBox.currentIndexChanged.connect(self.render_tree)
            self.phyloTreescomboBox.currentIndexChanged.connect(self.render_tree)
            self.StartSequenceAlignmentButton.clicked.connect(self.genetics.start_alignment_analysis)
            self.settingsButtonPage3.clicked.connect(self.open_result_settings_window)
            self.submitButtonPage3.clicked.connect(self.show_filtered_results)
            self.clearButtonPage4.clicked.connect(self.clear_result)
            self.statisticsButtonPage4.clicked.connect(self.display_phylogeographic_trees)
            self.clearButtonPage4.clicked.connect(self.clear_result_stat)
            self.downloadGraphButton.clicked.connect(self.genetics.download_genetic_tree_graph)
            self.preferencesButton.clicked.connect(self.open_climatic_tree_preferences_window)
            self.geneticSettingsButton.clicked.connect(self.open_genetic_settings_window)

            self.stackedWidget.setCurrentIndex(0)

            buttons = [
                self.geneticDataButton,
                self.climaticDataButton,
                self.helpButton,
                self.homeButton,
                self.resultsButton,
            ]
            buttons_Vertical = [
                self.fileBrowserButtonPage1,
                self.sequenceAlignmentButtonPage1,
                self.clearButtonPage1,
                self.statisticsButtonPage1,
                self.geneticTreeButtonPage1,
                self.fileBrowserButtonPage2,
                self.clearButtonPage2,
                self.climaticTreeButtonPage2,
                self.statisticsButtonPage2,
                self.settingsButtonPage3,
                self.settingsButtonPage4,
                self.submitButtonPage3,
                self.statisticsButtonPage3,
                self.submitButtonPage4,
                self.statisticsButtonPage4,
                self.StartSequenceAlignmentButton,
                self.clearButtonPage3,
                self.clearButtonPage4,
            ]

            # Define cursor and stylesheet for all buttons
            for button in buttons:
                button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
                button.setStyleSheet(
                    """
                    QPushButton {
                        padding: 10px 20px;
                        font-weight: bold;
                        background-color: #DEDDDA;
                        border-radius: 20px;
                        transition: background-color 0.3s ease;
                    }
                    QPushButton:hover {
                        background-color: #B7B7B6;
                    }
                    QPushButton:pressed {
                        background-color: #DEDDDA;
                    }
                """
                )
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 140))
                shadow_effect.setOffset(3, 3)
                button.setGraphicsEffect(shadow_effect)

            for buttonV in buttons_Vertical:
                buttonV.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
                buttonV.setStyleSheet(
                    """
                    QPushButton {
                        border-radius: 14px;
                        background-color: #EEEEEE;
                        padding: 10px 20px;
                        font-weight: bold;
                        transition: background-color 0.3s ease;
                    }
                    QPushButton:hover {
                        background-color: #D7D7D7;
                    }
                    QPushButton:pressed {
                        background-color: #EEEEEE;
                    }
                """
                )
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 110))
                shadow_effect.setOffset(3, 3)
                buttonV.setGraphicsEffect(shadow_effect)

            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
            self.darkModeButton.setStyleSheet(
                """
                QPushButton {
                    padding: 10px 20px;
                    font-weight: bold;
                    background-color: #DEDDDA;
                    border-radius: 14px;
                    transition: background-color 0.3s ease;
                }
                QPushButton:hover {
                    background-color: #B7B7B6;
                }
                QPushButton:pressed {
                    background-color: #DEDDDA;
                }
            """
            )
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)
            shadow_effect.setColor(QColor(0, 0, 0, 140))
            shadow_effect.setOffset(3, 3)
            self.darkModeButton.setGraphicsEffect(shadow_effect)

            QtCore.QMetaObject.connectSlotsByName(self)

        except AttributeError as e:
            self.show_error_dialog(f"An error occurred while setting up the UI: {e}", "Attribute Error")
        except Exception as e:
            self.show_error_dialog(
                f"An unexpected error occurred: {e}",
                "Unexpected Error",
            )

    def clear_results(self):
        self.textEditResults.clear()

    def open_genetic_settings_window(self):
        dialog = ParamDialog()
        if dialog.exec_() == QDialog.Accepted:
            self.geneticParam = dialog.params
            for property_name, new_value in self.geneticParam.items():
                update_yaml_param(Params, "scripts/utils/params.yaml", property_name, new_value)

    ################################

    def load_climate_statistics(self):
        """
        Load climate data from a CSV file, update the UI elements with the column names, and switch to the appropriate tab.

        This method reads the climate data from a CSV file (excluding the first column), updates combo boxes
        for X and Y axis settings with the column names, and switches to the climate data tab.

        Returns:
            None
        """
        # Load the data
        self.data = pd.read_csv(Params.file_name)
        self.columns = self.data.columns.tolist()

        # Remove the first column
        if self.columns:
            self.columns.pop(0)

        self.ClimaticChartSettingsAxisX.clear()
        self.ClimaticChartSettingsAxisX.addItems(self.columns)
        self.ClimaticChartSettingsAxisY.clear()
        self.ClimaticChartSettingsAxisY.addItems(self.columns)
        self.tabWidget2.setCurrentIndex(2)

    def generate_climate_graph(self):
        """
        Generate and display a graph based on the selected X and Y axis data and the chosen plot type.

        This method reads the selected data columns and plot type from the UI, generates the corresponding graph,
        and displays it in the specified QLabel widget.

        Returns:
            None
        """
        x_data = self.ClimaticChartSettingsAxisX.currentText()
        y_data = self.ClimaticChartSettingsAxisY.currentText()
        plot_type = self.PlotTypesCombobox.currentText()

        # Check if valid options are selected
        if plot_type == "" or x_data == "" or y_data == "":
            return

        # Check if selected columns are present in the DataFrame
        if x_data not in self.data.columns or y_data not in self.data.columns:
            return

        fig, ax = plt.subplots(figsize=(5.2, 5))  # Set figure size to 520x500 pixels (each inch is 100 pixels)

        # Identify the first column
        first_column_name = self.data.columns[0]

        # Replace underscores with spaces in the first column's data
        self.data[first_column_name] = self.data[first_column_name].str.replace("_", " ")

        # Round function for better readability
        def round_numbers(val, digits=3):
            return round(val, digits)

        if plot_type == "Bar Graph":
            self.data.plot(kind="bar", x=x_data, y=y_data, ax=ax)
            # Add first column values as labels
            for i, txt in enumerate(self.data[first_column_name]):
                ax.text(
                    i,
                    round_numbers(self.data[y_data][i]),
                    txt,
                    ha="center",
                    va="bottom",
                )
        elif plot_type == "Scatter Plot":
            self.data.plot(kind="scatter", x=x_data, y=y_data, ax=ax)
            # Add first column values as labels
            for i, txt in enumerate(self.data[first_column_name]):
                ax.annotate(
                    txt,
                    (
                        round_numbers(self.data[x_data][i]),
                        round_numbers(self.data[y_data][i]),
                    ),
                )
        elif plot_type == "Line Plot":
            self.data.plot(kind="line", x=x_data, y=y_data, ax=ax)
            # Add first column values as labels
            for i, txt in enumerate(self.data[first_column_name]):
                ax.text(
                    i,
                    round_numbers(self.data[y_data][i]),
                    txt,
                    ha="center",
                    va="bottom",
                )
        elif plot_type == "Pie Plot":
            self.data.set_index(x_data).plot(
                kind="pie",
                y=y_data,
                labels=self.data[first_column_name],
                ax=ax,
                legend=False,
            )
        elif plot_type == "Violin Plot":
            if pd.api.types.is_numeric_dtype(self.data[x_data]):
                # Bin the data and use midpoints for readability
                self.data["x_binned"] = pd.cut(self.data[x_data], bins=10)
                self.data["x_binned_mid"] = self.data["x_binned"].apply(lambda x: x.mid).astype(float)
                self.data["x_binned_mid"] = self.data["x_binned_mid"].round(1).astype(str)  # Round to 1 decimal place
                sns.violinplot(x="x_binned_mid", y=y_data, data=self.data, ax=ax)
                ax.set_xlabel(x_data)  # Set the X-axis label
            else:
                sns.violinplot(x=x_data, y=y_data, data=self.data, ax=ax)
                ax.set_xlabel(x_data)  # Set the X-axis label

        plot_path = os.path.join("scripts", "results", f"{plot_type.lower().replace(' ', '_')}.png")
        os.makedirs("scripts/results", exist_ok=True)
        plt.savefig(plot_path)
        plt.close(fig)

        # Display plot in QLabel
        pixmap = QPixmap(plot_path)
        self.ClimaticChart_2.setPixmap(pixmap)
        self.tabWidget2.setCurrentIndex(2)

    def download_climate_plot(self):
        """
        Download the generated plot.

        This method allows the user to download the currently displayed plot.

        Returns:
            None
        """
        plot_type = self.PlotTypesCombobox.currentText()
        if plot_type == "":
            self.show_error_dialog("Please generate a plot first.")
            return

        plot_path = os.path.join("results", f"{plot_type.lower().replace(' ', '_')}.png")
        if not os.path.exists(plot_path):
            self.show_error_dialog("No plot found to download.")
            return

        # Prompt the user to select a location to save the plot
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Plot As",
            os.path.basename(plot_path),
            "PNG Files (*.png);;All Files (*)",
            options=options,
        )
        if file_path:
            if not file_path.lower().endswith(".png"):
                file_path += ".png"
            shutil.copy(plot_path, file_path)

    # Pas certain que c'est utiliser...
    def stop_thread(self):
        if self.worker:
            self.worker.stop()
        if self.thread and self.thread.isRunning():
            self.thread.quit()
            self.thread.wait()

    # Pas certain que c'est utiliser...
    def closeEvent(self, event):
        self.stop_thread()
        event.accept()

    def retrieve_data_names(self, data_list):
        """
        Retrieve data from a list, excluding the first element.

        Args:
            data_list (list): The list to retrieve data from.

        Returns:
            list: A list of data excluding the first element.
        """
        try:
            if not data_list:
                raise ValueError("The provided list is empty.")
            return data_list[1:]
        except ValueError as e:
            self.show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def load_csv_climate_file(self):
        def create_sleek_table(df):
            num_rows, num_columns = df.shape
            table_widget = QTableWidget(num_rows, num_columns)
            table_widget.setStyleSheet(
                """
                QTableWidget {
                    background-color: #f2f2f2;
                    border: 1px solid #ddd;
                }
                QHeaderView::section {
                    background-color: #4CAF50;
                    color: white;
                    font-weight: bold;
                    border: none;
                    padding: 5px;
                }
                QTableWidget::item {
                    border: none;
                    padding: 5px;
                }
            """
            )
            table_widget.horizontalHeader().setStretchLastSection(True)
            table_widget.verticalHeader().setVisible(False)
            table_widget.horizontalHeader().setVisible(True)
            table_widget.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)
            table_widget.horizontalHeader().setDefaultSectionSize(150)

            # Set headers
            for col in range(num_columns):
                item = QTableWidgetItem(df.columns[col])
                item.setTextAlignment(Qt.AlignCenter)
                table_widget.setHorizontalHeaderItem(col, item)

            # Fill the table with data
            for row in range(num_rows):
                for col in range(num_columns):
                    value = str(df.iloc[row, col])
                    if re.search("^[0-9\\-]*\\.[0-9]*", value) is not None:
                        value = str(round(Decimal(value), 3))
                    item = QTableWidgetItem(value)
                    item.setTextAlignment(Qt.AlignCenter)
                    table_widget.setItem(row, col, item)

            return table_widget

        try:
            options = QFileDialog.Options()
            options |= QFileDialog.ReadOnly
            fullFilePath, _ = QFileDialog.getOpenFileName(
                None,
                "Select CSV file",
                "./datasets",
                "Comma Separated Values (*.csv)",
                options=options,
            )

            if fullFilePath:
                update_yaml_param(Params, "scripts/utils/params.yaml", "file_name", fullFilePath)
                self.statisticsButtonPage2.setEnabled(True)
                df = pd.read_csv(fullFilePath)
                columns = df.columns.tolist()

                if len(columns) < 2:
                    raise ValueError("The CSV file must contain at least two columns for latitude and longitude.")

                latitude_col = columns[-1]
                longitude_col = columns[-2]

                lat = df[latitude_col].tolist()
                long = df[longitude_col].tolist()

                # Replace underscores with spaces in species names (first column)
                df[columns[0]] = df[columns[0]].str.replace("_", " ")

                self.species = df[columns[0]].tolist()
                self.factors = df.drop(columns=[longitude_col, latitude_col]).values.tolist()
                clim_data_names = self.retrieve_data_names(columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "names", columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "data_names", clim_data_names)

                self.textEditClimData.clear()
                sleek_table = create_sleek_table(df)

                # Add the table widget to the textEditClimData layout
                layout = QVBoxLayout(self.textEditClimData)
                layout.addWidget(sleek_table)

                self.climaticTrees = utils.climaticPipeline(df)
                self.tree_keys = list(self.climaticTrees.keys())
                self.total_trees = len(self.tree_keys)
                self.current_index = 0
                self.climaticTreeButtonPage2.setEnabled(True)
                self.tabWidget2.setCurrentIndex(1)
                self.populate_map(lat, long)

                if self.statisticsButtonPage1.isEnabled():
                    self.resultsButton.setEnabled(True)
        except FileNotFoundError as e:
            self.show_error_dialog(f"File Not Found Error: {e}")
        except pd.errors.EmptyDataError as e:
            self.show_error_dialog(f"Empty Data Error: {e}")
        except ValueError as e:
            self.show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def populate_map(self, lat, long):
        """
        Create and display a folium map with given latitude and longitude.

        This method generates a map centered on the mean latitude and longitude of the provided coordinates.
        It places markers on the map for each coordinate pair and displays the map in a QWebEngineView.

        Args:
            lat (list): List of latitudes.
            long (list): List of longitudes.
        """
        try:
            if not lat or not long or len(lat) != len(long):
                raise ValueError("Latitude and longitude lists must be non-empty and of the same length.")

            lat = [float(x) for x in lat]  # Convert all elements in lat to float
            long = [float(x) for x in long]  # Convert all elements in long to float

            mean_lat = sum(lat) / len(lat)
            mean_long = sum(long) / len(long)

            m = folium.Map(location=[mean_lat, mean_long], zoom_start=4, tiles="OpenStreetMap")  # Adjusted zoom level for better visibility

            # Add markers to the map
            for latitude, longitude in zip(lat, long):
                folium.Marker([latitude, longitude]).add_to(m)

            # Adjust the map to fit all markers
            m.fit_bounds([[min(lat), min(long)], [max(lat), max(long)]])

            data = io.BytesIO()
            m.save(data, close_file=False)

            web_view = QWebEngineView(self.graphicsViewClimData)  # Embed the map inside graphicsViewClimData
            web_view.setHtml(data.getvalue().decode())
            layout = QVBoxLayout(self.graphicsViewClimData)
            layout.addWidget(web_view)
            self.graphicsViewClimData.setLayout(layout)

        except ValueError as e:
            self.show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def show_home_section(self):
        """
        Display the home page of the application.

        This method sets the icons for the climatic data and genetic data buttons to their inactive states
        and displays the home page by setting the stacked widget's current index to 0.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.homeButton.setIcon(QIcon(":active/home.png"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.stackedWidget.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def show_genetic_section(self):
        """
        Display the genetic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the genetic data page
        by setting the stacked widget's current index to 1, and sets the tab widget's current index to 0.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.geneticDataButton.setIcon(QIcon(":active/genetic.svg"))
            self.homeButton.setIcon(QIcon(":other/home.svg"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.stackedWidget.setCurrentIndex(1)
            self.tabWidget.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def show_climate_section(self):
        """
        Display the climatic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the climatic data page
        by setting the stacked widget's current index to 2, and sets the tab widget's current index to 0.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":active/climaticData.png"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.homeButton.setIcon(QIcon(":other/home.svg"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.stackedWidget.setCurrentIndex(2)
            self.tabWidget2.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def show_results_section(self):
        """
        Display the results page of the application.

        This method sets the stacked widget's current index to 3 to display the results page.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.homeButton.setIcon(QIcon(":other/home.svg"))
            self.resultsButton.setIcon(QIcon(":active/result.svg"))
            self.stackedWidget.setCurrentIndex(3)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def show_filtered_results(self):
        """
        Show the results filtered with a metric threshold provided by the user.

        This method reads the data from a CSV file, processes it through the climatic pipeline, filters the results,
        and displays the filtered results in an HTML table format within a QTextBrowser widget.
        It handles exceptions related to missing sequence alignment.

        Raises:
            AttributeError: If the sequence alignment has not been performed before attempting to generate the tree.
        """
        self.stackedWidget.setCurrentIndex(3)

        df = pd.read_csv(Params.file_name)
        utils.filterResults(self.climaticTrees, self.genetics.geneticTreeDict, df)
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
        layout = QVBoxLayout(self.textEditResults)
        for i in reversed(range(layout.count())):
            widget_to_remove = layout.itemAt(i).widget()
            layout.removeWidget(widget_to_remove)
            widget_to_remove.setParent(None)

        # Add the QWebEngineView to the QTextBrowser
        layout.addWidget(web_engine_view)

    def toggle_dark_mode(self):
        
        """
        Toggle the application's dark mode setting.

        This method switches the application's theme between dark mode and light mode. It updates the isDarkMode attribute,
        applies the corresponding style, and changes the icon of the darkModeButton.

        Attributes:
            isDarkMode (bool): A flag indicating whether dark mode is currently enabled.

        Actions:
            If dark mode is enabled, apply the dark style and set the darkModeButton icon to the 'light' icon.
            If dark mode is disabled, apply the light style and set the darkModeButton icon to the 'dark' icon.
        """
        try:
            self.isDarkMode = not self.isDarkMode
            buttons_vertical = [
                self.fileBrowserButtonPage1,
                self.sequenceAlignmentButtonPage1,
                self.clearButtonPage1,
                self.statisticsButtonPage1,
                self.geneticTreeButtonPage1,
                self.fileBrowserButtonPage2,
                self.clearButtonPage2,
                self.climaticTreeButtonPage2,
                self.statisticsButtonPage2,
                self.settingsButtonPage3,
                self.settingsButtonPage4,
                self.submitButtonPage3,
                self.statisticsButtonPage3,
                self.submitButtonPage4,
                self.statisticsButtonPage4,
                self.clearButtonPage3,
                self.StartSequenceAlignmentButton,
                self.clearButtonPage4,
            ]
            buttons = [
                self.geneticDataButton,
                self.climaticDataButton,
                self.helpButton,
                self.homeButton,
                self.resultsButton,
            ]

            style_buttons(buttons, self.isDarkMode)
            style_buttons(buttons_vertical, self.isDarkMode)

            if self.isDarkMode:
                qtmodern.styles.dark(app)
                self.top_frame.setStyleSheet("background-color: #646464;")
                self.darkModeButton.setIcon(QIcon(":other/light.png"))  # Set the 'light' icon for dark mode
            else:
                qtmodern.styles.light(app)
                self.top_frame.setStyleSheet("background-color: rgb(222, 221, 218);")
                self.darkModeButton.setIcon(QIcon(":other/dark.png"))  # Set the 'dark' icon

            # Common settings for both modes
            self.darkModeButton.setCursor(Qt.PointingHandCursor)
            self.darkModeButton.setStyleSheet(get_button_style(self.isDarkMode))
            self.darkModeButton.setGraphicsEffect(create_shadow_effect(10, 140))

        except AttributeError as e:
            self.show_error_dialog(f"An error occurred while setting attributes: {e}", "Attribute Error")

        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")


    def clear_climmatic_data(self):
        """
        Clear the text fields related to climatic data.

        This method disables buttons related to climatic data and clears the necessary fields.
        """
        try:
            self.statisticsButtonPage2.setEnabled(False)
            self.climaticTreeButtonPage2.setEnabled(False)
            self.resultsButton.setEnabled(False)
            self.climaticTrees = None
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}", "Error")

    # elle se fait override par celle qui suit, regarder si elle est utile au final (penses pas)
    def clear_result(self):
        """
        Clear the text fields related to climatic data.

        This method clears the content of the text edit widgets for climatic data, climatic statistics, and climatic trees.
        It also clears the graphics view for climatic data and resets the current index of the climatic statistics and chart lists to 0.
        """
        try:
            self.textEditResults.clear()
            self.textEditClimStats.clear()
            self.textEditClimTrees.clear()
            self.graphicsViewClimData.clear()
            self.ClimStatsListCondition.setCurrentIndex(0)
            self.ClimStatsListChart.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}", "Error")

    def clear_result_stat(self):
        """
        Clear the statistics result lists.

        This method resets the current index of the results statistics condition list and the results statistics chart list to 0.
        """
        try:
            self.ResultsStatsListCondition.setCurrentIndex(0)
            self.ResultsStatsListChart.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}", "Error")

    def open_climatic_tree_preferences_window(self):
        """
        Open the preferences dialog and update the application settings based on user input.

        This method opens a PreferencesDialog, updates it with the current preferences, and applies the new preferences if the user accepts the changes.

        Returns:
            None
        """
        try:
            dialog = PreferencesDialog(self)
            dialog.update_preferences(self.preferences)
            if dialog.exec_() == QDialog.Accepted:
                self.preferences = dialog.get_preferences()
                self.apply_preferences()
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while opening the preferences dialog: {e}")

    def apply_preferences(self):
        """
        Apply the updated preferences to the current plot.

        This method re-displays the climatic tree based on the current preferences.

        Returns:
            None
        """
        try:
            self.show_climatic_tree(self.current_index)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while applying preferences: {e}")

    def display_climatic_trees(self):
        """
        Display the climatic trees in the application.

        This method populates the climaticTreescomboBox with the keys of the climatic trees,
        sets the total number of trees, initializes the current index, and displays the first tree.
        It also switches to the climatic tree tab.

        Returns:
            None
        """
        try:
            self.climaticTreescomboBox.clear()
            self.tree_keys = list(self.climaticTrees.keys())
            self.total_trees = len(self.tree_keys)
            self.current_index = 0
            self.climaticTreescomboBox.addItems(self.tree_keys)
            self.show_climatic_tree(self.current_index)
            self.tabWidget2.setCurrentIndex(3)
        except KeyError as e:
            self.show_error_dialog(
                f"An error occurred while accessing the climatic trees: {e}",
                "Key Error",
            )
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def show_selected_climatic_tree(self, index):
        """
        Display the selected climatic tree based on the provided index.

        Args:
            index (int): The index of the selected climatic tree in the combo box.

        Returns:
            None
        """
        try:
            if index >= 0:
                self.show_climatic_tree(index)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while displaying the selected climatic tree: {e}")

    def show_climatic_tree(self, index):
        """
        Display the climatic tree at the specified index.

        This method retrieves the preferences and displays the climatic tree based on the selected view type (network or tree).

        Args:
            index (int): The index of the tree to display.

        Returns:
            None
        """
        try:
            if 0 <= index < self.total_trees:
                self.current_index = index
                key = self.tree_keys[index]
                tree = self.climaticTrees[key]

                # Get preferences
                preferences = self.preferences
                label_color = preferences.get("label_color", "black")
                edge_color = preferences.get("edge_color", "blue")
                reticulation_color = preferences.get("reticulation_color", "red")
                layout = preferences.get("layout", "horizontal")
                proportional_edge_lengths = preferences.get("proportional_edge_lengths", False)
                label_internal_vertices = preferences.get("label_internal_vertices", False)
                use_leaf_names = preferences.get("use_leaf_names", True)
                show_branch_length = preferences.get("show_branch_length", False)
                view_type = preferences.get("view_type", "network")

                if view_type == "network":
                    self.render_network_view(
                        tree,
                        label_color,
                        edge_color,
                        reticulation_color,
                        layout,
                        proportional_edge_lengths,
                        label_internal_vertices,
                        use_leaf_names,
                        show_branch_length,
                    )
                else:
                    self.render_tree_view(
                        tree,
                        label_color,
                        edge_color,
                        reticulation_color,
                        layout,
                        proportional_edge_lengths,
                        label_internal_vertices,
                        use_leaf_names,
                        show_branch_length,
                    )
        except KeyError as e:
            self.show_error_dialog(
                f"An error occurred while accessing the climatic tree data: {e}",
                "Key Error",
            )
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def render_network_view(
        self,
        tree,
        label_color,
        edge_color,
        reticulation_color,
        layout,
        proportional_edge_lengths,
        label_internal_vertices,
        use_leaf_names,
        show_branch_length,
    ):
        """
        Render the network view of the phylogenetic tree.

        This method applies user preferences to the tree and renders it using Plotly for interactive visualization.

        Args:
            tree (Phylo.BaseTree.Tree): The phylogenetic tree to render.
            label_color (str): Color for labels.
            edge_color (str): Color for edges.
            reticulation_color (str): Color for reticulation edges.
            layout (str): Layout for the network visualization (e.g., 'horizontal').
            proportional_edge_lengths (bool): Whether to use proportional edge lengths.
            label_internal_vertices (bool): Whether to label internal vertices.
            use_leaf_names (bool): Whether to use leaf names.
            show_branch_length (bool): Whether to show branch lengths.

        Returns:
            None
        """
        try:
            # Use seaborn for color palettes
            sns.set_palette("husl")

            # Apply additional preferences to the tree

            # Proportional Edge Lengths
            if proportional_edge_lengths:
                for clade in tree.find_clades():
                    if clade.branch_length is None:
                        clade.branch_length = 0.1  # Set a default length if not specified

            # Label Internal Vertices
            if label_internal_vertices:
                for clade in tree.find_clades():
                    if not clade.is_terminal() and clade.name is None:
                        clade.name = "Internal Node"

            # Render the tree using Plotly for interactive visualization
            graph = Phylo.to_networkx(tree)
            pos = self.get_layout(graph, layout)
            edge_trace, edge_annotations = self.create_edge_trace(tree, pos, edge_color, show_branch_length)
            node_trace = self.create_node_trace(graph, pos, label_color, use_leaf_names)

            fig = go.Figure(data=[edge_trace, node_trace])
            fig.update_layout(
                showlegend=False,
                xaxis=dict(showgrid=False, zeroline=False, visible=False),
                yaxis=dict(showgrid=False, zeroline=False, visible=False),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                width=911,
                height=441,
            )

            # Add edge annotations for branch lengths
            for annotation in edge_annotations:
                fig.add_annotation(annotation)

            results_dir = "results"
            os.makedirs(results_dir, exist_ok=True)
            img_path = os.path.join(results_dir, f"{self.tree_keys[self.current_index]}.png")
            pio.write_image(fig, img_path, format="png")

            pixmap = QPixmap(img_path)
            self.climaticTreesLabel.clear()
            self.climaticTreesLabel.setPixmap(pixmap)
            self.climaticTreesLabel.adjustSize()
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while rendering the network view: {e}")

    def render_tree_view(
        self,
        tree,
        label_color,
        proportional_edge_lengths,
        label_internal_vertices,
        use_leaf_names,
        show_branch_length,
    ):
        """
        Render the tree view of the phylogenetic tree.

        This method applies user preferences to the tree and renders it using Matplotlib.

        Args:
            tree (Phylo.BaseTree.Tree): The phylogenetic tree to render.
            label_color (str): Color for labels.
            proportional_edge_lengths (bool): Whether to use proportional edge lengths.
            label_internal_vertices (bool): Whether to label internal vertices.
            use_leaf_names (bool): Whether to use leaf names.
            show_branch_length (bool): Whether to show branch lengths.

        Returns:
            None
        """
        try:
            fig = plt.figure(figsize=(9.11, 4.41))  # Limit size to 911x441 pixels
            ax = fig.add_subplot(1, 1, 1)

            # Apply additional preferences to the tree

            # Proportional Edge Lengths
            if proportional_edge_lengths:
                for clade in tree.find_clades():
                    if clade.branch_length is None:
                        clade.branch_length = 0.1  # Set a default length if not specified

            # Label Internal Vertices
            if label_internal_vertices:
                for clade in tree.find_clades():
                    if not clade.is_terminal() and clade.name is None:
                        clade.name = "Internal Node"

            # Draw the tree using Matplotlib
            def label_func(clade):
                label = clade.name if use_leaf_names and clade.is_terminal() else ""
                return f"{label}\n{clade.branch_length:.2f}" if show_branch_length and label else label

            Phylo.draw(
                tree,
                do_show=False,
                axes=ax,
                label_func=label_func,
                label_colors={clade: label_color for clade in tree.find_clades()},
            )

            ax.axis("off")  # Remove axes

            results_dir = "results"
            os.makedirs(results_dir, exist_ok=True)
            img_path = os.path.join(results_dir, f"{self.tree_keys[self.current_index]}.png")
            plt.savefig(img_path, format="png")
            plt.close(fig)

            pixmap = QPixmap(img_path)
            self.climaticTreesLabel.clear()
            self.climaticTreesLabel.setPixmap(pixmap)
            self.climaticTreesLabel.adjustSize()
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while rendering the tree view: {e}")

    def create_node_trace(self, graph, pos, label_color, use_leaf_names):
        """
        Create a Plotly node trace for the phylogenetic tree network visualization.

        Args:
            graph (networkx.Graph): The networkx graph of the phylogenetic tree.
            pos (dict): A dictionary of positions keyed by node.
            label_color (str): Color for the node labels.
            use_leaf_names (bool): Whether to use leaf names as labels.

        Returns:
            go.Scatter: A Plotly scatter trace representing the nodes.
        """
        try:
            node_trace = go.Scatter(
                x=[],
                y=[],
                text=[],
                mode="markers+text" if use_leaf_names else "markers",
                textposition="top center",
                hoverinfo="text",
                marker=dict(
                    showscale=False,  # Disable color scale
                    colorscale="Viridis",  # Use a valid Plotly colorscale
                    size=10,
                    line_width=2,
                    color=label_color,
                ),
            )

            for node in graph.nodes():
                x, y = pos[node]
                node_trace["x"] += (x,)
                node_trace["y"] += (y,)
                if use_leaf_names:
                    name = node.name if hasattr(node, "name") and node.name else ""
                    node_trace["text"] += (name,)

            return node_trace
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while creating the node trace: {e}")

    def get_layout(self, graph, layout):
        """
        Get the layout for the phylogenetic tree network visualization.

        Args:
            graph (networkx.Graph): The networkx graph of the phylogenetic tree.
            layout (str): The layout type (e.g., 'horizontal', 'vertical', 'radial', 'axial').

        Returns:
            dict: A dictionary of positions keyed by node.
        """
        try:
            if layout == "horizontal":
                return nx.spring_layout(graph, scale=2)
            elif layout == "vertical":
                return nx.spring_layout(graph, scale=2, iterations=50)
            elif layout == "radial":
                return nx.shell_layout(graph)
            elif layout == "axial":
                return nx.spiral_layout(graph)
            else:
                raise ValueError(f"Unknown layout type: {layout}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while getting the layout: {e}")

    def create_edge_trace(self, tree, pos, edge_color, show_branch_length):
        """
        Create a Plotly edge trace for the phylogenetic tree network visualization.

        Args:
            tree (Phylo.BaseTree.Tree): The phylogenetic tree.
            pos (dict): A dictionary of positions keyed by node.
            edge_color (str): Color for the edges.
            show_branch_length (bool): Whether to show branch lengths as annotations.

        Returns:
            tuple: A tuple containing a Plotly scatter trace for edges and a list of edge annotations.
        """
        try:
            edge_trace = go.Scatter(
                x=[],
                y=[],
                line=dict(width=2, color=edge_color),
                hoverinfo="none",
                mode="lines",
            )
            edge_annotations = []

            for clade in tree.find_clades(order="level"):
                if clade.is_terminal():
                    continue
                for child in clade.clades:
                    x0, y0 = pos[clade]
                    x1, y1 = pos[child]
                    edge_trace["x"] += (x0, x1, None)
                    edge_trace["y"] += (y0, y1, None)

                    if show_branch_length:
                        mid_x = (x0 + x1) / 2
                        mid_y = (y0 + y1) / 2
                        branch_length = child.branch_length
                        edge_annotations.append(
                            dict(
                                x=mid_x,
                                y=mid_y,
                                text=f"{branch_length:.2f}",
                                showarrow=False,
                                xanchor="center",
                                yanchor="middle",
                                font=dict(color=edge_color),
                            )
                        )
            return edge_trace, edge_annotations
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while creating the edge trace: {e}")

    def download_climatic_tree_graph(self):
        """
        Download the current displayed climatic tree graph as a PNG file.

        This method prompts the user to select a location to save the current graph,
        and saves the graph to the specified location.

        Returns:
            None
        """
        try:
            current_key = self.tree_keys[self.current_index]
            default_file_name = f"{current_key}.png"

            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
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
                results_dir = "results"
                os.makedirs(results_dir, exist_ok=True)
                img_path = os.path.join(results_dir, f"{self.tree_keys[self.current_index]}.png")
                with open(img_path, "rb") as file:
                    with open(file_path, "wb") as dest_file:
                        dest_file.write(file.read())
        except FileNotFoundError as e:
            self.show_error_dialog(f"The temporary image file was not found: {e}", "File Not Found")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while downloading the climatic tree graph: {e}")

    ################################################
    def display_phylogeographic_trees(self):
        self.stackedWidget.setCurrentIndex(4)
        self.results_dir = "scripts/results"
        file_path = os.path.join(self.results_dir, "geneticTrees.json")

        with open(file_path, "r") as file:
            self.phylo_json = json.load(file)

        self.tree_keys = list(self.phylo_json.keys())
        self.total_trees = len(self.tree_keys)

        self.current_index1 = 0
        self.phyloTreescomboBox.clear()

        formatted_tree_keys = [self.rename_tree_key(key) for key in self.tree_keys]
        self.phyloTreescomboBox.addItems(formatted_tree_keys)

        self.climatic_data = pd.read_csv(Params.file_name)
        self.criteriaComboBox.clear()
        self.criteriaComboBox.addItems(self.climatic_data.columns[1:])

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

            selected_column = self.criteriaComboBox.currentText()

            # Match the original labels from the CSV
            ordered_climatic_data = self.climatic_data.set_index("id").reindex(tip_labels).reset_index()
            bar_values = ordered_climatic_data[selected_column].values

            # Ensure bar_values have no negative values
            bar_values = np.clip(bar_values, a_min=0, a_max=None)

            # Create a canvas for the plot
            canvas = toyplot.Canvas(width=921, height=450)

            # Add tree to canvas
            ax0 = canvas.cartesian(bounds=(50, 300, 50, 400), ymin=0, ymax=tree.ntips, padding=15)

            # Replace underscores with spaces for display purposes only
            tree.draw(axes=ax0, tip_labels=[label.replace("_", " ") for label in tip_labels])
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
            self.PhyloTreeLabel.clear()
            self.PhyloTreeLabel.setPixmap(pixmap)
            self.PhyloTreeLabel.adjustSize()

    def save_tree_graph(self):
        current_key = self.tree_keys[self.current_index1]
        default_file_name = f"{current_key}.png"

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
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


if __name__ == "__main__":
    app = QtWidgets.QApplication([])

    qtmodern.styles.light(app)

    window = UiMainWindow()

    mw = qtmodern.windows.ModernWindow(window)

    screen_geometry = app.primaryScreen().availableGeometry()

    center_point = screen_geometry.center()
    x = center_point.x() - mw.width() // 2
    y = center_point.y() - mw.height() // 2

    mw.move(x, y)

    mw.show()

    sys.exit(app.exec_())
