import json
import os
import sys


import numpy as np
import pandas as pd
import qtmodern.styles
import qtmodern.windows
import toyplot.png
import toytree
from aphylogeo import utils
from aphylogeo.params import Params
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QIcon, QPixmap
from PyQt6.QtWidgets import QFileDialog, QGraphicsDropShadowEffect, QVBoxLayout
from PyQt6.QtWebEngineWidgets import QWebEngineView
from Qt import main
from ui_helpers import create_shadow_effect, get_button_style, style_buttons
from utils import resources_rc  # noqa: F401  # Import the compiled resource module for resolving image resource path
from utils.result_settings_dialog import ResultSettingsDialog
from utils.error_dialog import show_error_dialog
from navigation import Navigation
from climat import Climat
from genetics import Genetics
from event_connector import connect_decorated_methods, QtEvents, connect_event

try:
    Params.load_from_file("./scripts/utils/params.yaml")
except FileNotFoundError:
    Params.validate_and_set_params(Params.PARAMETER_KEYS)
    

window_size = 50
starting_position = 1


class UiMainWindow(main.Ui_MainWindow, QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.genetics = Genetics(self)
        self.climat = Climat(self)
        self.navigation = Navigation(self)
        self.setupUi(self)
        self.setup_ui()
        connect_decorated_methods(self)

    @connect_event("settingsButtonPage3", QtEvents.clicked)
    def open_result_settings_window(self):
        """
        Initialize and display the parameters window.
        This method creates a new QMainWindow instance, sets up its UI using the UiDialog class, and displays the window.
        """

        dialog = ResultSettingsDialog()
        dialog.exec()

    def setup_ui(self):
        """
        Setup the UI components and initialize the main window.

        This method connects various UI buttons to their corresponding event handlers,
        sets up styles and effects for UI elements, and initializes the state of the application.
        """
        try:
            self.tree_keys = []
            self.total_trees = 0
            self.current_index = 0
            self.current_index1 = 0
            self.setObjectName("MainWindow")
            self.window_size_spinbox_2.setRange(1, 1000)
            self.starting_position_spinbox_2.setRange(1, 1000)
            self.isDarkMode = False  # Keep track of the state
            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
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
                button.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
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
                buttonV.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
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

            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
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
            show_error_dialog(f"An error occurred while setting up the UI: {e}", "Attribute Error")
        except Exception as e:
            show_error_dialog(
                f"An unexpected error occurred: {e}",
                "Unexpected Error",
            )

    @connect_event("clearButtonPage3", QtEvents.clicked)
    def clear_results(self):
        self.textEditResults.clear()


    #GENETICS
    
    @connect_event("geneticTreeButtonPage1", QtEvents.clicked)    
    def display_newick_trees_click(self):
        self.genetics.display_newick_trees()
        
    @connect_event(["starting_position_spinbox_2", "window_size_spinbox_2"], QtEvents.valueChanged)    
    def update_plot_click(self):
        self.genetics.update_plot()
        
    @connect_event("statisticsButtonPage1", QtEvents.clicked)        
    def initialize_species_list_click(self):
        self.genetics.initialize_species_list()
        
    @connect_event("downloadSimilarityButton", QtEvents.clicked)        
    def download_similarity_plot_chart_click(self):
        self.genetics.download_similarity_plot_chart()
    
    @connect_event("fileBrowserButtonPage1", QtEvents.clicked)   
    def select_fasta_file_click(self):
        self.genetics.select_fasta_file()
        
    @connect_event("StartSequenceAlignmentButton", QtEvents.clicked)        
    def start_alignment_analysis_click(self):
        self.genetics.start_alignment_analysis()
        
    @connect_event("sequenceAlignmentButtonPage1", QtEvents.clicked)          
    def show_sequence_alignment_page_click(self):
        self.genetics.show_sequence_alignment_page()
        
    @connect_event("clearButtonPage1", QtEvents.clicked)            
    def clear_genetic_data_click(self):
        self.genetics.clear_genetic_data()
        
    @connect_event("geneticTreescomboBox", QtEvents.currentIndexChanged) 
    def show_tree_click(self, index):
        self.genetics.show_tree(index)
        
    @connect_event("downloadGraphButton", QtEvents.clicked)           
    def download_genetic_tree_graph_click(self):
        self.genetics.download_genetic_tree_graph()
        
    @connect_event("geneticSettingsButton", QtEvents.clicked)     
    def open_genetic_settings_window_click(self):
        self.genetics.open_genetic_settings_window()
    
    
    
    
    #CLIMATE
    
    @connect_event("statisticsButtonPage2", QtEvents.clicked)      
    def load_climate_statistics_click(self):
        self.climat.load_climate_statistics()
        
    @connect_event(["ClimaticChartSettingsAxisX", "ClimaticChartSettingsAxisY", "PlotTypesCombobox"], QtEvents.currentIndexChanged)            
    def generate_climate_graph_click(self):
        self.climat.generate_climate_graph()  
          
    @connect_event("climatePlotDownloadButton", QtEvents.clicked)          
    def download_climate_plot_click(self):
        self.climat.download_climate_plot() 
        
    @connect_event("fileBrowserButtonPage2", QtEvents.clicked)        
    def load_csv_climate_file_click(self):
        self.climat.load_csv_climate_file()
        
    @connect_event("clearButtonPage2", QtEvents.clicked)        
    def clear_climmatic_data_click(self):
        self.climat.clear_climmatic_data()
    
    @connect_event("climaticTreeButtonPage2", QtEvents.clicked)      
    def display_climatic_trees_click(self):
        self.climat.display_climatic_trees()
        
    @connect_event("downloadGraphButton2", QtEvents.clicked)             
    def show_selected_climatic_tree_click(self):
        self.climat.show_selected_climatic_tree()
        
    @connect_event("preferencesButton", QtEvents.clicked)
    def open_climatic_tree_preferences_window_click(self):
        self.climat.open_climatic_tree_preferences_window()
    
    
    #NAVIGATION
    
    @connect_event("homeButton", QtEvents.clicked)      
    def show_home_section_click(self):
        self.navigation.show_home_section()
        
    @connect_event("geneticDataButton", QtEvents.clicked)        
    def show_genetic_section_click(self):
        self.navigation.show_genetic_section()
        
    @connect_event("climaticDataButton", QtEvents.clicked)        
    def show_climate_section_click(self):
        self.navigation.show_climate_section()
        
    @connect_event("resultsButton", QtEvents.clicked)        
    def show_results_section_click(self):
        self.navigation.show_results_section()
        
    @connect_event("helpButton", QtEvents.clicked)         
    def open_help_window_click(self):
        self.navigation.open_help_window()
    
    ################################


    def stop_thread(self):
        self.genetics.stopWorker()
        if self.thread and self.thread.isRunning():
            self.thread.quit()
            self.thread.wait()

    def closeEvent(self, event):
        self.stop_thread()
        event.accept()

    @connect_event("submitButtonPage3", QtEvents.clicked)
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
        utils.filterResults(self.climat.climaticTrees, self.genetics.geneticTreeDict, df)
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

    @connect_event("darkModeButton", QtEvents.clicked)
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
            show_error_dialog(f"An error occurred while setting attributes: {e}", "Attribute Error")

        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    # elle se fait override par celle qui suit, regarder si elle est utile au final (penses pas)
    @connect_event("clearButtonPage4", QtEvents.clicked)
    def clear_result(self):
        """
        Clear the text fields related to climatic data.

        This method clears the content of the text edit widgets for climatic data, climatic statistics, and climatic trees.
        It also clears the graphics view for climatic data and resets the current index of the climatic statistics and chart lists to 0.
        """
        try:
            self.textEditResults.clear()
            self.clearButtonPage3.clear()
            self.textEditClimTrees.clear()
            self.graphicsViewClimData.clear()
            self.ClimStatsListCondition.setCurrentIndex(0)
            self.ClimStatsListChart.setCurrentIndex(0)
            self.ResultsStatsListCondition.setCurrentIndex(0)
            self.ResultsStatsListChart.setCurrentIndex(0)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}", "Error")


   
    ################################################
    @connect_event(["statisticsButtonPage3", "statisticsButtonPage4"], QtEvents.clicked)
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

    @connect_event(["phyloTreescomboBox", "criteriaComboBox"], QtEvents.currentIndexChanged)
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

    @connect_event("downloadResultsPlotButton", QtEvents.clicked)
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

    sys.exit(app.exec())
