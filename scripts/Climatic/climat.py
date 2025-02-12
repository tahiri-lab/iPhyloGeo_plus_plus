import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from aphylogeo import utils
from aphylogeo.params import Params

from PyQt6.QtWidgets import QFileDialog

from utils.error_dialog import show_error_dialog
from utils.my_dumper import update_yaml_param
from utils.download_file import download_file_local, download_file_temporary_PLT
from event_connector import blocked_signals

from Climatic.climat_tree import ClimaticTree
from Climatic.climat_data import get_folium_data, create_sleek_table


class Climat:
    def __init__(self, main):
        self.main = main
        self.climaticTree = ClimaticTree(main)
        self.sleek_table = None

    def load_climate_statistics(self):
        """
        Load climate data from a CSV file, update the UI elements with the column names, and switch to the appropriate tab.

        This method reads the climate data from a CSV file (excluding the first column), updates combo boxes
        for X and Y axis settings with the column names, and switches to the climate data tab.

        Returns:
            None
        """
        self.data = pd.read_csv(Params.file_name)
        columns = self.data.columns.tolist()

        if columns:
            columns.pop(0)

        with blocked_signals(self.main.ClimaticChartSettingsAxisX, self.main.ClimaticChartSettingsAxisY):
            self.main.ClimaticChartSettingsAxisX.clear()
            self.main.ClimaticChartSettingsAxisX.addItems(columns)
            self.main.ClimaticChartSettingsAxisY.clear()
            self.main.ClimaticChartSettingsAxisY.addItems(columns)
            self.main.ClimaticChartSettingsAxisY.setCurrentIndex(1)
            self.main.tabWidget2.setCurrentIndex(2)
            
        self.generate_climate_graph()
        
    
    def generate_climate_graph(self):
        """
        Generate and display a graph based on the selected X and Y axis data and the chosen plot type.

        This method reads the selected data columns and plot type from the UI, generates the corresponding graph,
        and displays it in the specified QLabel widget.

        Returns:
            None
        """
        x_data = self.main.ClimaticChartSettingsAxisX.currentText()
        y_data = self.main.ClimaticChartSettingsAxisY.currentText()
        plot_type = self.main.PlotTypesCombobox.currentText()

        fig, ax = plt.subplots(figsize=(5.2, 5))  # Set figure size to 520x500 pixels (each inch is 100 pixels)

        # Identify the first column
        first_column_name = self.data.columns[0]

        # Replace underscores with spaces in the first column's data
        self.data[first_column_name] = self.data[first_column_name].str.replace("_", " ")
        
        match plot_type:
            case "Bar Graph":
                self.generate_bar_graph(x_data, y_data, ax, first_column_name)
            case "Scatter Plot":
                self.generate_scatter_plot(x_data, y_data, ax, first_column_name)
            case "Line Plot":
                self.generate_line_plot(x_data, y_data, ax, first_column_name)
            case "Pie Plot":
                self.generate_pie_plot(x_data, y_data, ax, first_column_name)
            case "Violin Plot":
                self.generate_violin_plot(x_data, y_data, ax)

        pixmap = download_file_temporary_PLT(plot_type, fig)

        self.main.ClimaticChart_2.setPixmap(pixmap)
        
        
    def generate_basic_graph(self, x_data, y_data, ax, first_column_name, kind):
            self.data.plot(kind=kind, x=x_data, y=y_data, ax=ax)
            # Add first column values as labels
            for i, txt in enumerate(self.data[first_column_name]):
                ax.text(
                    i,
                    round_numbers(self.data[y_data][i]),
                    txt,
                    ha="center",
                    va="bottom",
                )

    def generate_bar_graph(self, x_data, y_data, ax, first_column_name):
            self.generate_basic_graph(x_data, y_data, ax, first_column_name, "bar")
                
    def generate_scatter_plot(self, x_data, y_data, ax, first_column_name):
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
                
    def generate_line_plot(self, x_data, y_data, ax, first_column_name):
            self.generate_basic_graph(x_data, y_data, ax, first_column_name, "line")
                
    def generate_pie_plot(self, x_data, y_data, ax, first_column_name):
            self.data.set_index(x_data).plot(
                kind="pie",
                y=y_data,
                labels=self.data[first_column_name],
                ax=ax,
                legend=False,
            )
            
    def generate_violin_plot(self, x_data, y_data, ax):
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
                
                
    
    def download_climate_plot(self):
        """
        Download the generated plot.

        This method allows the user to download the currently displayed plot.

        Returns:
            None
        """
        plot_type = self.main.PlotTypesCombobox.currentText()
        download_file_local(plot_type, self.main)


    def load_csv_climate_file(self):
        try:
            options = QFileDialog.Option.ReadOnly
            fullFilePath, _ = QFileDialog.getOpenFileName(
                None,
                "Select CSV file",
                "./datasets",
                "Comma Separated Values (*.csv)",
                options=options,
            )

            if fullFilePath:
                update_yaml_param(Params, "scripts/utils/params.yaml", "file_name", fullFilePath)
                df = pd.read_csv(fullFilePath)
                columns = df.columns.tolist()

                if len(columns) < 2:
                    raise ValueError("The CSV file must contain at least two columns for latitude and longitude.")

                latitude_col = columns[-1]
                longitude_col = columns[-2]

                lat = df[latitude_col].tolist()
                long = df[longitude_col].tolist()

                df[columns[0]] = df[columns[0]].str.replace("_", " ")

                self.species = df[columns[0]].tolist()
                self.factors = df.drop(columns=[longitude_col, latitude_col]).values.tolist()
                clim_data_names = columns[1:]
                update_yaml_param(Params, "scripts/utils/params.yaml", "names", columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "data_names", clim_data_names)

                self.main.textEditClimData.clear()
                if self.sleek_table is not None:
                    self.main.climatTableLayout.removeWidget(self.sleek_table)
                    self.sleek_table.deleteLater()  
                    
                self.sleek_table = create_sleek_table(df)

                self.main.climatTableLayout.addWidget(self.sleek_table)

                self.climaticTree.climaticTrees = utils.climaticPipeline(df)
                self.tree_keys = list(self.climaticTree.climaticTrees.keys())
                self.total_trees = len(self.tree_keys)
                self.current_index = 0
                self.main.climaticTreeButtonPage2.setEnabled(True)
                self.main.statisticsButtonPage2.setEnabled(True)
                self.main.tabWidget2.setCurrentIndex(1)
                self.main.mapView.setHtml(get_folium_data(lat, long))

                if self.main.statisticsButtonPage1.isEnabled():
                    self.main.resultsButton.setEnabled(True)
        except FileNotFoundError as e:
            show_error_dialog(f"File Not Found Error: {e}")
        except pd.errors.EmptyDataError as e:
            show_error_dialog(f"Empty Data Error: {e}")
        except ValueError as e:
            show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def clear_climmatic_data(self):
        """
        Clear the text fields related to climatic data.

        This method disables buttons related to climatic data and clears the necessary fields.
        """
        try:
            with blocked_signals(self.main.ClimaticChartSettingsAxisX, self.main.ClimaticChartSettingsAxisY):
                self.main.ClimaticChartSettingsAxisX.clear()
                self.main.ClimaticChartSettingsAxisY.clear()
                
            self.climaticTree.climaticTrees.clear()
            self.main.ClimaticChart_2.clear()
            self.main.climaticTreesLabel.clear()
            self.main.climatTableLayout.removeWidget(self.sleek_table)
            self.sleek_table.deleteLater()
            self.sleek_table = None          
            self.main.mapView.setHtml("")          
            
            self.main.statisticsButtonPage2.setEnabled(False)
            self.main.climaticTreeButtonPage2.setEnabled(False)
            self.main.resultsButton.setEnabled(False)
            self.main.tabWidget2.setCurrentIndex(0)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}", "Error")
            
            
def round_numbers(val, digits=3):
    return round(val, digits)
