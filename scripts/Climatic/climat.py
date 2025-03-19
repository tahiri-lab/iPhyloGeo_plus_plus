import plotly.express as px
import pandas as pd
import numpy as np
from aphylogeo import utils
from aphylogeo.params import Params
from Climatic.climat_data import get_folium_data
from Climatic.climat_tree import ClimaticTree
from event_connector import blocked_signals
from PyQt6.QtWidgets import QFileDialog
from utils.custom_table import create_sleek_table
from utils.download_file import download_local_with_fig
from utils.error_dialog import show_error_dialog
from utils.my_dumper import update_yaml_param


class Climat:
    def __init__(self, main):
        self.main = main
        self.climaticTree = ClimaticTree(main)
        self.sleek_table = None
        self.customConfig ={
                'modeBarButtonsToRemove': ['toImage', 'autoScale', 'select2d', 'lasso2d'],
                'displaylogo' : False
            }

    def load_climate_statistics(self):
        """
        Load climate data from a CSV file, update the UI elements with the column names, and switch to the appropriate tab.

        This method reads the climate data from a CSV file (excluding the first column), updates combo boxes
        for X and Y axis settings with the column names, and switches to the climate data tab.

        Returns:
            None
        """        
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
        
        self.main.ClimaticChartSettingsAxisX.setEnabled(True)
        self.main.ClimaticChartSettingsAxisY.setEnabled(True)
        
        match plot_type:
            case "Scatter Plot":
                self.generate_scatter_plot(x_data, y_data)
            case "Line Plot":
                self.generate_line_plot(x_data, y_data)
            case "Bar Graph":
                self.generate_bar_graph(x_data, y_data)
            case "Violin Plot":
                self.generate_violin_plot(x_data, y_data)
            case "Pie Plot":
                self.generate_pie_plot(x_data, y_data)
            case "Correlation":
                self.generate_correlation_plot()
                self.main.ClimaticChartSettingsAxisX.setEnabled(False)
                self.main.ClimaticChartSettingsAxisY.setEnabled(False)
                  
        self.main.climatGraphView.setHtml(self.fig.to_html(include_plotlyjs='cdn', config= self.customConfig))

    def generate_bar_graph(self, x_data, y_data):
        self.fig = px.bar(
        data_frame = self.data,
        x=x_data,
        y=y_data,
        hover_name= self.first_column_name
        )

    def generate_scatter_plot(self, x_data, y_data):
        self.fig = px.scatter(
            data_frame = self.data,
            x = x_data,
            y = y_data,
            hover_name = self.first_column_name
        )

    def generate_line_plot(self, x_data, y_data):
        self.fig = px.line(
            data_frame = self.data,
            x = x_data,
            y = y_data,
            hover_data=[self.first_column_name],
            markers=True
        )

    def generate_pie_plot(self, x_data, y_data):
        self.fig = px.pie(
        data_frame = self.data,
        names=self.data[self.first_column_name],
        values=y_data,
        hover_data=[x_data]
        )

    def generate_violin_plot(self, x_data, y_data):   
        self.fig = px.violin(
            data_frame = self.data,
            x = x_data,
            y = y_data,
            hover_name = self.first_column_name
        )
           
    def generate_correlation_plot(self):
        matrix_correlation = self.data.iloc[:, 1:].corr(method="spearman").round(2)
        self.fig = px.imshow(
            matrix_correlation,
            aspect="auto"
        )
         
    def download_climate_plot(self):
        """
        Download the generated plot.

        This method allows the user to download the currently displayed plot.

        Returns:
            None
        """
        download_local_with_fig(self.fig)

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
                
                self.data = pd.read_csv(fullFilePath)
                columns = self.data.columns.tolist()
                self.first_column_name = self.data.columns[0]
                self.data[columns[0]] = self.data[columns[0]].str.replace("_", " ")
                
                latitude_col = columns[-1]
                longitude_col = columns[-2]
                lat = self.data[latitude_col].tolist()
                long = self.data[longitude_col].tolist()
                
                self.filter_data()
                columns = self.data.columns.tolist()
                
                if len(columns) < 3:
                    raise ValueError("The data has not enough columns. Either change the csv file or adjust the filters.")

                clim_data_names = columns[1:]
                update_yaml_param(Params, "scripts/utils/params.yaml", "names", columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "data_names", clim_data_names)

                self.main.textEditClimData.clear()
                if self.sleek_table is not None:
                    self.main.climatTableLayout.removeWidget(self.sleek_table)
                    self.sleek_table.deleteLater()

                self.sleek_table = create_sleek_table(self.data, True)

                self.main.climatTableLayout.addWidget(self.sleek_table)

                self.climaticTree.climaticTrees = utils.climaticPipeline(self.data)
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
            
    def filter_data(self):
        variance = self.data.var(numeric_only=True)
        numeric_to_keep = set(variance[variance > self.main.min_variance_climat.value()].index.tolist())
        
        final_columns = [col for col in self.data.columns if (col not in variance.index or col in numeric_to_keep)]
        
        self.data = self.data[final_columns]
        
        cols = list(self.data.columns[1:])
        data_processed = self.data.copy()
        threshold = self.main.max_correlation_climat.value()
        
        while True:
            corr_matrix = data_processed[cols].corr(method="spearman").round(2).abs()
            upper_tri = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
            
            max_corr = upper_tri.max().max()
            if max_corr < threshold:
                break
            
            col_correlation_scores = {}
            for col in upper_tri.columns:
                high_corr = upper_tri[col][upper_tri[col] >= threshold]
                if not high_corr.empty:
                    col_correlation_scores[col] = high_corr.mean()
            
            if col_correlation_scores:
                column_to_drop = max(col_correlation_scores, key=col_correlation_scores.get)
                cols.remove(column_to_drop)
            
        self.data = data_processed[[self.data.columns[0]] + cols]   

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
