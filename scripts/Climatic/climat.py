import io
import folium

import pandas as pd
import numpy as np
from aphylogeo import utils
from aphylogeo.params import Params
from Climatic.climat_tree import ClimaticTree
from Climatic.climat_statistics import ClimaticStatistics
from event_connector import blocked_signals
from PyQt6.QtWidgets import QFileDialog
from utils.custom_table import create_sleek_table
from utils.error_dialog import show_error_dialog
from utils.my_dumper import update_yaml_param


class Climat:
    def __init__(self, main):
        self.main = main
        self.climaticTree = ClimaticTree(main)
        self.climaticStatistics = ClimaticStatistics(main)
        self.sleek_table = None

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
                self.data[columns[0]] = self.data[columns[0]].str.replace("_", " ")
                
                self.setup_climatic_data_page(columns, True)
                
                columns = self.data.columns.tolist()
                
                if len(columns) < 3:
                    raise ValueError("The data has not enough columns. Either change the csv file or adjust the filters.")

                clim_data_names = columns[1:]
                update_yaml_param(Params, "scripts/utils/params.yaml", "names", columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "data_names", clim_data_names)

                self.climaticTree.climaticTrees = utils.climaticPipeline(self.data)
                self.main.climaticTreeButtonPage2.setEnabled(True)
                self.main.statisticsButtonPage2.setEnabled(True)
                self.main.tabWidget2.setCurrentIndex(1)

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
        self.climaticStatistics.data = self.data
        

    def setup_climatic_data_page(self, columns, filter = True):
        
        self.main.textEditClimData.clear()
        if self.sleek_table is not None:
            self.main.climatTableLayout.removeWidget(self.sleek_table)
            self.sleek_table.deleteLater()
            
        latitude_col = columns[-1]
        longitude_col = columns[-2]
        lat = self.data[latitude_col].tolist()
        long = self.data[longitude_col].tolist()
        
        if filter:
            self.filter_data()

        self.sleek_table = create_sleek_table(self.data, True)

        self.main.climatTableLayout.addWidget(self.sleek_table)
        
        self.main.mapView.setHtml(get_folium_data(lat, long))
        
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
            
def get_folium_data(lat, long):
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

        return data.getvalue().decode()

    except ValueError as e:
        show_error_dialog(f"Value Error: {e}")
    except Exception as e:
        show_error_dialog(f"An unexpected error occurred: {e}")