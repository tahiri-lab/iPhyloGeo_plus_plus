import io
import re
import folium

from decimal import Decimal
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from aphylogeo import utils
from aphylogeo.params import Params

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QFileDialog, QTableWidget, QTableWidgetItem, QVBoxLayout

from utils.error_dialog import show_error_dialog
from utils.my_dumper import update_yaml_param
from utils.download_file import download_file_local, download_file_temporary_PLT

from Climatic.climat_tree import ClimaticTree


class Climat:
    def __init__(self, main):
        self.main = main
        self.climaticTree = ClimaticTree(main)

    def load_climate_statistics(self):
        """
        Load climate data from a CSV file, update the UI elements with the column names, and switch to the appropriate tab.

        This method reads the climate data from a CSV file (excluding the first column), updates combo boxes
        for X and Y axis settings with the column names, and switches to the climate data tab.

        Returns:
            None
        """
        self.data = pd.read_csv(Params.file_name)
        self.columns = self.data.columns.tolist()

        if self.columns:
            self.columns.pop(0)

        self.main.ClimaticChartSettingsAxisX.clear()
        self.main.ClimaticChartSettingsAxisX.addItems(self.columns)
        self.main.ClimaticChartSettingsAxisY.clear()
        self.main.ClimaticChartSettingsAxisY.addItems(self.columns)
        self.main.tabWidget2.setCurrentIndex(2)

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

        pixmap = download_file_temporary_PLT(plot_type, fig)

        self.main.ClimaticChart_2.setPixmap(pixmap)
        self.main.tabWidget2.setCurrentIndex(2)

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
            if horizontal_header := table_widget.horizontalHeader():
                horizontal_header.setStretchLastSection(True)
                horizontal_header.setVisible(True)
                horizontal_header.setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)
                horizontal_header.setDefaultSectionSize(150)

            if vertical_header := table_widget.verticalHeader():
                vertical_header.setVisible(False)

            for col in range(num_columns):
                item = QTableWidgetItem(df.columns[col])
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                table_widget.setHorizontalHeaderItem(col, item)

            for row in range(num_rows):
                for col in range(num_columns):
                    value = str(df.iloc[row, col])
                    if re.search("^[0-9\\-]*\\.[0-9]*", value) is not None:
                        value = str(round(Decimal(value), 3))
                    item = QTableWidgetItem(value)
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    table_widget.setItem(row, col, item)

            return table_widget

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
                self.main.statisticsButtonPage2.setEnabled(True)
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
                clim_data_names = self.retrieve_data_names(columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "names", columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "data_names", clim_data_names)

                self.main.textEditClimData.clear()
                sleek_table = create_sleek_table(df)

                # Add the table widget to the textEditClimData layout
                layout = QVBoxLayout(self.main.textEditClimData)
                layout.addWidget(sleek_table)

                self.climaticTree.climaticTrees = utils.climaticPipeline(df)
                self.tree_keys = list(self.climaticTree.climaticTrees.keys())
                self.total_trees = len(self.tree_keys)
                self.current_index = 0
                self.main.climaticTreeButtonPage2.setEnabled(True)
                self.main.tabWidget2.setCurrentIndex(1)
                self.populate_map(lat, long)

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
            show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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

            self.main.mapView.setHtml(data.getvalue().decode())

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
            self.main.statisticsButtonPage2.setEnabled(False)
            self.main.climaticTreeButtonPage2.setEnabled(False)
            self.main.resultsButton.setEnabled(False)
            self.climaticTree.climaticTrees.clear()
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}", "Error")
