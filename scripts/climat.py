import io
import os
import re
import shutil
from decimal import Decimal

import folium
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import seaborn as sns
from aphylogeo import utils
from aphylogeo.params import Params
from Bio import Phylo
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPixmap
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import QDialog, QFileDialog, QTableWidget, QTableWidgetItem, QVBoxLayout
from utils.ClimaticGraphSettings import ClimaticGraphSettings
from utils.ClimaticPreferencesDialog import ClimaticPreferencesDialog
from utils.error_dialog import show_error_dialog
from utils.my_dumper import update_yaml_param

try:
    ClimaticGraphSettings.load_from_file("./scripts/utils/ClimaticGraphSettings.yaml")
except FileNotFoundError:
    # Use default settings if the file is not found
    pass


class Climat:
    def __init__(self, main):
        self.main = main
        self.climaticTrees = dict()

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

        plot_path = os.path.join("scripts", "results", f"{plot_type.lower().replace(' ', '_')}.png")
        os.makedirs("scripts/results", exist_ok=True)
        plt.savefig(plot_path)
        plt.close(fig)

        # Display plot in QLabel
        pixmap = QPixmap(plot_path)
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
        if plot_type == "":
            show_error_dialog("Please generate a plot first.")
            return

        plot_path = os.path.join("results", f"{plot_type.lower().replace(' ', '_')}.png")
        if not os.path.exists(plot_path):
            show_error_dialog("No plot found to download.")
            return

        # Prompt the user to select a location to save the plot
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self.main,
            "Save Plot As",
            os.path.basename(plot_path),
            "PNG Files (*.png);;All Files (*)",
            options=options,
        )
        if file_path:
            if not file_path.lower().endswith(".png"):
                file_path += ".png"
            shutil.copy(plot_path, file_path)

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

                self.climaticTrees = utils.climaticPipeline(df)
                self.tree_keys = list(self.climaticTrees.keys())
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

            web_view = QWebEngineView(self.main.graphicsViewClimData)  # Embed the map inside graphicsViewClimData

            web_view.setHtml(data.getvalue().decode())
            layout = QVBoxLayout(self.main.graphicsViewClimData)
            layout.addWidget(web_view)
            self.main.graphicsViewClimData.setLayout(layout)

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
            self.climaticTrees = None
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}", "Error")

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
            show_error_dialog(f"An unexpected error occurred while applying preferences: {e}")

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
            self.main.climaticTreescomboBox.clear()
            self.tree_keys = list(self.climaticTrees.keys())
            self.total_trees = len(self.tree_keys)
            self.current_index = 0
            self.main.climaticTreescomboBox.addItems(self.tree_keys)
            self.show_climatic_tree(self.current_index)
            self.main.tabWidget2.setCurrentIndex(3)
        except KeyError as e:
            show_error_dialog(
                f"An error occurred while accessing the climatic trees: {e}",
                "Key Error",
            )
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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
            show_error_dialog(f"An unexpected error occurred while displaying the selected climatic tree: {e}")

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

                label_color = ClimaticGraphSettings.label_color
                edge_color = ClimaticGraphSettings.edge_color
                reticulation_color = ClimaticGraphSettings.reticulation_color
                layout = ClimaticGraphSettings.layout
                proportional_edge_lengths = ClimaticGraphSettings.proportional_edge_lengths
                label_internal_vertices = ClimaticGraphSettings.label_internal_vertices
                use_leaf_names = ClimaticGraphSettings.use_leaf_names
                show_branch_length = ClimaticGraphSettings.show_branch_length
                view_type = ClimaticGraphSettings.view_type
                if view_type == "network":
                    self.render_network_view(
                        tree,
                        label_color,
                        edge_color,
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
                        proportional_edge_lengths,
                        label_internal_vertices,
                        use_leaf_names,
                        show_branch_length,
                    )
        except KeyError as e:
            show_error_dialog(
                f"An error occurred while accessing the climatic tree data: {e}",
                "Key Error",
            )
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def render_network_view(
        self,
        tree,
        label_color,
        edge_color,
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
            layout (str): Layout for the network visualization (e.g., 'horizontal').
            proportional_edge_lengths (bool): Whether to use proportional edge lengths.
            label_internal_vertices (bool): Whether to label internal vertices.
            use_leaf_names (bool): Whether to use leaf names.
            show_branch_length (bool): Whether to show branch lengths.

        Returns:
            None
        """
        try:
            sns.set_palette("husl")

            if proportional_edge_lengths:
                for clade in tree.find_clades():
                    if clade.branch_length is None:
                        clade.branch_length = 0.1  # Set a default length if not specified

            if label_internal_vertices:
                for clade in tree.find_clades():
                    if not clade.is_terminal() and clade.name is None:
                        clade.name = "Internal Node"

            graph = Phylo.to_networkx(tree)
            pos = self.get_layout(graph, layout)
            edge_trace_result = self.create_edge_trace(tree, pos, edge_color, show_branch_length)
            if edge_trace_result is None:
                raise TypeError("The edge_trace is None")
            edge_trace, edge_annotations = edge_trace_result
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
            self.main.climaticTreesLabel.clear()
            self.main.climaticTreesLabel.setPixmap(pixmap)
            self.main.climaticTreesLabel.adjustSize()
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while rendering the network view: {e}")

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
                label_colors=dict.fromkeys(tree.find_clades(), label_color),
            )

            ax.axis("off")  # Remove axes

            results_dir = "results"
            os.makedirs(results_dir, exist_ok=True)
            img_path = os.path.join(results_dir, f"{self.tree_keys[self.current_index]}.png")
            plt.savefig(img_path, format="png")
            plt.close(fig)

            pixmap = QPixmap(img_path)
            self.main.climaticTreesLabel.clear()
            self.main.climaticTreesLabel.setPixmap(pixmap)
            self.main.climaticTreesLabel.adjustSize()
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while rendering the tree view: {e}")

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
            show_error_dialog(f"An unexpected error occurred while creating the node trace: {e}")

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
            show_error_dialog(f"An unexpected error occurred while getting the layout: {e}")

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
            show_error_dialog(f"An unexpected error occurred while creating the edge trace: {e}")

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

            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getSaveFileName(
                self.main,
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
            show_error_dialog(f"The temporary image file was not found: {e}", "File Not Found")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while downloading the climatic tree graph: {e}")

    def open_climatic_tree_preferences_window(self):
        """
        Open the preferences dialog and update the application settings based on user input.

        This method opens a PreferencesDialog, updates it with the current preferences, and applies the new preferences if the user accepts the changes.

        Returns:
            None
        """
        dialog = ClimaticPreferencesDialog()
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.apply_preferences()
