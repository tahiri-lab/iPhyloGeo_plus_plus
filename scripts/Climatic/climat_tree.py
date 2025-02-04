import plotly.graph_objs as go

from Bio import Phylo
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

from PyQt6.QtWidgets import QDialog

from Climatic.climatic_graph_settings import ClimaticGraphSettings
from Climatic.climatic_preferences_dialog import ClimaticPreferencesDialog

from utils.download_file import download_file_local, download_file_temporary_PLT, download_file_temporary_PIO
from utils.error_dialog import show_error_dialog

try:
    ClimaticGraphSettings.load_from_file("./scripts/utils/ClimaticGraphSettings.yaml")
except FileNotFoundError:
    # Use default settings if the file is not found
    pass


class ClimaticTree():
    def __init__(self, main):
        self.climaticTrees = dict()
        self.main = main
        
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
                # reticulation_color = ClimaticGraphSettings.reticulation_color
                layout = ClimaticGraphSettings.layout
                proportional_edge_lengths = ClimaticGraphSettings.proportional_edge_lengths
                label_internal_vertices = ClimaticGraphSettings.label_internal_vertices
                use_leaf_names = ClimaticGraphSettings.use_leaf_names
                show_branch_length = ClimaticGraphSettings.show_branch_length
                view_type = ClimaticGraphSettings.view_type
                
                tree = self.check_render_properties_tree(tree, proportional_edge_lengths, label_internal_vertices)
                
                if view_type == "network":
                    self.render_network_view(
                        tree,
                        label_color,
                        edge_color,
                        layout,
                        use_leaf_names,
                        show_branch_length,
                    )
                else:
                    self.render_tree_view(
                        tree,
                        label_color,
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
            
    def check_render_properties_tree(self, tree, proportional_edge_lengths, label_internal_vertices):
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
                    
        return tree

    def render_network_view(
        self,
        tree,
        label_color,
        edge_color,
        layout,
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

            pixmap = download_file_temporary_PIO(self.tree_keys[self.current_index], fig)

            self.main.climaticTreesLabel.clear()
            self.main.climaticTreesLabel.setPixmap(pixmap)
            self.main.climaticTreesLabel.adjustSize()
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while rendering the network view: {e}")

    def render_tree_view(
        self,
        tree,
        label_color,
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

            pixmap = download_file_temporary_PLT(self.tree_keys[self.current_index], fig)

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

#NOSELF
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

#NOSELF
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
            download_file_local(current_key, self.main)
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