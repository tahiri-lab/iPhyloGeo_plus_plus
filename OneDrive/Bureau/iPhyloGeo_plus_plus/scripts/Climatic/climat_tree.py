from Climatic.climat_tree_graph import get_network_graph, get_tree_graph
from Climatic.climatic_graph_settings import ClimaticGraphSettings
from Climatic.climatic_preferences_dialog import ClimaticPreferencesDialog
from event_connector import blocked_signals
from PyQt6.QtWidgets import QDialog
from utils.download_file import download_file_local, download_file_temporary_PIO, download_file_temporary_PLT
from utils.error_dialog import show_error_dialog

try:
    ClimaticGraphSettings.load_from_file("ClimaticGraphSettings.yaml")
except FileNotFoundError:
    # Use default settings if the file is not found
    pass


class ClimaticTree:
    def __init__(self, main):
        self.climaticTrees = dict()
        self.main = main

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
            self.tree_keys = list(self.climaticTrees.keys())
            self.total_trees = len(self.tree_keys)

            with blocked_signals(self.main.climaticTreescomboBox):
                self.main.climaticTreescomboBox.clear()
                self.main.climaticTreescomboBox.addItems(self.tree_keys)

            self.current_key = self.tree_keys[0]
            self.show_climatic_tree(0)

            self.main.tabWidget2.setCurrentIndex(3)

        except KeyError as e:
            show_error_dialog(
                f"An error occurred while accessing the climatic trees: {e}",
                "Key Error",
            )
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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
                self.current_key = self.tree_keys[index]
                tree = self.climaticTrees[self.current_key]
                tree = check_render_properties_tree(tree)

                view_type = ClimaticGraphSettings.view_type
                if view_type == "network":
                    self.render_network_view(tree)
                else:
                    self.render_tree_view(tree)
        except KeyError as e:
            show_error_dialog(
                f"An error occurred while accessing the climatic tree data: {e}",
                "Key Error",
            )
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def render_network_view(self, tree):
        """
        Render the network view of the phylogenetic tree.

        This method applies user preferences to the tree and renders it using Plotly for interactive visualization.

        Args:
            tree (Phylo.BaseTree.Tree): The phylogenetic tree to render.
        Returns:
            None
        """
        try:
            fig = get_network_graph(tree)

            pixmap = download_file_temporary_PIO(self.current_key, fig)

            self.main.climaticTreesLabel.clear()
            self.main.climaticTreesLabel.setPixmap(pixmap)
            self.main.climaticTreesLabel.adjustSize()
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while rendering the network view: {e}")

    def render_tree_view(self, tree):
        """
        Render the tree view of the phylogenetic tree.

        This method applies user preferences to the tree and renders it using Matplotlib.

        Args:
            tree (Phylo.BaseTree.Tree): The phylogenetic tree to render.

        Returns:
            None
        """
        try:
            fig = get_tree_graph(tree)

            pixmap = download_file_temporary_PLT(self.current_key, fig)

            self.main.climaticTreesLabel.clear()
            self.main.climaticTreesLabel.setPixmap(pixmap)
            self.main.climaticTreesLabel.adjustSize()
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while rendering the tree view: {e}")

    def download_climatic_tree_graph(self):
        """
        Download the current displayed climatic tree graph as a PNG file.

        This method prompts the user to select a location to save the current graph,
        and saves the graph to the specified location.

        Returns:
            None
        """
        try:
            download_file_local(self.current_key)
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
            self.show_climatic_tree(self.main.climaticTreescomboBox.currentIndex())


def check_render_properties_tree(tree):
    label_internal_vertices = ClimaticGraphSettings.label_internal_vertices
    proportional_edge_lengths = ClimaticGraphSettings.proportional_edge_lengths
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
