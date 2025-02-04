import json
import toyplot.png
import toytree
import os
from PyQt6.QtGui import QPixmap

from utils.error_dialog import show_error_dialog
from utils.download_file import download_file_local

class GeneticTree:
    def __init__(self, main):
        self.main = main
        
    def display_newick_trees(self):
        """
        Display Newick format trees in the application using Toytree.

        This method loads Newick format trees from a JSON file, formats the tree names,
        and updates the UI to display the trees using Toytree.

        Actions:
            - Loads Newick format trees from 'results/geneticTrees.json'.
            - Formats the tree names by replacing underscores with ' nt '.
            - Updates the tree combo box with formatted tree names.
            - Displays the first tree in the list.
        """

        self.main.tabWidget.setCurrentIndex(4)
        file_path = "results/geneticTrees.json"
        with open(file_path, "r") as file:
            self.newick_json = json.load(file)

        self.tree_keys = list(self.newick_json.keys())
        self.total_trees = len(self.tree_keys)
        self.current_index = 0
        self.main.geneticTreescomboBox.clear()

        # Format the tree keys to replace underscore with ' nt '
        formatted_tree_keys = [self.format_tree_name(key) for key in self.tree_keys]
        self.main.geneticTreescomboBox.addItems(formatted_tree_keys)

        self.show_tree(self.current_index)

    def format_tree_name(self, tree_name):
        """
        Format the tree name by replacing underscores with ' nt '.

        Args:
            tree_name (str): The original tree name.

        Returns:
            str: The formatted tree name.
        """
        parts = tree_name.split("_")
        if len(parts) == 2:
            return f"{parts[0]} nt {parts[1]} nt"
        return tree_name

    def show_tree(self, index):
        """
        Display the phylogenetic tree at the specified index using Toytree.

        Args:
            index (int): The index of the tree to display.

        Returns:
            None
        """
        if index is None or index < 0 or index >= self.total_trees:
            return

        self.current_index = index  # Keep track of the current index
        key = self.tree_keys[index]  # This is the key with underscores
        newick_str = self.newick_json[key]

        # Read the tree using Toytree
        tree = toytree.tree(newick_str)

        # Replace underscores with spaces in tip labels
        for node in tree.treenode.traverse():
            if node.is_leaf():
                node.name = node.name.replace("_", " ")  # Replace underscores with spaces

        # Customize the tip labels and their style
        tip_labels = tree.get_tip_labels()

        # Draw the tree with customized style
        canvas, _, _ = tree.draw(
            width=921,
            height=450,
            tip_labels=tip_labels,  # These labels now have spaces
            tip_labels_style={"font-size": "15px"},
            fixed_order=tip_labels,
            tip_labels_colors="black",
            edge_type="c",
        )

        # Adjust the canvas size to ensure it fits within the specified dimensions
        canvas = toyplot.Canvas(width=921, height=450)
        ax = canvas.cartesian(bounds=(50, 870, 50, 400), padding=15)
        tree.draw(axes=ax)

        # Save the canvas to a permanent file in the .results/ directory
        self.tree_img_path = os.path.join("results", f"{key}.png")
        os.makedirs(os.path.dirname(self.tree_img_path), exist_ok=True)
        toyplot.png.render(canvas, self.tree_img_path)

        # Create a QPixmap from the saved image file
        pixmap = QPixmap(self.tree_img_path)

        # Clear the QLabel before setting the new QPixmap
        self.main.GeneticTreeLabel.clear()
        self.main.GeneticTreeLabel.setPixmap(pixmap)
        self.main.GeneticTreeLabel.adjustSize()

    def download_genetic_tree_graph(self):
        """
        Download the current displayed tree graph as a PNG file.

        This method prompts the user to select a location to save the current tree graph,
        and saves the graph to the specified location.

        Returns:
            None
        """
        try:
            current_key = self.tree_keys[self.current_index]
            download_file_local(current_key, self.main)
        except FileNotFoundError as e:
            show_error_dialog(f"The tree image file was not found: {e}", "File Not Found")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while downloading the tree image: {e}")
