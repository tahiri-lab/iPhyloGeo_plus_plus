from utils.error_dialog import show_error_dialog
from utils.download_file import download_file_local
from utils.tree_graph import init_tree, generate_tree
from event_connector  import blocked_signals

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
        
        self.jsonData, self.tree_keys, formatted_tree_keys = init_tree()

        with blocked_signals(self.main.geneticTreescomboBox):
            self.main.geneticTreescomboBox.clear()
            self.main.geneticTreescomboBox.addItems(formatted_tree_keys)

        self.show_tree(0)

    def show_tree(self, index):
        """
        Display the phylogenetic tree at the specified index using Toytree.

        Args:
            index (int): The index of the tree to display.

        Returns:
            None
        """
        self.key = self.tree_keys[index]
        
        pixmap = generate_tree(self.key, self.jsonData)

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
            download_file_local(self.key, self.main)
        except FileNotFoundError as e:
            show_error_dialog(f"The tree image file was not found: {e}", "File Not Found")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while downloading the tree image: {e}")
