import json
from ete3 import Tree

def process_newick_trees(file_path):
    """
    Reads Newick trees from a JSON file, converts them into ete3 Tree objects,
    and provides options for visualization and analysis.

    Args:
        file_path (str): Path to the JSON file containing the Newick trees.

    Returns:
        list: A list of ete3 Tree objects.
    """
    with open(file_path, "r") as f:
        tree_data = json.load(f)

    ete_trees = []
    for key, newick_string in tree_data.items():
        t = Tree(newick_string, format=1)  # format=1 for Newick
        t.name = key  # Store the JSON key as the tree name for reference
        ete_trees.append(t)

    # Example usage:
    for t in ete_trees:
        print(f"Tree '{t.name}':")
        print(t)  # Print the Newick representation
        # Uncomment for visualization:
        t.show()
        # Further analysis or modification using ete3 functions can be done here

    return ete_trees

# Example usage:
file_path = "./results/geneticTrees.json"  # Replace with your file path
trees = process_newick_trees(file_path)
