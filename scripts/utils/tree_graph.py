import os
import json
import toytree
import toyplot.png
import toyplot.color

import numpy as np
from numpy.typing import ArrayLike

from typing import cast
from PyQt6.QtGui import QPixmap


def init_tree():
     
    file_path = "results/geneticTrees.json"
    with open(file_path, "r") as file:
        newick_json = json.load(file)

    tree_keys = list(newick_json.keys())

    # Format the tree keys to replace underscore with ' nt '
    formatted_tree_keys = [format_tree_name(key) for key in tree_keys]
    
    return newick_json, tree_keys, formatted_tree_keys
    
def format_tree_name(tree_name):
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


def generate_tree(key, jsonData, isDarkMode):
    """
    Display the phylogenetic tree at the specified index using Toytree.

    Args:
       key (string): The json key
       jsonData (json): The json Data

    Returns:
        Pixmap
    """

    tree, canvas, _ = create_tree_variables(jsonData, key)
    
    color = "white" if isDarkMode else "black"
    
    palette = toyplot.color.Palette([color])
    ax = canvas.cartesian(bounds=(50, 870, 50, 400), padding=15, palette=palette)
    tree.draw(axes=ax, edge_colors=color, node_colors=color, tip_labels_colors=color,)
    
    
    ax.x.label.style = {"fill": color}
    ax.y.label.style = {"fill": color}
    ax.x.ticks.labels.style = {"fill": color}
    ax.y.ticks.labels.style = {"fill": color}
    ax.x.spine.style = {"stroke": color}
    ax.y.spine.style = {"stroke": color}

    pixmap = download_and_render(canvas, key)
    
    return pixmap

def generate_tree_with_bar(key, jsonData, selected_column, climatic_data, isDarkMode ):
    """
    Display the phylogenetic tree at the specified index using Toytree with bars on his right.

    Args:
       key (string): The json key
       jsonData (json): The json Data
       selected_column (string): Column name for the bar
       climatic_data (DataFrame): Data related to Climate 

    Returns:
        Pixmap
    """
    tree, canvas, tip_labels = create_tree_variables(jsonData, key)

    color = "white" if isDarkMode else "black"
    
    ax = canvas.cartesian(bounds=(50, 300, 50, 400), ymin=0, ymax=tree.ntips, padding=15)
    tree.draw(axes=ax, tip_labels=[label.replace("_","") for label in tip_labels], tip_labels_colors=color, edge_colors=color, node_colors=color)
    ax.show = False
    
    climatic_bar(canvas, tree, tip_labels, selected_column, climatic_data, color)

    pixmap = download_and_render(canvas, key)
    
    return pixmap

def create_tree_variables(jsonData, key):
    
    newick_str = jsonData[key]

    tree = toytree.tree(newick_str)   

    tip_labels = tree.get_tip_labels()

    canvas = toyplot.Canvas(width=921, height=450)
    
    return tree, canvas, tip_labels


def download_and_render(canvas, key):
    tree_img_path = os.path.join("results", f"{key}.png")
    os.makedirs(os.path.dirname(tree_img_path), exist_ok=True)
    toyplot.png.render(canvas, tree_img_path)

    pixmap = QPixmap(tree_img_path)
    return pixmap
  

def climatic_bar(canvas, tree, tip_labels, selected_column, climatic_data, color):

    # Match the original labels from the CSV
    ordered_climatic_data = climatic_data.set_index("id").reindex(tip_labels).reset_index()
    bar_values = ordered_climatic_data[selected_column].values

    # Ensure bar_values have no negative values
    bar_values = np.clip(cast(ArrayLike, bar_values), a_min=0, a_max=None)

    # Add bar plot to canvas
    ax1 = canvas.cartesian(bounds=(325, 900, 50, 400), ymin=0, ymax=tree.ntips, padding=15)
    ax1.bars(
        np.arange(tree.ntips),
        bar_values,
        along="y",
    )
    
    ax1.show = True
    ax1.y.show = False
    ax1.x.ticks.show = True
    ax1.x.label.text = selected_column
    ax1.x.label.style = {"fill": color}
    ax1.x.ticks.labels.style = {"fill": color}
    ax1.x.spine.style = {"stroke": color}