import matplotlib
matplotlib.use('Agg')  # Prevents matplotlib from opening a new window

import os
from arrow import get
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, ConnectionStyle
import matplotlib.colors as mcolors  
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from Bio import Phylo


def generate_tree_map(treeData, gps):
    # --------------------------------------
    # ------------  Phylogenetic MAP -------
    # --------------------------------------
    # Create a new map with PlateCarree projection
    fig = plt.figure(figsize=(26, 11))

    # Calculate positions for all nodes
    y_step = 1
    tree = treeData
    calc_node_positions(tree.root, 0, 0, y_step)

    # Create a figure for the subplot
    ax_tree = fig.add_subplot(121)

    # --------------------------------------
    # ------------  Color Groups -----------
    # --------------------------------------
    # Grouping by same coordinates and associating the same color

    gps_colors = gps.groupby(['LONG', 'LAT']).size().reset_index(name='count')

    # Use a colormap to generate distinct colors
    cmap = plt.get_cmap('tab20')
    num_colors = len(gps_colors)
    colors = [mcolors.to_hex(cmap(i / num_colors)) for i in range(num_colors)]

    gps_colors['group_color'] = colors

    def get_color(longitude, latitude):
        return gps_colors[(gps_colors['LONG'] == longitude) & (gps_colors['LAT'] == latitude)]['group_color'].values[0]

    gps['color'] = gps.apply(lambda row: get_color(row['LONG'], row['LAT']), axis=1)
    text_colors_id = dict(zip(gps['id'], gps['color']))
    text_colors_coords = dict(zip(zip(gps['LONG'], gps['LAT']), gps['color']))

    # --------------------------------------
    # ------------ Drawing Tree ------------
    # --------------------------------------

    # Plot the tree
    Phylo.draw(tree, do_show=False, axes=ax_tree, label_func=custom_label, label_colors=text_colors_id)
    ax_tree.set_title("", fontsize=18)

    # Set axes limits to verify the data range
    ax_tree.set_xlim(0, max(node.position[0] for node in tree.get_terminals()) + 1)
    ax_tree.set_ylim(0, max(node.position[1] for node in tree.get_terminals()) + 2)

    # Generate DataFrame with node coordinates
    rows = []
    for clade in tree.find_clades():
        if clade.is_terminal():
            label = clade.name
            color = text_colors_id.get(label)
            x, y = plot_adjusted_node(ax_tree, clade, y_step, color)
            rows.append([label, (x, y)])

    # Create DataFrame with node coordinates
    df = pd.DataFrame(rows, columns=["ID", "Coordinates"])

    # --------------------------------------
    # ------------  GRAPH MAP --------------
    # --------------------------------------

    extent = calculate_extent(gps)

    # Create subplot 2 with the map plot
    ax2 = fig.add_subplot(122, projection=ccrs.PlateCarree())
    ax2.set_extent(extent)

    # Plot points from GPS dataframe on the map
    for index, row in gps.iterrows():
        ax2.plot(
            row["LONG"], row["LAT"], 
            'o',  # Circle marker
            markersize=8,  # Same size as in ax.plot
            markerfacecolor=row['color'],  # Use the color from the 'color' column
            markeredgewidth=2,  # Same edge width
            markeredgecolor='black'  # Same edge color
        )

    # Add coastlines and country borders for context
    ax2.coastlines(resolution='10m')
    ax2.add_feature(cfeature.LAND)
    ax2.add_feature(cfeature.OCEAN)
    ax2.add_feature(cfeature.COASTLINE)
    ax2.add_feature(cfeature.BORDERS)

    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")
    ax2.set_title("Species Coordinates", fontsize=18)

    # --------------------------------------
    # ------------  Line mapping -----------
    # --------------------------------------

    # Group gps DataFrame by ID and create a dictionary of lists of coordinates
    gps_grouped = gps.groupby('id')[['LONG', 'LAT']].apply(
        lambda x: list(zip(x['LONG'], x['LAT']))).to_dict()

    for index, row in df.iterrows():
        # Get corresponding list of coordinates and color from gps DataFrame
        if row['ID'] in gps_grouped:
            species_coords_list = gps_grouped[row['ID']]

            num_coords = len(df)

            # Create connection patches for each coordinate in the list
            for species_coords in species_coords_list:
                longitude, latitude = species_coords  # Unpack coordinates and color
                rad_value = 0.15 - (0.3 * index / (num_coords - 1))
                con = ConnectionPatch(
                    xyA=row['Coordinates'], coordsA="data",
                    xyB=(longitude, latitude), coordsB="data",
                    axesA=ax_tree, axesB=ax2,
                    color=text_colors_coords.get((longitude, latitude)),
                    linewidth=3, linestyle="-", alpha=0.5,
                    zorder=2,
                    connectionstyle= ConnectionStyle("Arc3", rad=rad_value)
                )
                fig.add_artist(con)

    plt.tight_layout()
    output_file = 'results/figure_tree2map.png'
    plt.savefig(output_file, format='png')
    plt.close(fig)



def custom_label(clade):
    if clade.is_terminal():
        return clade.name
    else:
        return None

def calc_node_positions(tree, x_start, y_start, y_step):
    if tree.is_terminal():
        MINOFFSET = 10
        x_pos = x_start + tree.branch_length if tree.branch_length else x_start
        name_length = len(tree.name) if tree.name else 0
        if name_length <= MINOFFSET:
            x_pos += 0.01 * name_length
        else:
            x_pos += 0.009 * name_length
        y_pos = y_start
        y_start += y_step
    else:
        x_pos = x_start + tree.branch_length if tree.branch_length else x_start
        y_pos = y_start

        child_y_start = y_start
        for child in tree.clades:
            child_x_pos, child_y_pos, y_start = calc_node_positions(child, x_pos, y_start, y_step)

        y_pos = (y_start + child_y_start) / 2

    tree.position = (x_pos, y_pos)
    return x_pos, y_pos, y_start

def plot_adjusted_node(ax, node, y_offset, color):
    x, y = node.position
    y += y_offset

    ax.plot(x, y, 'o', markersize=8, markeredgewidth=2, markeredgecolor='black', markerfacecolor=color)
    return x, y

def calculate_extent(gps, margin_factor=0.2):
    """
    Calculate the extent of the map based on GPS data and a margin factor.

    Args:
        gps (DataFrame): DataFrame containing GPS data with 'LONG' and 'LAT' columns.
        margin_factor (float): Factor to extend the map boundaries.

    Returns:
        list: Extent of the map [min_lon, max_lon, max_lat, min_lat].
    """
    lon_margin = margin_factor * (gps['LONG'].max() - gps['LONG'].min())
    lat_margin = margin_factor * (gps['LAT'].max() - gps['LAT'].min())
    min_lon = gps['LONG'].min() - lon_margin
    max_lon = gps['LONG'].max() + lon_margin
    min_lat = gps['LAT'].min() - lat_margin
    max_lat = gps['LAT'].max() + lat_margin
    return [min_lon, max_lon, max_lat, min_lat]
