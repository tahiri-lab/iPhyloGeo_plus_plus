import matplotlib.pyplot as plt
import networkx as nx
import plotly.graph_objs as go
import seaborn as sns
from Bio import Phylo
from Climatic.climatic_graph_settings import ClimaticGraphSettings
from utils.error_dialog import show_error_dialog


def get_network_graph(tree):
    try:
        sns.set_palette("husl")

        graph = Phylo.to_networkx(tree)
        pos = get_layout(graph)

        edge_trace_result = create_edge_trace(tree, pos)
        if edge_trace_result is None:
            raise TypeError("The edge_trace is None")
        edge_trace, edge_annotations = edge_trace_result

        node_trace = create_node_trace(graph, pos)

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

        return fig

    except Exception as e:
        show_error_dialog(f"An unexpected error occurred while rendering the network view: {e}")


def get_tree_graph(tree):
    try:
        label_color = ClimaticGraphSettings.label_color
        use_leaf_names = ClimaticGraphSettings.use_leaf_names
        show_branch_length = ClimaticGraphSettings.show_branch_length
        edge_color = ClimaticGraphSettings.edge_color

        fig = plt.figure(figsize=(9.11, 4.41))  # Limit size to 911x441 pixels
        ax = fig.add_subplot(1, 1, 1)

        for clade in tree.find_clades():
            clade.color = edge_color

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

        return fig
    except Exception as e:
        show_error_dialog(f"An unexpected error occurred while rendering the tree view: {e}")


def create_node_trace(graph, pos):
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
        use_leaf_names = ClimaticGraphSettings.use_leaf_names
        label_color = ClimaticGraphSettings.label_color
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
            textfont=dict(color=label_color),
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


def get_layout(graph):
    """
    Get the layout for the phylogenetic tree network visualization.

    Args:
        graph (networkx.Graph): The networkx graph of the phylogenetic tree.
        layout (str): The layout type (e.g., 'horizontal', 'vertical', 'radial', 'axial').

    Returns:
        dict: A dictionary of positions keyed by node.
    """
    try:
        layout = ClimaticGraphSettings.layout

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


def create_edge_trace(tree, pos):
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
        edge_color = ClimaticGraphSettings.edge_color
        show_branch_length = ClimaticGraphSettings.show_branch_length

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
