import io
import json
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from decimal import Decimal

import folium
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import qtmodern.styles
import qtmodern.windows
import seaborn as sns
import toyplot.png
import toytree
from aphylogeo import utils
from aphylogeo.alignement import AlignSequences
from aphylogeo.genetic_trees import GeneticTrees
from aphylogeo.params import Params
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.ticker import MaxNLocator
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import QObject, Qt, QThread, pyqtSignal
from PyQt6.QtGui import QColor, QIcon, QMovie, QPixmap
from PyQt6.QtWidgets import QApplication, QDialog, QFileDialog, QGraphicsDropShadowEffect, QTableWidget, QTableWidgetItem, QVBoxLayout
from PyQt6.QtWebEngineWidgets import QWebEngineView
from Qt import main, loading
from ui_helpers import create_shadow_effect, get_button_style, style_buttons
from utils import resources_rc  # noqa: F401  # Import the compiled resource module for resolving image resource path
from utils.genetic_params_dialog import ParamDialog
from utils.help import HelpDialog
from utils.my_dumper import update_yaml_param
from utils.PreferencesDialog import PreferencesDialog  # Import PreferencesDialog
from utils.result_settings_dialog import ResultSettingsDialog
from utils.settings import Params2

from event_connector import QtEvents, connect_event, connect_decorated_methods

try:
    Params.load_from_file("./scripts/utils/params.yaml")
except FileNotFoundError:
    Params.validate_and_set_params(Params2.PARAMETER_KEYS)


class Worker(QObject):
    progress = pyqtSignal(int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, filepath):
        super().__init__()
        self.filepath = filepath
        self.running = True

    def run(self):
        try:
            # Step 1: Load sequences
            self.progress.emit(0)
            if not self.running:
                return
            sequenceFile = utils.loadSequenceFile(self.filepath)  # noqa: N806

            # Step 2: Align sequences
            self.progress.emit(1)
            if not self.running:
                return
            align_sequence = AlignSequences(sequenceFile)
            alignments = align_sequence.align()

            # Step 3: Generate genetic trees
            self.progress.emit(2)
            if not self.running:
                return
            geneticTrees = utils.geneticPipeline(alignments.msa)  # noqa: N806

            trees = GeneticTrees(trees_dict=geneticTrees, format="newick")

            # Step 4: Preparing results
            msa = alignments.to_dict().get("msa")

            # Step 5: Save results
            alignments.save_to_json(f"./scripts/results/aligned_{Params.reference_gene_file}.json")
            trees.save_trees_to_json("./scripts/results/geneticTrees.json")

            # Emit finished signal with the genetic trees dictionary
            result = {"msa": msa, "geneticTrees": geneticTrees}
            self.progress.emit(3)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(Exception(f"{e} consider to change the Tree Type in the alignment settings")))

    def stop(self):
        self.running = False


window_size = 50
starting_position = 1


class UiMainWindow(main.Ui_MainWindow, QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setup_ui()
        connect_decorated_methods(self)

    @connect_event("helpButton", QtEvents.clicked)
    def open_help_window(self):
        """
        Initialize and display the 'How to Use' window.

        This method creates a new QMainWindow instance, sets up its UI using the UiHowToUse class, and displays the window.
        """

        dialog = HelpDialog()
        dialog.exec()

    @connect_event("settingsButtonPage3", QtEvents.clicked)
    def open_result_settings_window(self):
        """
        Initialize and display the parameters window.
        This method creates a new QMainWindow instance, sets up its UI using the UiDialog class, and displays the window.
        """

        dialog = ResultSettingsDialog()
        dialog.exec()

    def show_error_dialog(self, message, title="error"):
        """
        Display a professional error dialog with the given title and message.

        Args:
            title (str): The title of the error dialog.
            message (str): The error message to display.
        """
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setWindowTitle(title)
        msgBox.setText(message)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.Ok)
        msgBox.exec()

    def setup_ui(self):
        """
        Setup the UI components and initialize the main window.

        This method connects various UI buttons to their corresponding event handlers,
        sets up styles and effects for UI elements, and initializes the state of the application.
        """
        try:
            self.preferences = {
                "label_color": "black",
                "edge_color": "blue",
                "reticulation_color": "red",
                "layout": "horizontal",
                "proportional_edge_lengths": False,
                "label_internal_vertices": False,
                "use_leaf_names": True,
                "show_branch_length": False,
                "view_type": "network",
            }
            self.tree_keys = []
            self.total_trees = 0
            self.current_index = 0
            self.current_index1 = 0
            self.setObjectName("MainWindow")
            self.window_size_spinbox_2.setRange(1, 1000)
            self.starting_position_spinbox_2.setRange(1, 1000)
            self.isDarkMode = False  # Keep track of the state
            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
            self.stackedWidget.setCurrentIndex(0)

            buttons = [
                self.geneticDataButton,
                self.climaticDataButton,
                self.helpButton,
                self.homeButton,
                self.resultsButton,
            ]
            buttons_Vertical = [
                self.fileBrowserButtonPage1,
                self.sequenceAlignmentButtonPage1,
                self.clearButtonPage1,
                self.statisticsButtonPage1,
                self.geneticTreeButtonPage1,
                self.fileBrowserButtonPage2,
                self.clearButtonPage2,
                self.climaticTreeButtonPage2,
                self.statisticsButtonPage2,
                self.settingsButtonPage3,
                self.settingsButtonPage4,
                self.submitButtonPage3,
                self.statisticsButtonPage3,
                self.submitButtonPage4,
                self.statisticsButtonPage4,
                self.StartSequenceAlignmentButton,
                self.clearButtonPage3,
                self.clearButtonPage4,
            ]

            # Define cursor and stylesheet for all buttons
            for button in buttons:
                button.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
                button.setStyleSheet(
                    """
                    QPushButton {
                        padding: 10px 20px;
                        font-weight: bold;
                        background-color: #DEDDDA;
                        border-radius: 20px;
                        transition: background-color 0.3s ease;
                    }
                    QPushButton:hover {
                        background-color: #B7B7B6;
                    }
                    QPushButton:pressed {
                        background-color: #DEDDDA;
                    }
                """
                )
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 140))
                shadow_effect.setOffset(3, 3)
                button.setGraphicsEffect(shadow_effect)

            for buttonV in buttons_Vertical:
                buttonV.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
                buttonV.setStyleSheet(
                    """
                    QPushButton {
                        border-radius: 14px;
                        background-color: #EEEEEE;
                        padding: 10px 20px;
                        font-weight: bold;
                        transition: background-color 0.3s ease;
                    }
                    QPushButton:hover {
                        background-color: #D7D7D7;
                    }
                    QPushButton:pressed {
                        background-color: #EEEEEE;
                    }
                """
                )
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 110))
                shadow_effect.setOffset(3, 3)
                buttonV.setGraphicsEffect(shadow_effect)

            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
            self.darkModeButton.setStyleSheet(
                """
                QPushButton {
                    padding: 10px 20px;
                    font-weight: bold;
                    background-color: #DEDDDA;
                    border-radius: 14px;
                    transition: background-color 0.3s ease;
                }
                QPushButton:hover {
                    background-color: #B7B7B6;
                }
                QPushButton:pressed {
                    background-color: #DEDDDA;
                }
            """
            )
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)
            shadow_effect.setColor(QColor(0, 0, 0, 140))
            shadow_effect.setOffset(3, 3)
            self.darkModeButton.setGraphicsEffect(shadow_effect)

            QtCore.QMetaObject.connectSlotsByName(self)

        except AttributeError as e:
            self.show_error_dialog(f"An error occurred while setting up the UI: {e}", "Attribute Error")
        except Exception as e:
            self.show_error_dialog(
                f"An unexpected error occurred: {e}",
                "Unexpected Error",
            )

    @connect_event("clearButtonPage3", QtEvents.clicked)
    def clear_results(self):
        self.textEditResults.clear()

    @connect_event("geneticSettingsButton", QtEvents.clicked)
    def open_genetic_settings_window(self):
        dialog = ParamDialog()
        if dialog.exec() == QDialog.accepted:
            self.geneticParam = dialog.params
            for property_name, new_value in self.geneticParam.items():
                update_yaml_param(Params, "scripts/utils/params.yaml", property_name, new_value)

    # Update_plot_start

    def read_msa(self, msa_data):
        """
        Reads multiple sequence alignment (MSA) data and organizes it into a dictionary.

        Args:
            msa_data (dict): A dictionary containing MSA data where keys are identifiers and values are sequences.

        Returns:
            dict: A dictionary with sequence identifiers as keys and concatenated sequences as values.
        """
        try:
            genetic_data = {}
            for key, value in msa_data.items():
                lines = value.strip().split("\n")
                current_id = None
                for line in lines:
                    if line.startswith(">"):
                        current_id = line[1:].strip()
                        if current_id not in genetic_data:
                            genetic_data[current_id] = []
                    else:
                        genetic_data[current_id].append(line.strip())

            genetic_data = {sequence_id: "".join(sequences) for sequence_id, sequences in genetic_data.items()}
            return genetic_data

        except KeyError as e:
            self.show_error_dialog(f"Key Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def standardize_sequence_lengths(self, genetic_data):
        """
        Standardizes the lengths of genetic sequences by padding shorter sequences with '-'.

        Args:
            genetic_data (dict): A dictionary with sequence identifiers as keys and sequences as values.

        Returns:
            dict: A dictionary with sequence identifiers as keys and standardized-length sequences as values.
        """
        try:
            max_length = max(len(seq) for seq in genetic_data.values())
            standardized_data = {key: seq.ljust(max_length, "-") for key, seq in genetic_data.items()}
            return standardized_data

        except ValueError as e:
            self.show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def plot_alignment_chart(self, genetic_data, starting_position, window_size, output_path):
        """
        Plots an alignment chart with conservation and sequence alignment.

        Args:
            genetic_data (dict): A dictionary with sequence identifiers as keys and sequences as values.
            starting_position (int): The starting position for the alignment window.
            window_size (int): The size of the alignment window.
            output_path (str): The path to save the output plot.

        Returns:
            None
        """
        try:
            # Replace underscores with spaces in the keys of genetic_data
            genetic_data = {key.replace("_", " "): value for key, value in genetic_data.items()}

            end_position = starting_position + window_size
            truncated_data = {key: value[starting_position:end_position] for key, value in genetic_data.items()}
            alignment = MultipleSeqAlignment([SeqRecord(Seq(seq), id=key) for key, seq in truncated_data.items()])

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 4), gridspec_kw={"height_ratios": [1, 8]})
            ax1.set_axis_off()
            ax2.set_axis_off()

            # Calculate conservation
            conservation, _ = self.calculate_conservation_and_gaps(alignment)

            # Plot conservation
            bar_width = 1.0
            ax1.bar(
                range(len(conservation)),
                conservation,
                color="#4CAF50",
                width=bar_width,
                align="edge",
            )
            ax1.set_xlim(0, len(conservation))
            ax1.set_ylim(0, 1)
            ax1.set_title("CONSERVATION", fontsize=12, pad=20)

            # Plot alignment
            seqs = [str(record.seq) for record in alignment]
            ids = [record.id for record in alignment]
            colors = {
                "A": "#4CAF50",
                "T": "#F44336",
                "G": "#2196F3",
                "C": "#FFEB3B",
                "-": "grey",
                "N": "black",  # Handling 'N' or any other characters if present
            }

            font_size = 10
            rect_height = 0.8
            rect_width = 1.0

            for i, seq in enumerate(seqs):
                ax2.text(
                    -1,
                    len(seqs) - i - 1 + rect_height / 2,
                    ids[i],
                    ha="right",
                    va="center",
                    fontsize=font_size,
                )
                for j, nucleotide in enumerate(seq):
                    color = colors.get(nucleotide, "black")  # Default to black if not found
                    rect = mpatches.Rectangle((j, len(seqs) - i - 1), rect_width, rect_height, color=color)
                    ax2.add_patch(rect)
                    ax2.text(
                        j + rect_width / 2,
                        len(seqs) - i - 1 + rect_height / 2,
                        nucleotide,
                        ha="center",
                        va="center",
                        fontsize=font_size,
                    )

            consensus_seq = self.calculate_consensus(alignment)
            for j, nucleotide in enumerate(consensus_seq):
                color = colors.get(nucleotide, "black")  # Default to black if not found
                rect = mpatches.Rectangle((j, -1), rect_width, rect_height, color=color)
                ax2.add_patch(rect)
                ax2.text(
                    j + rect_width / 2,
                    -1 + rect_height / 2,
                    nucleotide,
                    ha="center",
                    va="center",
                    fontsize=font_size,
                    fontweight="bold",
                )

            # Add Consensus title
            ax2.text(
                -1,
                -1 + rect_height / 2,
                "Consensus",
                ha="right",
                va="center",
                fontsize=font_size,
                fontweight="bold",
            )

            ax2.set_xlim(0, len(consensus_seq))
            ax2.set_ylim(-2, len(seqs))
            ax2.set_yticks(range(len(ids)))
            ax2.set_yticklabels(ids, fontsize=font_size)

            # Ensure number of ticks matches the length of the sequences
            ax2.set_xticks(range(len(consensus_seq)))
            ax2.set_xticklabels(
                range(starting_position, starting_position + len(consensus_seq)),
                fontsize=font_size,
            )

            plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.05)
            plt.savefig(output_path)
            plt.close()

        except KeyError as e:
            self.show_error_dialog(f"Key Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def calculate_conservation_and_gaps(self, alignment):
        """
        Calculate conservation and gap frequencies for a given multiple sequence alignment.

        Args:
            alignment (MultipleSeqAlignment): A BioPython MultipleSeqAlignment object containing the aligned sequences.

        Returns:
            tuple: A tuple containing two lists:
                - conservation (list): A list of conservation scores for each column in the alignment.
                - gaps (list): A list of gap frequencies for each column in the alignment.
        """
        try:
            length = alignment.get_alignment_length()
            conservation = []
            gaps = []

            for i in range(length):
                column = alignment[:, i]
                counter = Counter(column)
                most_common = counter.most_common(1)[0][1]
                conservation.append(most_common / len(column))
                gaps.append(counter["-"] / len(column))

            return conservation, gaps

        except KeyError as e:
            self.show_error_dialog(f"Key Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def calculate_consensus(self, alignment):
        """
        Calculate the consensus sequence for a given multiple sequence alignment.

        Args:
            alignment (MultipleSeqAlignment): A BioPython MultipleSeqAlignment object containing the aligned sequences.

        Returns:
            str: A string representing the consensus sequence.
        """
        try:
            length = alignment.get_alignment_length()
            consensus = []

            for i in range(length):
                column = alignment[:, i]
                counter = {}
                for base in column:
                    if base in counter:
                        counter[base] += 1
                    else:
                        counter[base] = 1
                most_common_base = max(counter, key=counter.get)
                consensus.append(most_common_base)

            return "".join(consensus)

        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event(["starting_position_spinbox_2", "window_size_spinbox_2"], QtEvents.valueChanged)
    def update_plot(self):
        """
        Update the plot based on the current starting position and window size.

        This method reads the MSA data, standardizes the sequence lengths, and plots the alignment chart.
        The resulting plot is displayed in the specified label widget and the tab widget is updated.

        Returns:
            None
        """
        try:
            starting_position = self.starting_position_spinbox_2.value()
            window_size = self.window_size_spinbox_2.value()
            output_path = "scripts/results/sequence_alignment_plot.png"

            genetic_data = self.read_msa(self.msa)
            standardized_data = self.standardize_sequence_lengths(genetic_data)
            self.plot_alignment_chart(standardized_data, starting_position, window_size, output_path)

            pixmap = QPixmap(output_path)
            self.seqAlignLabel.setPixmap(pixmap)
            self.tabWidget.setCurrentIndex(2)

        except AttributeError as e:
            self.show_error_dialog(f"Attribute Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    # Update_plot_end

    ################################

    # bug que plus que tu click sur le boutton, plus il y a de page qui ouvre
    @connect_event("statisticsButtonPage1", QtEvents.clicked)
    def initialize_species_list(self):
        # Load species names into the combo box
        self.referenceComboBox.clear()
        unique_species = set()
        for key, value in self.msa.items():
            parts = value.strip().split("\n")
            for i in range(0, len(parts), 2):
                header = parts[i].strip(">").replace("_", " ")  # Replace underscores with spaces
                unique_species.add(header)

        # Add unique species to the combo box
        for species in unique_species:
            self.referenceComboBox.addItem(species)

        # Optionally set the first species as the default selected item
        if self.referenceComboBox.count() > 0:
            self.referenceComboBox.setCurrentIndex(0)

        # Connect the value change signals to update the plot
        self.similarityWindowSizeSpinBox.valueChanged.connect(self.update_similarity_plot)
        self.startingPositionSimilaritySpinBox.valueChanged.connect(self.update_similarity_plot)
        self.referenceComboBox.currentIndexChanged.connect(self.update_similarity_plot)

        # Generate the initial plot
        self.update_similarity_plot()

    def update_similarity_plot(self):
        try:
            window_size = self.similarityWindowSizeSpinBox.value()
            start_pos = self.startingPositionSimilaritySpinBox.value()
            reference_species = self.referenceComboBox.currentText().replace(" ", "_")  # Convert back to original format

            sequences = defaultdict(str)
            for key, value in self.msa.items():
                parts = value.strip().split("\n")
                for i in range(0, len(parts), 2):
                    header = parts[i].strip(">")
                    sequence = parts[i + 1]
                    sequences[header.replace("_", " ")] += sequence  # Replace underscores with spaces

            if reference_species.replace("_", " ") not in sequences:
                print(f"Reference species {reference_species} not found in MSA data.")
                return

            max_len = max(len(seq) for seq in sequences.values())
            padded_records = []
            for header, sequence in sequences.items():
                padded_seq = sequence.ljust(max_len, "-")
                padded_records.append(SeqRecord(Seq(padded_seq), id=header))

            alignment = MultipleSeqAlignment(padded_records)
            reference_index = [record.id for record in alignment].index(reference_species.replace("_", " "))
            reference_sequence = str(alignment[reference_index].seq)
            similarities = []

            for record in alignment:
                similarity = [1 if ref == res else 0 for ref, res in zip(reference_sequence[start_pos:], str(record.seq)[start_pos:])]
                similarities.append(similarity)

            similarities = np.array(similarities)

            def sliding_window_avg(arr, window_size, step_size):
                return [np.mean(arr[i : i + window_size]) for i in range(0, len(arr) - window_size + 1, step_size)]

            step_size = 10  # Adjusted step size for better plotting

            windowed_similarities = []
            for sim in similarities:
                windowed_similarities.append(sliding_window_avg(sim, window_size, step_size))

            windowed_similarities = np.array(windowed_similarities)

            sns.set(style="whitegrid")
            fig, ax = plt.subplots(figsize=(12, 8))
            x = np.arange(start_pos, len(reference_sequence) - window_size + 1, step_size)

            for idx, record in enumerate(alignment):
                ax.plot(x, windowed_similarities[idx], label=record.id, linewidth=2.0)

            ax.set_xlabel("Position", fontsize=14)
            ax.set_ylabel("Similarity", fontsize=14)
            ax.set_title("Sequence Similarity Plot", fontsize=16, weight="bold")
            ax.legend(
                title="Species",
                fontsize=12,
                title_fontsize="13",
                loc="upper right",
                bbox_to_anchor=(1.2, 1),
            )
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.grid(True, linestyle="--", alpha=0.6)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_linewidth(1.2)
            ax.spines["bottom"].set_linewidth(1.2)

            self.plot_path = "./scripts/results/similarity_plot.png"
            fig.savefig(self.plot_path, bbox_inches="tight", dpi=300)

            pixmap = QPixmap(self.plot_path)
            self.textEditGenStats_2.setPixmap(pixmap)
            self.textEditGenStats_2.setFixedSize(900, 400)
            self.textEditGenStats_2.setScaledContents(True)
            self.tabWidget.setCurrentIndex(3)
        except Exception as e:
            print(f"Error updating similarity plot: {e}")

    @connect_event("downloadSimilarityButton", QtEvents.clicked)
    def download_similarity_plot_chart(self):
        file_url = "scripts/results/similarity_plot.png"  # The file path
        save_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", "PNG Files (*.png);;All Files (*)")
        if not save_path:
            return  # User cancelled th
        shutil.copy(file_url, save_path)  # e save dialog

    @connect_event("statisticsButtonPage2", QtEvents.clicked)
    def load_climate_statistics(self):
        """
        Load climate data from a CSV file, update the UI elements with the column names, and switch to the appropriate tab.

        This method reads the climate data from a CSV file (excluding the first column), updates combo boxes
        for X and Y axis settings with the column names, and switches to the climate data tab.

        Returns:
            None
        """
        # Load the data
        self.data = pd.read_csv(Params.file_name)
        self.columns = self.data.columns.tolist()

        # Remove the first column
        if self.columns:
            self.columns.pop(0)

        self.ClimaticChartSettingsAxisX.clear()
        self.ClimaticChartSettingsAxisX.addItems(self.columns)
        self.ClimaticChartSettingsAxisY.clear()
        self.ClimaticChartSettingsAxisY.addItems(self.columns)
        self.tabWidget2.setCurrentIndex(2)

    @connect_event(["ClimaticChartSettingsAxisX", "ClimaticChartSettingsAxisY", "PlotTypesCombobox"], QtEvents.currentIndexChanged)
    def generate_climate_graph(self):
        """
        Generate and display a graph based on the selected X and Y axis data and the chosen plot type.

        This method reads the selected data columns and plot type from the UI, generates the corresponding graph,
        and displays it in the specified QLabel widget.

        Returns:
            None
        """
        x_data = self.ClimaticChartSettingsAxisX.currentText()
        y_data = self.ClimaticChartSettingsAxisY.currentText()
        plot_type = self.PlotTypesCombobox.currentText()

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
        self.ClimaticChart_2.setPixmap(pixmap)
        self.tabWidget2.setCurrentIndex(2)

    @connect_event("climatePlotDownloadButton", QtEvents.clicked)
    def download_climate_plot(self):
        """
        Download the generated plot.

        This method allows the user to download the currently displayed plot.

        Returns:
            None
        """
        plot_type = self.PlotTypesCombobox.currentText()
        if plot_type == "":
            self.show_error_dialog("Please generate a plot first.")
            return

        plot_path = os.path.join("results", f"{plot_type.lower().replace(' ', '_')}.png")
        if not os.path.exists(plot_path):
            self.show_error_dialog("No plot found to download.")
            return

        # Prompt the user to select a location to save the plot
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Plot As",
            os.path.basename(plot_path),
            "PNG Files (*.png);;All Files (*)",
            options=options,
        )
        if file_path:
            if not file_path.lower().endswith(".png"):
                file_path += ".png"
            shutil.copy(plot_path, file_path)

    @connect_event("fileBrowserButtonPage1", QtEvents.clicked)
    def select_fasta_file(self):
        """
        Open a dialog to select a FASTA file, update parameters, and display the content with color-coded sequences.

        This method allows the user to select a FASTA file from the file system. It updates the relevant YAML parameters
        with the file's path and name, reads the file content, and displays the sequences with color-coded nucleotides
        in a text edit widget.

        Actions:
            - Opens a file dialog to select a FASTA file.
            - Updates 'reference_gene_file' and 'reference_gene_dir' parameters in the YAML file.
            - Reads and displays the content of the selected FASTA file.
            - Color-codes nucleotides (A, C, G, T) in the displayed sequence.
            - Enables the sequence alignment button and updates icons.
        """
        try:
            options = QFileDialog.Options()
            options |= QFileDialog.ReadOnly
            fullFileName, _ = QFileDialog.getOpenFileName(
                None,
                "Select FASTA file",
                "./datasets",
                "FASTA Files (*.fasta);;All Files (*)",
                options=options,
            )
            if fullFileName:
                update_yaml_param(
                    Params,
                    "scripts/utils/params.yaml",
                    "reference_gene_file",
                    os.path.basename(fullFileName),
                )
                update_yaml_param(
                    Params,
                    "scripts/utils/params.yaml",
                    "reference_gene_dir",
                    os.path.dirname(fullFileName),
                )

                with open(fullFileName, "r") as f:
                    self.clear_genetic_data()
                    content = f.read()
                    sequence = ""
                    for line in content.splitlines():
                        if line.startswith(">"):
                            species_title = line[1:].replace("_", " ")  # Replace underscores with whitespace
                            line = f'<span style="color: darkgreen; font-size: 24px;">{species_title}</span>'
                            sequence += "<br>" + line + "<br>"
                        else:
                            nucleotide_colors = {
                                "A": "green",
                                "C": "blue",
                                "G": "red",
                                "T": "black",
                            }
                            colored_line = ""
                            for char in line:
                                color = nucleotide_colors.get(char, "")
                                if color:
                                    colored_line += f'<span style="color: {color}; font-size: 20px;">{char}</span>'
                                else:
                                    colored_line += char
                            sequence += colored_line + " "

                    self.textEditFasta.setHtml(
                        f"<div style='background-color: #ffffff; color: #000000; padding: 10px; white-space: pre-wrap; word-wrap: break-word;'>{sequence}</div>"
                    )
                    self.sequenceAlignmentButtonPage1.setEnabled(True)
                    self.sequenceAlignmentButtonPage1.setIcon(QIcon(":inactive/sequence.svg"))
                    self.tabWidget.setCurrentIndex(1)
        except FileNotFoundError as e:
            self.show_error_dialog(f"File Not Found Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("StartSequenceAlignmentButton", QtEvents.clicked)
    def start_alignment_analysis(self):
        """
        Perform sequence alignment and store the resulting genetic tree dictionary.

        This method calls the `callSeqAlign` method to perform sequence alignment and stores the resulting genetic tree dictionary in the `geneticTreeDict` attribute.
        """
        self.starting_position_spinbox_2.setEnabled(True)
        self.window_size_spinbox_2.setEnabled(True)
        self.geneticTreeDict = self.call_seq_align()

    def call_seq_align(self):
        """
        Execute the sequence alignment pipeline and display progress using a worker thread.

        This method performs the following steps:
        1. Loads sequences from the reference gene file.
        2. Aligns the loaded sequences.
        3. Generates genetic trees based on the aligned sequences.
        4. Prepares and displays the results in the UI.
        5. Saves the alignment and genetic tree results to JSON files.

        Returns:
            None
        """

        def update_progress(step):
            if step < loading_screen.checkListWidget.count():
                item = loading_screen.checkListWidget.item(step)
                item.setCheckState(Qt.Checked)
                progress_value = int((step + 1) * (100 / loading_screen.checkListWidget.count()))
                loading_screen.progressBar.setValue(progress_value)
                QApplication.processEvents()
            else:
                loading_screen.progressBar.setValue(100)
                QApplication.processEvents()

        def handle_finished(result):
            loading_screen.close()
            self.msa = result["msa"]
            self.geneticTrees = result["geneticTrees"]
            self.geneticTreeButtonPage1.setEnabled(True)
            self.update_plot()
            self.statisticsButtonPage1.setEnabled(True)
            if self.climaticTreeButtonPage2.isEnabled():
                self.resultsButton.setEnabled(True)

        def handle_error(error_message):
            loading_screen.close()
            self.show_error_dialog(f"An unexpected error occurred: {error_message}")

        if loading_screen := loading.Ui_LoadingDialog():
            loading_screen.setupUi(loading_screen)
            # loading_screen.setWindowFlags(Qt.WindowType.FramelessWindowHint)  # Remove the title bar and frame

            # loading_screen.setWindowModality(Qt.WindowModality.ApplicationModal)

            # Set the QMovie for the movieLabel
            movie = QMovie(":active/dna.gif")  # Use the resource path for the gif
            loading_screen.movieLabel.setMovie(movie)

            # Resize the movie to fit within the QLabel
            movie.setScaledSize(QtCore.QSize(100, 100))  # Set the desired size here

            # Ensure the QLabel is centered and the GIF is properly displayed
            loading_screen.movieLabel.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            movie.start()

            # Show the loading screen
            loading_screen.show()

        QtWidgets.QApplication.processEvents()

        self.workerThread = QThread()
        self.worker = Worker(Params.reference_gene_filepath)
        self.worker.moveToThread(self.workerThread)

        self.worker.progress.connect(update_progress)
        self.worker.finished.connect(handle_finished)
        self.worker.error.connect(handle_error)

        self.workerThread.started.connect(self.worker.run)
        self.worker.finished.connect(self.workerThread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.workerThread.finished.connect(self.workerThread.deleteLater)

        self.workerThread.start()

        # Use a loop to wait until the thread finishes and the result is set
        while self.workerThread.isRunning():
            QApplication.processEvents()

        return self.geneticTrees

    # Pas certain que c'est utiliser...
    def stop_thread(self):
        if self.worker:
            self.worker.stop()
        if self.thread and self.thread.isRunning():
            self.thread.quit()
            self.thread.wait()

    # Pas certain que c'est utiliser...
    def closeEvent(self, event):
        self.stop_thread()
        event.accept()

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
            self.show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("fileBrowserButtonPage2", QtEvents.clicked)
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
                horizontal_header.setDefaultAlignment(Qt.AlignmentFlag.AlignHCenter)
                horizontal_header.setDefaultSectionSize(150)

            if vertical_header := table_widget.verticalHeader():
                vertical_header.setVisible(False)

            # Set headers
            for col in range(num_columns):
                item = QTableWidgetItem(df.columns[col])
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                table_widget.setHorizontalHeaderItem(col, item)

            # Fill the table with data
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
                self.statisticsButtonPage2.setEnabled(True)
                df = pd.read_csv(fullFilePath)
                columns = df.columns.tolist()

                if len(columns) < 2:
                    raise ValueError("The CSV file must contain at least two columns for latitude and longitude.")

                latitude_col = columns[-1]
                longitude_col = columns[-2]

                lat = df[latitude_col].tolist()
                long = df[longitude_col].tolist()

                # Replace underscores with spaces in species names (first column)
                df[columns[0]] = df[columns[0]].str.replace("_", " ")

                self.species = df[columns[0]].tolist()
                self.factors = df.drop(columns=[longitude_col, latitude_col]).values.tolist()
                clim_data_names = self.retrieve_data_names(columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "names", columns)
                update_yaml_param(Params, "scripts/utils/params.yaml", "data_names", clim_data_names)

                self.textEditClimData.clear()
                sleek_table = create_sleek_table(df)

                # Add the table widget to the textEditClimData layout
                layout = QVBoxLayout(self.textEditClimData)
                layout.addWidget(sleek_table)

                self.climaticTrees = utils.climaticPipeline(df)
                self.tree_keys = list(self.climaticTrees.keys())
                self.total_trees = len(self.tree_keys)
                self.current_index = 0
                self.climaticTreeButtonPage2.setEnabled(True)
                self.tabWidget2.setCurrentIndex(1)
                self.populate_map(lat, long)

                if self.statisticsButtonPage1.isEnabled():
                    self.resultsButton.setEnabled(True)
        except FileNotFoundError as e:
            self.show_error_dialog(f"File Not Found Error: {e}")
        except pd.errors.EmptyDataError as e:
            self.show_error_dialog(f"Empty Data Error: {e}")
        except ValueError as e:
            self.show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

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

            web_view = QWebEngineView(self.graphicsViewClimData)  # Embed the map inside graphicsViewClimData
            web_view.setHtml(data.getvalue().decode())
            layout = QVBoxLayout(self.graphicsViewClimData)
            layout.addWidget(web_view)
            self.graphicsViewClimData.setLayout(layout)

        except ValueError as e:
            self.show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("homeButton", QtEvents.clicked)
    def show_home_section(self):
        """
        Display the home page of the application.

        This method sets the icons for the climatic data and genetic data buttons to their inactive states
        and displays the home page by setting the stacked widget's current index to 0.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.homeButton.setIcon(QIcon(":active/home.png"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.stackedWidget.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("geneticDataButton", QtEvents.clicked)
    def show_genetic_section(self):
        """
        Display the genetic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the genetic data page
        by setting the stacked widget's current index to 1, and sets the tab widget's current index to 0.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.geneticDataButton.setIcon(QIcon(":active/genetic.svg"))
            self.homeButton.setIcon(QIcon(":other/home.svg"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.stackedWidget.setCurrentIndex(1)
            self.tabWidget.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("climaticDataButton", QtEvents.clicked)
    def show_climate_section(self):
        """
        Display the climatic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the climatic data page
        by setting the stacked widget's current index to 2, and sets the tab widget's current index to 0.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":active/climaticData.png"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.homeButton.setIcon(QIcon(":other/home.svg"))
            self.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.stackedWidget.setCurrentIndex(2)
            self.tabWidget2.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("resultsButton", QtEvents.clicked)
    def show_results_section(self):
        """
        Display the results page of the application.

        This method sets the stacked widget's current index to 3 to display the results page.
        """
        try:
            self.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.homeButton.setIcon(QIcon(":other/home.svg"))
            self.resultsButton.setIcon(QIcon(":active/result.svg"))
            self.stackedWidget.setCurrentIndex(3)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("sequenceAlignmentButtonPage1", QtEvents.clicked)
    def show_sequence_alignment_page(self):
        """
        Display the sequence alignment page.

        This method sets the stacked widget's current index to 1 and the tab widget's current index to 2
        to display the sequence alignment page.
        """
        try:
            self.stackedWidget.setCurrentIndex(1)
            self.tabWidget.setCurrentIndex(2)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("submitButtonPage3", QtEvents.clicked)
    def show_filtered_results(self):
        """
        Show the results filtered with a metric threshold provided by the user.

        This method reads the data from a CSV file, processes it through the climatic pipeline, filters the results,
        and displays the filtered results in an HTML table format within a QTextBrowser widget.
        It handles exceptions related to missing sequence alignment.

        Raises:
            AttributeError: If the sequence alignment has not been performed before attempting to generate the tree.
        """
        self.stackedWidget.setCurrentIndex(3)

        df = pd.read_csv(Params.file_name)
        utils.filterResults(self.climaticTrees, self.geneticTreeDict, df)
        df_results = pd.read_csv("./scripts/results/output.csv")
        df_results["Name of species"] = df_results["Name of species"].str.replace("_", " ")
        # Replace the first column values with Params.file_name just before visualization
        df_results.iloc[:, 0] = Params.reference_gene_file

        # Convert to HTML table with sleek, modern styling
        html_table = df_results.to_html(index=False, border=0, classes="dataframe", table_id="styled-table")

        # Add CSS styles for a professional look with smooth hover animations
        html_style = """
        <style>
            #styled-table {
                font-family: 'Trebuchet MS', Arial, Helvetica, sans-serif;
                border-collapse: collapse;
                width: 100%;
                border-radius: 5px;
                box-shadow: 0 2px 5px rgba(0, 0, 0, 0.15);
                overflow: hidden;
            }
            #styled-table td, #styled-table th {
                border: 1px solid #ddd;
                padding: 12px;
                transition: all 0.3s ease-in-out;
            }
            #styled-table tr:nth-child(even) {
                background-color: #f9f9f9;
            }
            #styled-table tr:hover {
                background-color: #f1f1f1;
                transform: scale(1.01);
            }
            #styled-table th {
                padding-top: 12px;
                padding-bottom: 12px;
                text-align: left;
                background-color: #4CAF50;
                color: white;
            }
            #styled-table td {
                padding-left: 12px;
                padding-right: 12px;
            }
        </style>
        """

        # Combine the HTML style and table
        html_content = html_style + html_table

        # Set up the QWebEngineView
        web_engine_view = QWebEngineView()
        web_engine_view.setHtml(html_content)

        # Clear the existing content and layout of the QTextBrowser (if any)
        layout = QVBoxLayout(self.textEditResults)
        for i in reversed(range(layout.count())):
            widget_to_remove = layout.itemAt(i).widget()
            layout.removeWidget(widget_to_remove)
            widget_to_remove.setParent(None)

        # Add the QWebEngineView to the QTextBrowser
        layout.addWidget(web_engine_view)

    @connect_event("darkModeButton", QtEvents.clicked)
    def toggle_dark_mode(self):
        """
        Toggle the application's dark mode setting.

        This method switches the application's theme between dark mode and light mode. It updates the isDarkMode attribute,
        applies the corresponding style, and changes the icon of the darkModeButton.

        Attributes:
            isDarkMode (bool): A flag indicating whether dark mode is currently enabled.

        Actions:
            If dark mode is enabled, apply the dark style and set the darkModeButton icon to the 'light' icon.
            If dark mode is disabled, apply the light style and set the darkModeButton icon to the 'dark' icon.
        """
        try:
            self.isDarkMode = not self.isDarkMode
            buttons_vertical = [
                self.fileBrowserButtonPage1,
                self.sequenceAlignmentButtonPage1,
                self.clearButtonPage1,
                self.statisticsButtonPage1,
                self.geneticTreeButtonPage1,
                self.fileBrowserButtonPage2,
                self.clearButtonPage2,
                self.climaticTreeButtonPage2,
                self.statisticsButtonPage2,
                self.settingsButtonPage3,
                self.settingsButtonPage4,
                self.submitButtonPage3,
                self.statisticsButtonPage3,
                self.submitButtonPage4,
                self.statisticsButtonPage4,
                self.clearButtonPage3,
                self.StartSequenceAlignmentButton,
                self.clearButtonPage4,
            ]
            buttons = [
                self.geneticDataButton,
                self.climaticDataButton,
                self.helpButton,
                self.homeButton,
                self.resultsButton,
            ]

            style_buttons(buttons, self.isDarkMode)
            style_buttons(buttons_vertical, self.isDarkMode)

            if self.isDarkMode:
                qtmodern.styles.dark(app)
                self.top_frame.setStyleSheet("background-color: #646464;")
                self.darkModeButton.setIcon(QIcon(":other/light.png"))  # Set the 'light' icon for dark mode
            else:
                qtmodern.styles.light(app)
                self.top_frame.setStyleSheet("background-color: rgb(222, 221, 218);")
                self.darkModeButton.setIcon(QIcon(":other/dark.png"))  # Set the 'dark' icon

            # Common settings for both modes
            self.darkModeButton.setCursor(Qt.PointingHandCursor)
            self.darkModeButton.setStyleSheet(get_button_style(self.isDarkMode))
            self.darkModeButton.setGraphicsEffect(create_shadow_effect(10, 140))

        except AttributeError as e:
            self.show_error_dialog(f"An error occurred while setting attributes: {e}", "Attribute Error")

        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("clearButtonPage1", QtEvents.clicked)
    def clear_genetic_data(self):
        """
        Clear the genetic data fields.

        This method clears the content of the text edit widgets related to FASTA sequences, sequence alignments, and genetic trees.
        It also resets the current index of the genetic statistics list to 0 and disables relevant buttons.
        """
        try:
            self.textEditFasta.clear()
            self.seqAlignLabel.clear()
            self.textEditGenStats_2.clear()
            self.sequenceAlignmentButtonPage1.setEnabled(False)
            self.statisticsButtonPage1.setEnabled(False)
            self.geneticTreeButtonPage1.setEnabled(False)
            self.GeneticTreeLabel.clear()
            self.resultsButton.setEnabled(False)
            self.geneticTreeDict = None
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("clearButtonPage2", QtEvents.clicked)
    def clear_climmatic_data(self):
        """
        Clear the text fields related to climatic data.

        This method disables buttons related to climatic data and clears the necessary fields.
        """
        try:
            self.statisticsButtonPage2.setEnabled(False)
            self.climaticTreeButtonPage2.setEnabled(False)
            self.resultsButton.setEnabled(False)
            self.climaticTrees = None
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}", "Error")

    # elle se fait override par celle qui suit, regarder si elle est utile au final (penses pas)
    @connect_event("clearButtonPage4", QtEvents.clicked)
    def clear_result(self):
        """
        Clear the text fields related to climatic data.

        This method clears the content of the text edit widgets for climatic data, climatic statistics, and climatic trees.
        It also clears the graphics view for climatic data and resets the current index of the climatic statistics and chart lists to 0.
        """
        try:
            self.textEditResults.clear()
            self.textEditClimStats.clear()
            self.textEditClimTrees.clear()
            self.graphicsViewClimData.clear()
            self.ClimStatsListCondition.setCurrentIndex(0)
            self.ClimStatsListChart.setCurrentIndex(0)
            self.ResultsStatsListCondition.setCurrentIndex(0)
            self.ResultsStatsListChart.setCurrentIndex(0)
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}", "Error")

    @connect_event("geneticTreeButtonPage1", QtEvents.clicked)
    def display_newick_trees(self):
        """
        Display Newick format trees in the application using Toytree.

        This method loads Newick format trees from a JSON file, formats the tree names,
        and updates the UI to display the trees using Toytree.

        Actions:
            - Loads Newick format trees from 'scripts/results/geneticTrees.json'.
            - Formats the tree names by replacing underscores with ' nt '.
            - Updates the tree combo box with formatted tree names.
            - Displays the first tree in the list.
        """

        self.tabWidget.setCurrentIndex(4)
        file_path = "scripts/results/geneticTrees.json"
        with open(file_path, "r") as file:
            self.newick_json = json.load(file)

        self.tree_keys = list(self.newick_json.keys())
        self.total_trees = len(self.tree_keys)
        self.current_index = 0
        self.geneticTreescomboBox.clear()

        # Format the tree keys to replace underscore with ' nt '
        formatted_tree_keys = [self.format_tree_name(key) for key in self.tree_keys]
        self.geneticTreescomboBox.addItems(formatted_tree_keys)

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

    @connect_event("geneticTreescomboBox", QtEvents.currentIndexChanged)
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
        canvas, axes, mark = tree.draw(
            width=921,
            height=450,
            tip_labels=tip_labels,  # These labels now have spaces
            tip_labels_style={"font-size": "15px"},
            fixed_order=tip_labels,
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
        self.GeneticTreeLabel.clear()
        self.GeneticTreeLabel.setPixmap(pixmap)
        self.GeneticTreeLabel.adjustSize()

    @connect_event("downloadGraphButton", QtEvents.clicked)
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
            default_file_name = f"{current_key}.png"

            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            file_path, _ = QFileDialog.getSaveFileName(
                self,
                "Save Tree Image As",
                default_file_name,
                "PNG Files (*.png);;All Files (*)",
                options=options,
            )
            if file_path:
                if not file_path.lower().endswith(".png"):
                    file_path += ".png"
                shutil.copy(self.tree_img_path, file_path)
        except FileNotFoundError as e:
            self.show_error_dialog(f"The tree image file was not found: {e}", "File Not Found")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while downloading the tree image: {e}")

    @connect_event("preferencesButton", QtEvents.clicked)
    def open_climatic_tree_preferences_window(self):
        """
        Open the preferences dialog and update the application settings based on user input.

        This method opens a PreferencesDialog, updates it with the current preferences, and applies the new preferences if the user accepts the changes.

        Returns:
            None
        """
        try:
            dialog = PreferencesDialog(self)
            dialog.update_preferences(self.preferences)
            if dialog.exec() == QDialog.Accepted:
                self.preferences = dialog.get_preferences()
                self.apply_preferences()
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while opening the preferences dialog: {e}")

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
            self.show_error_dialog(f"An unexpected error occurred while applying preferences: {e}")

    @connect_event("climaticTreeButtonPage2", QtEvents.clicked)
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
            self.climaticTreescomboBox.clear()
            self.tree_keys = list(self.climaticTrees.keys())
            self.total_trees = len(self.tree_keys)
            self.current_index = 0
            self.climaticTreescomboBox.addItems(self.tree_keys)
            self.show_climatic_tree(self.current_index)
            self.tabWidget2.setCurrentIndex(3)
        except KeyError as e:
            self.show_error_dialog(
                f"An error occurred while accessing the climatic trees: {e}",
                "Key Error",
            )
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    @connect_event("climaticTreescomboBox", QtEvents.currentIndexChanged)
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
            self.show_error_dialog(f"An unexpected error occurred while displaying the selected climatic tree: {e}")

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
                preferences = self.preferences
                label_color = preferences.get("label_color", "black")
                edge_color = preferences.get("edge_color", "blue")
                reticulation_color = preferences.get("reticulation_color", "red")
                layout = preferences.get("layout", "horizontal")
                proportional_edge_lengths = preferences.get("proportional_edge_lengths", False)
                label_internal_vertices = preferences.get("label_internal_vertices", False)
                use_leaf_names = preferences.get("use_leaf_names", True)
                show_branch_length = preferences.get("show_branch_length", False)
                view_type = preferences.get("view_type", "network")

                if view_type == "network":
                    self.render_network_view(
                        tree,
                        label_color,
                        edge_color,
                        reticulation_color,
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
                        edge_color,
                        reticulation_color,
                        layout,
                        proportional_edge_lengths,
                        label_internal_vertices,
                        use_leaf_names,
                        show_branch_length,
                    )
        except KeyError as e:
            self.show_error_dialog(
                f"An error occurred while accessing the climatic tree data: {e}",
                "Key Error",
            )
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred: {e}")

    def render_network_view(
        self,
        tree,
        label_color,
        edge_color,
        reticulation_color,
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
            reticulation_color (str): Color for reticulation edges.
            layout (str): Layout for the network visualization (e.g., 'horizontal').
            proportional_edge_lengths (bool): Whether to use proportional edge lengths.
            label_internal_vertices (bool): Whether to label internal vertices.
            use_leaf_names (bool): Whether to use leaf names.
            show_branch_length (bool): Whether to show branch lengths.

        Returns:
            None
        """
        try:
            # Use seaborn for color palettes
            sns.set_palette("husl")

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

            # Render the tree using Plotly for interactive visualization
            graph = Phylo.to_networkx(tree)
            pos = self.get_layout(graph, layout)
            edge_trace, edge_annotations = self.create_edge_trace(tree, pos, edge_color, show_branch_length)
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
            self.climaticTreesLabel.clear()
            self.climaticTreesLabel.setPixmap(pixmap)
            self.climaticTreesLabel.adjustSize()
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while rendering the network view: {e}")

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
                label_colors={clade: label_color for clade in tree.find_clades()},
            )

            ax.axis("off")  # Remove axes

            results_dir = "results"
            os.makedirs(results_dir, exist_ok=True)
            img_path = os.path.join(results_dir, f"{self.tree_keys[self.current_index]}.png")
            plt.savefig(img_path, format="png")
            plt.close(fig)

            pixmap = QPixmap(img_path)
            self.climaticTreesLabel.clear()
            self.climaticTreesLabel.setPixmap(pixmap)
            self.climaticTreesLabel.adjustSize()
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while rendering the tree view: {e}")

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
            self.show_error_dialog(f"An unexpected error occurred while creating the node trace: {e}")

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
            self.show_error_dialog(f"An unexpected error occurred while getting the layout: {e}")

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
            self.show_error_dialog(f"An unexpected error occurred while creating the edge trace: {e}")

    @connect_event("downloadGraphButton2", QtEvents.clicked)
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

            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            file_path, _ = QFileDialog.getSaveFileName(
                self,
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
            self.show_error_dialog(f"The temporary image file was not found: {e}", "File Not Found")
        except Exception as e:
            self.show_error_dialog(f"An unexpected error occurred while downloading the climatic tree graph: {e}")

    ################################################
    @connect_event(["statisticsButtonPage3", "statisticsButtonPage4"], QtEvents.clicked)
    def display_phylogeographic_trees(self):
        self.stackedWidget.setCurrentIndex(4)
        self.results_dir = "scripts/results"
        file_path = os.path.join(self.results_dir, "geneticTrees.json")

        with open(file_path, "r") as file:
            self.phylo_json = json.load(file)

        self.tree_keys = list(self.phylo_json.keys())
        self.total_trees = len(self.tree_keys)

        self.current_index1 = 0
        self.phyloTreescomboBox.clear()

        formatted_tree_keys = [self.rename_tree_key(key) for key in self.tree_keys]
        self.phyloTreescomboBox.addItems(formatted_tree_keys)

        self.climatic_data = pd.read_csv(Params.file_name)
        self.criteriaComboBox.clear()
        self.criteriaComboBox.addItems(self.climatic_data.columns[1:])

        self.render_tree(self.current_index1)

    def rename_tree_key(self, tree_key):
        parts = tree_key.split("_")
        if len(parts) == 2:
            return f"{parts[0]} nt {parts[1]} nt"
        return tree_key

    @connect_event(["phyloTreescomboBox", "criteriaComboBox"], QtEvents.currentIndexChanged)
    def render_tree(self, index):
        if 0 <= index < self.total_trees:
            self.current_index1 = index
            key = self.tree_keys[index]
            newick_str = self.phylo_json[key]

            tree = toytree.tree(newick_str)
            tip_labels = tree.get_tip_labels()

            if not hasattr(self, "climatic_data"):
                return

            selected_column = self.criteriaComboBox.currentText()

            # Match the original labels from the CSV
            ordered_climatic_data = self.climatic_data.set_index("id").reindex(tip_labels).reset_index()
            bar_values = ordered_climatic_data[selected_column].values

            # Ensure bar_values have no negative values
            bar_values = np.clip(bar_values, a_min=0, a_max=None)

            # Create a canvas for the plot
            canvas = toyplot.Canvas(width=921, height=450)

            # Add tree to canvas
            ax0 = canvas.cartesian(bounds=(50, 300, 50, 400), ymin=0, ymax=tree.ntips, padding=15)

            # Replace underscores with spaces for display purposes only
            tree.draw(axes=ax0, tip_labels=[label.replace("_", " ") for label in tip_labels])
            ax0.show = False

            # Add bar plot to canvas
            ax1 = canvas.cartesian(bounds=(325, 900, 50, 400), ymin=0, ymax=tree.ntips, padding=15)
            ax1.bars(
                np.arange(tree.ntips),
                bar_values,
                along="y",
            )

            # Adjust the appearance for better visibility
            ax1.show = True
            ax1.y.show = False
            ax1.x.ticks.show = True
            ax1.x.label.text = selected_column

            # Save and display the plot
            self.temp_img_path = os.path.join(self.results_dir, f"{key}.png")
            toyplot.png.render(canvas, self.temp_img_path)

            pixmap = QPixmap(self.temp_img_path)
            self.PhyloTreeLabel.clear()
            self.PhyloTreeLabel.setPixmap(pixmap)
            self.PhyloTreeLabel.adjustSize()

    @connect_event("downloadResultsPlotButton", QtEvents.clicked)
    def save_tree_graph(self):
        current_key = self.tree_keys[self.current_index1]
        default_file_name = f"{current_key}.png"

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Graph As",
            default_file_name,
            "PNG Files (*.png);;All Files (*)",
            options=options,
        )

        if file_path:
            if not file_path.lower().endswith(".png"):
                file_path += ".png"

            with open(self.temp_img_path, "rb") as temp_file:
                with open(file_path, "wb") as file:
                    file.write(temp_file.read())


if __name__ == "__main__":
    # Create the application instance
    app = QtWidgets.QApplication([])

    # Apply modern styles
    qtmodern.styles.light(app)

    # Create the main window instance
    window = UiMainWindow()

    # Wrap the main window with the ModernWindow style
    mw = qtmodern.windows.ModernWindow(window)

    # Get screen geometry to determine the available screen space
    screen_geometry = app.primaryScreen().availableGeometry()

    # Calculate the center position of the screen
    center_point = screen_geometry.center()
    x = center_point.x() - mw.width() // 2
    y = center_point.y() - mw.height() // 2

    # Move the main window to the center of the screen
    mw.move(x, y)

    # Show the main window
    mw.show()

    # Execute the application's event loop
    sys.exit(app.exec())
