import io
import toytree
import tempfile
import json
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtGui import QPixmap
import matplotlib.pyplot as plt

import tempfile
import json
import toytree
import matplotlib.pyplot as plt
from PyQt5.QtGui import QPixmap
import json
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from PyQt5.QtWidgets import QVBoxLayout, QTextBrowser
from PyQt5.QtWebEngineWidgets import QWebEngineView
import sys
import toyplot.png
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from PyQt5.QtGui import QPixmap
from io import BytesIO
import tempfile
from collections import Counter
from decimal import Decimal
from io import BytesIO
import toyplot
import toytree
import resources_rc
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
import qtmodern.windows
import seaborn as sns
import toytree
import yaml
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QThread, Qt
from PyQt5.QtCore import pyqtSignal, QObject
from PyQt5.QtGui import QIcon, QColor, QPixmap, QImage
from PyQt5.QtGui import QMovie
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import (QDialog)
from PyQt5.QtWidgets import QFileDialog, QGraphicsDropShadowEffect
from aphylogeo import utils
from aphylogeo.alignement import AlignSequences
from aphylogeo.genetic_trees import GeneticTrees
from aphylogeo.params import Params

from PreferencesDialog import PreferencesDialog  # Import PreferencesDialog
from help import UiHowToUse
from settings import HoverLabel

Params.load_from_file("params.yaml")


class Worker(QObject):
    progress = pyqtSignal(int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, filepath):
        super().__init__()
        self.filepath = filepath

    def run(self):
        try:
            # Step 1: Load sequences
            self.progress.emit(0)
            sequenceFile = utils.loadSequenceFile(self.filepath)

            # Step 2: Align sequences
            self.progress.emit(1)
            align_sequence = AlignSequences(sequenceFile)
            alignments = align_sequence.align()

            # Step 3: Generate genetic trees
            self.progress.emit(2)
            geneticTrees = utils.geneticPipeline(alignments.msa)
            trees = GeneticTrees(trees_dict=geneticTrees, format="newick")

            # Step 4: Preparing results
            msa = alignments.to_dict().get("msa")

            # Step 5: Save results
            alignments.save_to_json(f"./results/aligned_{Params.reference_gene_file}.json")
            trees.save_trees_to_json("./results/geneticTrees.json")

            # Emit finished signal with the genetic trees dictionary
            result = {"msa": msa, "geneticTrees": geneticTrees}
            self.progress.emit(3)
            self.finished.emit(result)

        except Exception as e:
            self.error.emit(str(e))


class MyDumper(yaml.Dumper):
    """
     Custom YAML Dumper to modify the default indentation and list representation behavior.

     Methods:
         increase_indent(flow=False, indentless=False):
             Increase the indentation level in the YAML output.

         represent_list(data):
             Represent Python lists in a flow style in the YAML output.
     """

    def increase_indent(self, flow=False, indentless=False):
        """
        Increase the indentation level in the YAML output.

        Args:
            flow (bool): Indicates whether the current context is a flow style. Defaults to False.
            indentless (bool): Indicates whether to use an indentless format. This argument is ignored. Defaults to False.

        Returns:
            The result from the superclass's increase_indent method with modified behavior.
        """
        return super(MyDumper, self).increase_indent(flow, False)

    def represent_list(self, data):
        """
        Represent Python lists in a flow style in the YAML output.

        Args:
            data (list): The list to represent in the YAML output.

        Returns:
            The YAML representation of the list in a flow style.
        """
        return self.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.add_representer(list, MyDumper.represent_list, Dumper=MyDumper)
"""
Add a representer for the list type to the YAML dumper.

This ensures that lists are represented in a flow style using the MyDumper class.

Args:
    list (type): The Python list type to represent.
    MyDumper.represent_list (method): The method that defines how to represent lists.
    Dumper (yaml.Dumper): The custom dumper class to use, in this case, MyDumper.
"""


def update_yaml_param(params, file_path, property_name, new_value):
    """
    Updates a specified property within a YAML file with a new value.

    Args:
        params: An object with an update_from_dict method, typically used for updating parameters.
        file_path (str): The path to the YAML file.
        property_name (str): The name of the property to modify (e.g., 'file_name').
        new_value: The new value to set for the property (can be any valid YAML type).

    Raises:
        FileNotFoundError: If the specified YAML file does not exist.
        KeyError: If the specified property name is not found in the YAML file.
    """
    if isinstance(new_value, list):
        new_value = [element.strip() for element in new_value]
    params.update_from_dict({property_name: new_value})

    # 1. Load existing YAML data
    with open(file_path, "r") as yaml_file:
        data = yaml.safe_load(yaml_file)  # Use safe_load for security

    # 2. Update the specified property
    if property_name in data:
        data[property_name] = new_value
    # 3. Write the updated data back to the file
    with open(file_path, "w") as yaml_file:
        yaml.dump(data, yaml_file, default_flow_style=None, Dumper=MyDumper, sort_keys=False)


window_size = 50
starting_position = 1


class UiMainWindow(QtWidgets.QMainWindow):

    def useWindow(self):
        """
        Initialize and display the 'How to Use' window.

        This method creates a new QMainWindow instance, sets up its UI using the UiHowToUse class, and displays the window.
        """
        self.window = QtWidgets.QMainWindow()
        self.ui = UiHowToUse()
        self.ui.initUI()
        self.ui.show()

    def paramWin(self):
        """
        Initialize and display the parameters window.
        This method creates a new QMainWindow instance, sets up its UI using the UiDialog class, and displays the window.
        """
        dialog = QtWidgets.QDialog()
        ui = HoverLabel.Settings()
        ui.setupUi(dialog)
        dialog.exec_()

    def openClimTree(self):
        """
        Initialize and display the climatic tree window.

        This method imports the Ui_ct class, creates a new QMainWindow instance,
        sets up its UI using the Ui_ct class, and displays the window.
        It also sets the current index for stackedWidget and tabWidget2.
        """
        try:
            from cltree import Ui_ct

            self.window = QtWidgets.QMainWindow()
            self.ui = Ui_ct()
            self.ui.setupUi(self.window)
            self.window.show()

            self.stackedWidget.setCurrentIndex(2)
            self.tabWidget2.setCurrentIndex(3)

        except ImportError as e:
            self.showErrorDialog(f"An error occurred while importing: {e}", "Import Error")
        except AttributeError as e:
            self.showErrorDialog(f"An error occurred while setting attributes: {e}", "Attribute Error")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}", "Unexpected Error")

    def showErrorDialog(self, message, title="error"):
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
        msgBox.exec_()

    def __init__(self):
        super(UiMainWindow, self).__init__()
        uic.loadUi("Qt/main.ui", self)
        self.setupUi()

    def setupUi(self):
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
                "view_type": "network"
            }

            self.setObjectName("MainWindow")
            self.window_size_spinbox_2.setRange(1, 1000)
            self.starting_position_spinbox_2.setRange(1, 1000)
            self.starting_position_spinbox_2.valueChanged.connect(self.update_plot)
            self.window_size_spinbox_2.valueChanged.connect(self.update_plot)
            self.homeButton.clicked.connect(self.showHomePage)
            self.geneticDataButton.clicked.connect(self.showGenDatPage)
            self.clearButtonPage3.clicked.connect(self.clearResults)
            self.climaticDataButton.clicked.connect(self.showClimDatPage)
            self.helpButton.clicked.connect(self.useWindow)
            self.darkModeButton.clicked.connect(self.toggleDarkMode)
            self.climaticTreeButtonPage2.clicked.connect(self.displayClimaticTrees)
            self.climaticTreescomboBox.currentIndexChanged.connect(self.show_selected_climatic_tree)
            self.downloadGraphButton2.clicked.connect(self.download_climatic_tree_graph)
            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
            self.isDarkMode = False  # Keep track of the state
            self.fileBrowserButtonPage1.clicked.connect(self.pressItFasta)
            self.geneticTreeButtonPage1.clicked.connect(self.display_newick_trees)
            self.sequenceAlignmentButtonPage1.clicked.connect(self.showSequencePage)
            self.clearButtonPage1.clicked.connect(self.clearGen)
            self.statisticsButtonPage1.clicked.connect(self.plot_sequence_similarity)

            self.clearButtonPage2.clicked.connect(self.clearClim)
            self.fileBrowserButtonPage2.clicked.connect(self.pressItCSV)
            self.statisticsButtonPage2.clicked.connect(self.load_data_climate)
            self.resultsButton.clicked.connect(self.showResultsPage)
            self.ClimaticChartSettingsAxisX.currentIndexChanged.connect(self.generate_graph)
            self.ClimaticChartSettingsAxisY.currentIndexChanged.connect(self.generate_graph)
            self.radioButtonPiePlot.toggled.connect(self.generate_graph)
            self.radioButtonLinePlot.toggled.connect(self.generate_graph)
            self.radioButtonScatterPlot.toggled.connect(self.generate_graph)
            self.radioButtonBarGraph.toggled.connect(self.generate_graph)
            self.geneticTreescomboBox.currentIndexChanged.connect(self.show_selected_tree)
            self.StartSequenceAlignmentButton.clicked.connect(self.SeqAlign)
            self.settingsButtonPage3.clicked.connect(self.paramWin)
            self.submitButtonPage3.clicked.connect(self.showFilteredResults)
            self.clearButtonPage4.clicked.connect(self.clearResult)
            self.statisticsButtonPage4.clicked.connect(self.showResultsStatsPage)
            self.clearButtonPage4.clicked.connect(self.clearResultStat)
            self.downloadGraphButton.clicked.connect(self.download_graph)
            self.preferencesButton.clicked.connect(self.open_preferences_dialog)
            self.stackedWidget.setCurrentIndex(0)
            self.tree_keys = []
            self.total_trees = 0
            self.current_index = 0

            buttons = [self.geneticDataButton, self.climaticDataButton, self.helpButton, self.homeButton,
                       self.resultsButton]
            buttons_Vertical = [
                self.fileBrowserButtonPage1, self.sequenceAlignmentButtonPage1, self.clearButtonPage1,
                self.statisticsButtonPage1, self.geneticTreeButtonPage1, self.fileBrowserButtonPage2,
                self.clearButtonPage2, self.climaticTreeButtonPage2, self.statisticsButtonPage2,
                self.settingsButtonPage3, self.settingsButtonPage4, self.submitButtonPage3,
                self.statisticsButtonPage3, self.submitButtonPage4, self.statisticsButtonPage4,
                self.StartSequenceAlignmentButton, self.clearButtonPage3, self.clearButtonPage4
            ]

            # Define cursor and stylesheet for all buttons
            for button in buttons:
                button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
                button.setStyleSheet("""
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
                """)
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 140))
                shadow_effect.setOffset(3, 3)
                button.setGraphicsEffect(shadow_effect)

            for buttonV in buttons_Vertical:
                buttonV.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
                buttonV.setStyleSheet("""
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
                """)
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 110))
                shadow_effect.setOffset(3, 3)
                buttonV.setGraphicsEffect(shadow_effect)

            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
            self.darkModeButton.setStyleSheet("""
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
            """)
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)
            shadow_effect.setColor(QColor(0, 0, 0, 140))
            shadow_effect.setOffset(3, 3)
            self.darkModeButton.setGraphicsEffect(shadow_effect)

            QtCore.QMetaObject.connectSlotsByName(self)

        except AttributeError as e:
            self.showErrorDialog(f"An error occurred while setting up the UI: {e}", "Attribute Error")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}", "Unexpected Error", )


    def clearResults(self):
        self.textEditResults.clear()



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
                lines = value.strip().split('\n')
                current_id = None
                for line in lines:
                    if line.startswith('>'):
                        current_id = line[1:].strip()
                        if current_id not in genetic_data:
                            genetic_data[current_id] = []
                    else:
                        genetic_data[current_id].append(line.strip())

            genetic_data = {sequence_id: ''.join(sequences) for sequence_id, sequences in genetic_data.items()}
            return genetic_data

        except KeyError as e:
            self.showErrorDialog(f"Key Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

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
            standardized_data = {key: seq.ljust(max_length, '-') for key, seq in genetic_data.items()}
            return standardized_data

        except ValueError as e:
            self.showErrorDialog(f"Value Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

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
            end_position = starting_position + window_size
            truncated_data = {key: value[starting_position:end_position] for key, value in genetic_data.items()}
            alignment = MultipleSeqAlignment([
                SeqRecord(Seq(seq), id=key)
                for key, seq in truncated_data.items()
            ])

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 4), gridspec_kw={'height_ratios': [1, 8]})
            ax1.set_axis_off()
            ax2.set_axis_off()

            # Calculate conservation
            conservation, _ = self.calculate_conservation_and_gaps(alignment)

            # Plot conservation
            bar_width = 1.0
            ax1.bar(range(len(conservation)), conservation, color='#4CAF50', width=bar_width, align='edge')
            ax1.set_xlim(0, len(conservation))
            ax1.set_ylim(0, 1)
            ax1.set_title('CONSERVATION', fontsize=12, pad=20)

            # Plot alignment
            seqs = [str(record.seq) for record in alignment]
            ids = [record.id for record in alignment]
            colors = {
                'A': '#4CAF50',
                'T': '#F44336',
                'G': '#2196F3',
                'C': '#FFEB3B',
                '-': 'grey',
                'N': 'black'  # Handling 'N' or any other characters if present
            }

            font_size = 10
            rect_height = 0.8
            rect_width = 1.0

            for i, seq in enumerate(seqs):
                ax2.text(-1, len(seqs) - i - 1 + rect_height / 2, ids[i], ha='right', va='center', fontsize=font_size)
                for j, nucleotide in enumerate(seq):
                    color = colors.get(nucleotide, 'black')  # Default to black if not found
                    rect = mpatches.Rectangle((j, len(seqs) - i - 1), rect_width, rect_height, color=color)
                    ax2.add_patch(rect)
                    ax2.text(j + rect_width / 2, len(seqs) - i - 1 + rect_height / 2, nucleotide, ha='center',
                             va='center', fontsize=font_size)

            consensus_seq = self.calculate_consensus(alignment)
            for j, nucleotide in enumerate(consensus_seq):
                color = colors.get(nucleotide, 'black')  # Default to black if not found
                rect = mpatches.Rectangle((j, -1), rect_width, rect_height, color=color)
                ax2.add_patch(rect)
                ax2.text(j + rect_width / 2, -1 + rect_height / 2, nucleotide, ha='center', va='center',
                         fontsize=font_size, fontweight='bold')

            # Add Consensus title
            ax2.text(-1, -1 + rect_height / 2, "Consensus", ha='right', va='center', fontsize=font_size,
                     fontweight='bold')

            ax2.set_xlim(0, len(consensus_seq))
            ax2.set_ylim(-2, len(seqs))
            ax2.set_yticks(range(len(ids)))
            ax2.set_yticklabels(ids, fontsize=font_size)

            # Ensure number of ticks matches the length of the sequences
            ax2.set_xticks(range(len(consensus_seq)))
            ax2.set_xticklabels(range(starting_position, starting_position + len(consensus_seq)), fontsize=font_size)

            plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.05)
            plt.savefig(output_path)
            plt.close()

        except KeyError as e:
            self.showErrorDialog(f"Key Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

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
                gaps.append(counter['-'] / len(column))

            return conservation, gaps

        except KeyError as e:
            self.showErrorDialog(f"Key Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

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

            return ''.join(consensus)

        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

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
            output_path = "sequence_alignment_plot.png"

            genetic_data = self.read_msa(self.msa)
            standardized_data = self.standardize_sequence_lengths(genetic_data)
            self.plot_alignment_chart(standardized_data, starting_position, window_size, output_path)

            pixmap = QPixmap(output_path)
            self.seqAlignLabel.setPixmap(pixmap)
            self.tabWidget.setCurrentIndex(2)

        except AttributeError as e:
            self.showErrorDialog(f"Attribute Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def plot_sequence_similarity(self):
        """
        Plot the sequence similarity based on a multiple sequence alignment (MSA).

        This method reads the MSA data, combines sequences for each species, pads sequences to the same length, computes similarity scores,
        and plots the similarity across the alignment. The plot is then displayed in the specified label widget.

        Returns:
            None
        """
        try:
            from collections import defaultdict

            # Dictionary to hold combined sequences for each species
            sequences = defaultdict(str)

            # Combine sequences across all ranges for each species
            for key, value in self.msa.items():
                parts = value.strip().split('\n')
                for i in range(0, len(parts), 2):
                    header = parts[i].strip('>')
                    sequence = parts[i + 1]
                    sequences[header] += sequence

            # Find the maximum sequence length
            max_len = max(len(seq) for seq in sequences.values())

            # Pad sequences to the same length
            padded_records = []
            for header, sequence in sequences.items():
                padded_seq = sequence.ljust(max_len, '-')
                padded_records.append(SeqRecord(Seq(padded_seq), id=header))

            # Create a MultipleSeqAlignment object
            alignment = MultipleSeqAlignment(padded_records)

            # Compute the similarity score for each position
            reference_sequence = str(alignment[0].seq)
            similarities = []

            for record in alignment:
                similarity = [1 if ref == res else 0 for ref, res in zip(reference_sequence, str(record.seq))]
                similarities.append(similarity)

            similarities = np.array(similarities)

            # Compute sliding window averages
            def sliding_window_avg(arr, window_size, step_size):
                return [np.mean(arr[i:i + window_size]) for i in range(0, len(arr) - window_size + 1, step_size)]

            # Use Params.window_size directly
            window_size = Params.window_size
            step_size = Params.step_size  # You can also parameterize this if needed

            windowed_similarities = []
            for sim in similarities:
                windowed_similarities.append(sliding_window_avg(sim, window_size, step_size))

            windowed_similarities = np.array(windowed_similarities)

            # Plot the similarities with advanced styling
            sns.set(style="whitegrid")
            fig, ax = plt.subplots(figsize=(12, 8))
            x = np.arange(0, len(reference_sequence) - window_size + 1, step_size)

            for idx, record in enumerate(alignment):
                ax.plot(x, windowed_similarities[idx], label=record.id, linewidth=2.0)

            ax.set_xlabel('Position', fontsize=14)
            ax.set_ylabel('Similarity', fontsize=14)
            ax.set_title('Sequence Similarity Plot', fontsize=16, weight='bold')
            ax.legend(title='Species', fontsize=12, title_fontsize='13', loc='upper right', bbox_to_anchor=(1.2, 1))
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.grid(True, linestyle='--', alpha=0.6)

            # Customize the look of the plot
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(1.2)
            ax.spines['bottom'].set_linewidth(1.2)

            # Save the plot to a temporary file
            plot_path = 'similarity_plot.png'
            fig.savefig(plot_path, bbox_inches='tight', dpi=300)

            # Load the plot into the QLabel
            pixmap = QPixmap(plot_path)
            self.textEditGenStats_2.setPixmap(pixmap)
            self.textEditGenStats_2.setFixedSize(900, 400)
            self.textEditGenStats_2.setScaledContents(True)

            # Optionally, remove the temporary file
            os.remove(plot_path)
            self.tabWidget.setCurrentIndex(3)

        except AttributeError as e:
            self.showErrorDialog(f"Attribute Error: {e}")
        except KeyError as e:
            self.showErrorDialog(f"Key Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")
    def load_data_climate(self):
        """
        Load climate data from a CSV file, update the UI elements with the column names, and switch to the appropriate tab.

        This method reads the climate data from a CSV file (excluding the 'id' column), updates combo boxes
        for X and Y axis settings with the column names, and switches to the climate data tab.

        Returns:
            None
        """
        try:
            # Load the data, including the 'id' column
            self.data = pd.read_csv(Params.file_name)
            self.columns = self.data.columns.tolist()
            self.columns.remove('id')  # Remove 'id' from columns for axis selection
            self.ClimaticChartSettingsAxisX.clear()
            self.ClimaticChartSettingsAxisX.addItems(self.columns)
            self.ClimaticChartSettingsAxisY.clear()
            self.ClimaticChartSettingsAxisY.addItems(self.columns)
            self.tabWidget2.setCurrentIndex(2)

        except FileNotFoundError as e:
            self.showErrorDialog(f"File Not Found Error: {e}")
        except pd.errors.EmptyDataError as e:
            self.showErrorDialog(f"Empty Data Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def generate_graph(self):
        """
        Generate and display a graph based on the selected X and Y axis data and the chosen plot type.

        This method reads the selected data columns and plot type from the UI, generates the corresponding graph,
        and displays it in the specified QLabel widget.

        Returns:
            None
        """
        try:
            x_data = self.ClimaticChartSettingsAxisX.currentText()
            y_data = self.ClimaticChartSettingsAxisY.currentText()

            # Check if the user has selected a radio button and one or both combo boxes are not properly set
            if (self.radioButtonBarGraph.isChecked() or
                    self.radioButtonScatterPlot.isChecked() or
                    self.radioButtonLinePlot.isChecked() or
                    self.radioButtonPiePlot.isChecked() or
                    self.radioButtonViolinPlot.isChecked()):

                if x_data == "Select an option" or y_data == "Select an option":
                    self.showErrorDialog("Please select valid options for both X and Y axes.")
                    return

            fig, ax = plt.subplots(figsize=(5.2, 5))  # Set figure size to 520x500 pixels (each inch is 100 pixels)

            if self.radioButtonBarGraph.isChecked():
                plot_type = 'Bar Graph'
                self.data.plot(kind='bar', x=x_data, y=y_data, ax=ax)
                # Add 'id' as labels
                for i, txt in enumerate(self.data['id']):
                    ax.text(i, self.data[y_data][i], txt, ha='center', va='bottom')
            elif self.radioButtonScatterPlot.isChecked():
                plot_type = 'Scatter Plot'
                self.data.plot(kind='scatter', x=x_data, y=y_data, ax=ax)
                # Add 'id' as labels
                for i, txt in enumerate(self.data['id']):
                    ax.annotate(txt, (self.data[x_data][i], self.data[y_data][i]))
            elif self.radioButtonLinePlot.isChecked():
                plot_type = 'Line Plot'
                self.data.plot(kind='line', x=x_data, y=y_data, ax=ax)
                # Add 'id' as labels
                for i, txt in enumerate(self.data['id']):
                    ax.text(i, self.data[y_data][i], txt, ha='center', va='bottom')
            elif self.radioButtonPiePlot.isChecked():
                plot_type = 'Pie Plot'
                self.data.set_index(x_data).plot(kind='pie', y=y_data, labels=self.data['id'], ax=ax, legend=False)
            elif self.radioButtonViolinPlot.isChecked():
                plot_type = 'Violin Plot'
                if pd.api.types.is_numeric_dtype(self.data[x_data]):
                    self.data['x_binned'] = pd.cut(self.data[x_data], bins=10)
                    self.data['x_binned'] = self.data['x_binned'].astype(str)  # Convert intervals to strings
                    sns.violinplot(x='x_binned', y=y_data, data=self.data, ax=ax)
                else:
                    sns.violinplot(x=x_data, y=y_data, data=self.data, ax=ax)

            buf = BytesIO()
            plt.savefig(buf, format='png')
            buf.seek(0)
            pixmap = QPixmap()
            pixmap.loadFromData(buf.getvalue())
            self.ClimaticChart_2.setPixmap(pixmap)
            self.tabWidget2.setCurrentIndex(2)

        except KeyError as e:
            self.showErrorDialog(f"Key Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def pressItFasta(self):
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
            fullFileName, _ = QFileDialog.getOpenFileName(None, "Select FASTA file", "../datasets",
                                                          "FASTA Files (*.fasta);;All Files (*)",
                                                          options=options)
            if fullFileName:
                update_yaml_param(Params, "params.yaml", "reference_gene_file", os.path.basename(fullFileName))
                update_yaml_param(Params, "params.yaml", "reference_gene_dir", os.path.dirname(fullFileName))

                with open(fullFileName, "r") as f:
                    self.clearGen()
                    content = f.read()
                    sequence = ""
                    for line in content.splitlines():
                        if line.startswith('>'):
                            line = f'<span style="color: green; font-weight: bold; font-size: 20px;">{line}</span>'
                            sequence += "<br>" + line + "<br>"
                        else:
                            nucleotide_colors = {
                                'A': 'yellow',
                                'C': 'blue',
                                'G': 'red',
                                'T': 'orange'
                            }
                            colored_line = ''
                            for char in line:
                                color = nucleotide_colors.get(char, '')
                                if color:
                                    colored_line += f'<span style="color: {color}; font-weight: bold; font-size: 20px;">{char}</span>'
                                else:
                                    colored_line += char
                            sequence += colored_line

                    self.textEditFasta.setHtml(
                        f"<div style='background-color: #000000; color: #ffffff; padding: 10px; white-space: pre-wrap; word-wrap: break-word;'>{sequence}</div>"
                    )
                    self.sequenceAlignmentButtonPage1.setEnabled(True)
                    self.sequenceAlignmentButtonPage1.setIcon(QIcon(":inactive/sequence.svg"))
                    self.tabWidget.setCurrentIndex(1)
        except FileNotFoundError as e:
            self.showErrorDialog(f"File Not Found Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def SeqAlign(self):
        """
        Perform sequence alignment and store the resulting genetic tree dictionary.

        This method calls the `callSeqAlign` method to perform sequence alignment and stores the resulting genetic tree dictionary in the `geneticTreeDict` attribute.
        """
        try:
            align = str(self.SequenceAlignmentMethod.currentIndex() + 1)
            fit = str(self.SequenceFitMethod.currentIndex() + 1)
            update_yaml_param(Params, "params.yaml", "alignment_method", align)
            update_yaml_param(Params, "params.yaml", "fit_method", fit)
            self.starting_position_spinbox_2.setEnabled(True)
            self.window_size_spinbox_2.setEnabled(True)
            self.geneticTreeDict = self.callSeqAlign()
        except AttributeError as e:
            self.showErrorDialog(f"Attribute Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def callSeqAlign(self):
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
            self.showErrorDialog(f"An unexpected error occurred: {error_message}")

        loading_screen = uic.loadUi("Qt/loading.ui")
        loading_screen.setWindowFlags(Qt.FramelessWindowHint)  # Remove the title bar and frame

        loading_screen.setWindowModality(Qt.ApplicationModal)

        # Set the QMovie for the movieLabel
        movie = QMovie(":active/dna.gif")  # Use the resource path for the gif
        loading_screen.movieLabel.setMovie(movie)

        # Resize the movie to fit within the QLabel
        movie.setScaledSize(QtCore.QSize(50, 50))  # Set the desired size here

        # Ensure the QLabel is centered and the GIF is properly displayed
        loading_screen.movieLabel.setAlignment(QtCore.Qt.AlignCenter)
        movie.start()

        # Show the loading screen
        loading_screen.show()
        QtWidgets.QApplication.processEvents()

        self.thread = QThread()
        self.worker = Worker(Params.reference_gene_filepath)
        self.worker.moveToThread(self.thread)

        self.worker.progress.connect(update_progress)
        self.worker.finished.connect(handle_finished)
        self.worker.error.connect(handle_error)

        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)

        self.thread.start()

        # Use a loop to wait until the thread finishes and the result is set
        while self.thread.isRunning():
            QApplication.processEvents()

        return self.geneticTrees

    def update_climate_chart(self):
        """
        Update the climate chart based on the selected condition and chart type.

        This method reads climate data from a CSV file, filters the data based on the selected condition,
        and generates a chart based on the selected chart type. The generated chart is then displayed
        in the specified QLabel widget.

        Returns:
            None
        """
        try:
            data = pd.read_csv(Params.file_name)
            condition = self.ClimStatsListCondition.currentText()
            chart_type = self.ClimStatsListChart.currentText()

            if condition == 'Temperature':
                values = data['T2M']
            elif condition == 'Wind':
                values = data['WS10M']
            elif condition == 'Humidity':
                values = data['QV2M']
            elif condition == 'Altitude':
                values = data['ALLSKY_SFC_SW_DWN']  # Assuming altitude is represented by this column
            else:
                raise ValueError(f"Unknown condition: {condition}")

            plt.figure(figsize=(9.11, 3.91))

            if chart_type == 'Bar Chart':
                values.plot(kind='bar')
            elif chart_type == 'Line Chart':
                values.plot(kind='line')
            elif chart_type == 'Pie Chart':
                values.value_counts().plot(kind='pie')
            elif chart_type == 'Area Chart':
                values.plot(kind='area')
            elif chart_type == 'Scatter Chart':
                plt.scatter(data.index, values)
            else:
                raise ValueError(f"Unknown chart type: {chart_type}")

            plt.title(f'{condition} - {chart_type}')
            plt.savefig('chart.png')
            pixmap = QPixmap('chart.png')
            self.ClimaticChart.setPixmap(pixmap)
            self.tabWidget2.setCurrentIndex(2)

        except FileNotFoundError as e:
            self.showErrorDialog(f"File Not Found Error: {e}")
        except pd.errors.EmptyDataError as e:
            self.showErrorDialog(f"Empty Data Error: {e}")
        except KeyError as e:
            self.showErrorDialog(f"Key Error: {e}")
        except ValueError as e:
            self.showErrorDialog(f"Value Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def retrieveDataNames(self, data_list):
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
            self.showErrorDialog(f"Value Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def populateMap(self, lat, long):
        """
        Create a folium map with markers based on provided latitude and longitude coordinates,
        and render it to an image.

        Args:
            lat (list): List of latitude coordinates.
            long (list): List of longitude coordinates.

        Returns:
            None
        """
        try:
            if not lat or not long or len(lat) != len(long):
                raise ValueError("Latitude and longitude lists must be non-empty and of the same length.")

            mean_lat = sum(Decimal(y) for y in lat) / len(lat)
            mean_long = sum(Decimal(x) for x in long) / len(long)

            m = folium.Map(location=[mean_lat, mean_long], zoom_start=14, tiles="OpenStreetMap")
            for i in range(len(lat)):
                folium.Marker([Decimal(lat[i]), Decimal(long[i])]).add_to(m)

            # Save the map to HTML and render it to an image
            data = io.BytesIO()
            m.save(data, close_file=False)
            html = data.getvalue().decode()

            # Use a temporary file to store the HTML content
            temp_file_path = 'temp_map.html'
            with open(temp_file_path, 'w') as f:
                f.write(html)

            # Load the HTML file in QWebEngineView and capture as an image
            self.webview = QWebEngineView()
            self.webview.setHtml(html)
            self.webview.loadFinished.connect(self.capture_image)

        except ValueError as e:
            self.showErrorDialog(f"Value Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def capture_image(self):
        """
        Capture the image from the QWebEngineView and process it for display.

        This method grabs the content of the web page displayed in the QWebEngineView
        and passes it to the `display_image` method for further processing.

        Returns:
            None
        """
        try:
            self.webview.page().grab().then(self.display_image)
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while capturing the image: {e}")

    def display_image(self, image):
        """
        Display the captured image in a QLabel or any other suitable widget.

        Args:
            image (QImage): The captured image from the QWebEngineView.

        Returns:
            None
        """
        try:
            # Convert QImage to QPixmap and display it in a QLabel (example)
            pixmap = QPixmap.fromImage(image)
            self.imageLabel.setPixmap(pixmap)
            self.imageLabel.setScaledContents(True)
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while displaying the image: {e}")

    def pressItCSV(self):
        """
        Retrieve data from a climatic file and display it in a table.

        This method allows the user to select a CSV file from the file system. It updates the relevant YAML parameters
        with the file's path, reads the file content, and displays the data in a table widget. It also processes
        location data (latitude and longitude) if available and populates a map.

        Actions:
            - Opens a file dialog to select a CSV file.
            - Updates 'file_name' parameter in the YAML file.
            - Reads and displays the content of the selected CSV file.
            - Processes and stores species and factor data.
            - Populates a map with location data if available.
            - Updates the UI to reflect the loaded data.
        """

        def is_valid_decimal(value):
            try:
                Decimal(value)
                return True
            except:
                return False

        try:
            options = QFileDialog.Options()
            options |= QFileDialog.ReadOnly
            fullFilePath, _ = QFileDialog.getOpenFileName(None, "Select CSV file", "../datasets",
                                                          "Comma Separated Values (*.csv)",
                                                          options=options)

            if fullFilePath:
                update_yaml_param(Params, "params.yaml", "file_name", fullFilePath)
                self.statisticsButtonPage2.setEnabled(True)
                with open(fullFilePath, "r") as c:
                    lines = c.readlines()
                    num_rows = len(lines)
                    first_line = lines[0].split(",")
                    lat = []
                    long = []
                    self.species = []
                    self.factors = []
                    loc = False
                    num_columns = len(first_line)
                    if first_line[-2] == 'LAT':
                        first_line_without_loc = first_line[:-2]
                        clim_data_names = self.retrieveDataNames(first_line_without_loc)
                        update_yaml_param(Params, "params.yaml", "names", first_line_without_loc)
                        loc = True
                    else:
                        clim_data_names = self.retrieveDataNames(first_line)
                        update_yaml_param(Params, "params.yaml", "names", first_line)
                    update_yaml_param(Params, "params.yaml", "data_names", clim_data_names)
                    self.textEditClimData.clear()
                    cursor = QtGui.QTextCursor(self.textEditClimData.textCursor())
                    clim_data_table = cursor.insertTable(num_rows, num_columns)
                    fmt = clim_data_table.format()
                    fmt.setWidth(QtGui.QTextLength(QtGui.QTextLength.PercentageLength, 98))
                    clim_data_table.setFormat(fmt)
                    format = QtGui.QTextCharFormat()

                    table_format = clim_data_table.format()
                    table_format.setBorder(1.5)  # Set border width
                    table_format.setBorderBrush(QtGui.QBrush(QtGui.QColor("gray")))  # Set border color
                    clim_data_table.setFormat(table_format)

                    # Create a QTextBlockFormat for center alignment
                    center_format = QtGui.QTextBlockFormat()
                    center_format.setAlignment(QtCore.Qt.AlignCenter)

                    header_format = QtGui.QTextCharFormat()
                    header_format.setFontWeight(QtGui.QFont.Bold)

                    for i, line in enumerate(lines):
                        line_split = line.split(",")
                        if line != lines[0]:
                            line_data = []
                            self.species.append(line_split[0])
                            for j in range(1, len(line_split) - (2 if loc else 0)):
                                line_data.append(line_split[j])
                            self.factors.append(line_data)

                            if loc:
                                lat.append(line_split[-2])
                                long.append(line_split[-1])

                        for j, value in enumerate(line_split):
                            if i == 0:  # If it's the first line, apply center alignment and bold format
                                cursor.setBlockFormat(center_format)
                                cursor.setCharFormat(header_format)  # Apply bold and black color format to header
                            else:
                                cursor.setCharFormat(format)  # Apply black color format to all text
                            if re.search("^[0-9\\-]*\\.[0-9]*", value) is not None:
                                cursor.insertText(str(round(Decimal(value), 3)))
                            else:
                                cursor.insertText(value)
                            cursor.movePosition(QtGui.QTextCursor.NextCell)

                    df = pd.read_csv(Params.file_name)
                    self.climaticTrees = utils.climaticPipeline(df)
                    self.tree_keys = list(self.climaticTrees.keys())
                    self.total_trees = len(self.tree_keys)
                    self.current_index = 0
                    self.climaticTreeButtonPage2.setEnabled(True)
                    self.tabWidget2.setCurrentIndex(1)

                    if loc and lat and long:
                        self.populateMap(lat, long)
                    if self.statisticsButtonPage1.isEnabled():
                        self.resultsButton.setEnabled(True)
        except FileNotFoundError as e:
            self.showErrorDialog(f"File Not Found Error: {e}")
        except pd.errors.EmptyDataError as e:
            self.showErrorDialog(f"Empty Data Error: {e}")
        except ValueError as e:
            self.showErrorDialog(f"Value Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def populateMap(self, lat, long):
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

            m = folium.Map(location=[mean_lat, mean_long],
                           zoom_start=14,
                           tiles="OpenStreetMap")
            for latitude, longitude in zip(lat, long):
                folium.Marker([latitude, longitude]).add_to(m)

            data = io.BytesIO()
            m.save(data, close_file=False)

            web_view = QWebEngineView(self.graphicsViewClimData)  # Embed the map inside graphicsViewClimData
            web_view.setHtml(data.getvalue().decode())
            layout = QtWidgets.QVBoxLayout(self.graphicsViewClimData)
            layout.addWidget(web_view)
            self.graphicsViewClimData.setLayout(layout)

        except ValueError as e:
            self.showErrorDialog(f"Value Error: {e}")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def showHomePage(self):
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
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def showGenDatPage(self):
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
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def showClimDatPage(self):
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
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def showResultsPage(self):
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
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def showResultsStatsPage(self):
        """
        Display the results statistics page of the application.

        This method sets the stacked widget's current index to 4 to display the results statistics page.
        """
        try:
            self.stackedWidget.setCurrentIndex(4)
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def showSequencePage(self):
        """
        Display the sequence alignment page.

        This method sets the stacked widget's current index to 1 and the tab widget's current index to 2
        to display the sequence alignment page.
        """
        try:
            self.stackedWidget.setCurrentIndex(1)
            self.tabWidget.setCurrentIndex(2)
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    from PyQt5.QtWidgets import QVBoxLayout, QTextBrowser

    def showFilteredResults(self):
        """
        Show the results filtered with a metric threshold provided by the user.

        This method reads the data from a CSV file, processes it through the climatic pipeline, filters the results,
        and displays the filtered results in an HTML table format within a QTextBrowser widget.
        It handles exceptions related to missing sequence alignment.

        Raises:
            AttributeError: If the sequence alignment has not been performed before attempting to generate the tree.
        """
        try:
            # df = pd.read_csv(Params.file_name)
            # utils.filterResults(self.climaticTrees, self.geneticTreeDict, df)
            df_results = pd.read_csv("./results/output.csv")

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
                    position: relative;
                    padding-left: 12px;
                    padding-right: 12px;
                }
                #styled-table td:before {
                    content: "";
                    position: absolute;
                    left: 0;
                    top: 50%;
                    transform: translateY(-50%);
                    width: 6px;
                    height: 6px;
                    border-radius: 50%;
                    background-color: #4CAF50;
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

        except AttributeError as e:
            self.textEditClimTree.setText("Please do the sequence alignment before attempting to generate the tree!")
            self.stackedWidget.setCurrentIndex(2)
            self.tabWidget2.setCurrentIndex(2)
            self.showErrorDialog(f"AttributeError: {e}")

        except FileNotFoundError as e:
            self.showErrorDialog(f"FileNotFoundError: {e}")

        except pd.errors.EmptyDataError as e:
            self.showErrorDialog(f"EmptyDataError: {e}")

        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def toggleDarkMode(self):
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
            buttons_Vertical = [self.fileBrowserButtonPage1, self.sequenceAlignmentButtonPage1,
                                self.clearButtonPage1, self.statisticsButtonPage1, self.geneticTreeButtonPage1,
                                self.fileBrowserButtonPage2, self.clearButtonPage2,
                                self.climaticTreeButtonPage2, self.statisticsButtonPage2,
                                self.settingsButtonPage3, self.settingsButtonPage4,
                                self.submitButtonPage3, self.statisticsButtonPage3,
                                self.submitButtonPage4, self.statisticsButtonPage4,
                                self.clearButtonPage3, self.StartSequenceAlignmentButton,
                                self.clearButtonPage4]
            buttons = [self.geneticDataButton, self.climaticDataButton, self.helpButton, self.homeButton,
                       self.resultsButton]

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
            self.darkModeButton.setStyleSheet(self.getButtonStyle(self.isDarkMode, True))
            self.darkModeButton.setGraphicsEffect(self.createShadowEffect(10, 140))

            for button in buttons:
                button.setCursor(Qt.PointingHandCursor)
                button.setStyleSheet(self.getButtonStyle(self.isDarkMode))
                button.setGraphicsEffect(self.createShadowEffect(10, 140))

            for buttonV in buttons_Vertical:
                buttonV.setCursor(Qt.PointingHandCursor)
                buttonV.setStyleSheet(self.getButtonStyle(self.isDarkMode, False, True))
                buttonV.setGraphicsEffect(self.createShadowEffect(10, 110))

        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def getButtonStyle(self, darkMode, isMainButton=False, isVertical=False):
        """
        Generate the appropriate style for buttons based on the mode (dark/light) and type (main/vertical).

        Args:
            darkMode (bool): A flag indicating if dark mode is enabled.
            isMainButton (bool): A flag indicating if the button is a main button.
            isVertical (bool): A flag indicating if the button is a vertical button.

        Returns:
            str: The stylesheet string for the button.
        """
        if darkMode:
            background_color = "#646464" if isMainButton else "#464645"
            hover_color = "#B7B7B6" if isMainButton else "#9F9F9F"
        else:
            background_color = "#DEDDDA" if isMainButton else "#EEEEEE"
            hover_color = "#B7B7B6" if isMainButton else "#D7D7D7"

        return f"""
            QPushButton {{
                padding: 10px 20px;
                font-weight: bold;
                background-color: {background_color};
                border-radius: 20px;
                transition: background-color 0.3s ease; /* Add transition */
            }}
            QPushButton:hover {{
                background-color: {hover_color}; 
            }}
            QPushButton:pressed {{
                background-color: {background_color};
            }}
        """

    def createShadowEffect(self, blur_radius, alpha):
        """
        Create a shadow effect for the buttons.

        Args:
            blur_radius (int): The blur radius for the shadow.
            alpha (int): The alpha (transparency) value for the shadow color.

        Returns:
            QGraphicsDropShadowEffect: The configured shadow effect.
        """
        shadow_effect = QGraphicsDropShadowEffect()
        shadow_effect.setBlurRadius(blur_radius)
        shadow_effect.setColor(QColor(0, 0, 0, alpha))
        shadow_effect.setOffset(3, 3)
        return shadow_effect

    # press the button to delete data
    def clearGen(self):
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
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def clearClim(self):
        """
        Clear the text fields related to climatic data.

        This method disables buttons related to climatic data and clears the necessary fields.
        """
        try:
            self.statisticsButtonPage2.setEnabled(False)
            self.climaticTreeButtonPage2.setEnabled(False)
            self.resultsButton.setEnabled(False)
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}", "Error")

    def clearResult(self):
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
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}", "Error")

    def clearResultStat(self):
        """
        Clear the statistics result lists.

        This method resets the current index of the results statistics condition list and the results statistics chart list to 0.
        """
        try:
            self.ResultsStatsListCondition.setCurrentIndex(0)
            self.ResultsStatsListChart.setCurrentIndex(0)
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}", "Error")



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
        try:
            self.tabWidget.setCurrentIndex(4)
            file_path = "results/geneticTrees.json"
            with open(file_path, 'r') as file:
                self.newick_json = json.load(file)

            self.tree_keys = list(self.newick_json.keys())
            self.total_trees = len(self.tree_keys)
            self.current_index = 0
            self.geneticTreescomboBox.clear()

            # Format the tree keys to replace underscore with ' nt '
            formatted_tree_keys = [self.format_tree_name(key) for key in self.tree_keys]
            self.geneticTreescomboBox.addItems(formatted_tree_keys)


            self.show_tree(self.current_index)

        except FileNotFoundError as e:
            self.showErrorDialog(f"The file {file_path} was not found: {e}", "File Not Found")
        except json.JSONDecodeError as e:
            self.showErrorDialog(f"An error occurred while decoding the JSON file: {e}", "JSON Decode Error")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def format_tree_name(self, tree_name):
        """
        Format the tree name by replacing underscores with ' nt '.

        Args:
            tree_name (str): The original tree name.

        Returns:
            str: The formatted tree name.
        """
        try:
            parts = tree_name.split('_')
            if len(parts) == 2:
                return f"{parts[0]} nt {parts[1]} nt"
            return tree_name
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while formatting the tree name: {e}")

    def show_selected_tree(self, index):
        """
        Display the selected tree based on the provided index.

        Args:
            index (int): The index of the selected tree in the combo box.

        Returns:
            None
        """
        try:
            if index >= 0:
                self.show_tree(index)
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while displaying the selected tree: {e}")

    def show_tree(self, index):
        """
        Display the phylogenetic tree at the specified index using Toytree.

        Args:
            index (int): The index of the tree to display.

        Returns:
            None
        """
        try:
            if 0 <= index < self.total_trees:
                self.current_index = index  # Keep track of the current index
                key = self.tree_keys[index]
                newick_str = self.newick_json[key]

                # Read the tree using Toytree
                tree = toytree.tree(newick_str)

                # Customize the tip labels and their style
                tip_labels = tree.get_tip_labels()
                custom_tip_labels = [
                    '{}. {}'.format(i[0], i[1:]) for i in tip_labels
                ]

                # Draw the tree with customized style
                canvas, axes, mark = tree.draw(
                    width=921,
                    height=450,
                    tip_labels=custom_tip_labels,
                    tip_labels_style={"font-size": "15px"},
                    fixed_order=tip_labels,
                    edge_type='c'  # This line sets the edge type to 'c' as specified in the image
                )

                # Save the canvas to a temporary file
                temp_img_file = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
                temp_img_path = temp_img_file.name
                toyplot.png.render(canvas, temp_img_path)

                # Create a QPixmap from the temporary image file
                pixmap = QPixmap(temp_img_path)

                # Clear the QLabel before setting the new QPixmap
                self.GeneticTreeLabel.clear()
                self.GeneticTreeLabel.setPixmap(pixmap)
                self.GeneticTreeLabel.adjustSize()
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while rendering the tree: {e}")

    def download_graph(self):
        """
        Download the current displayed graph as a PNG file.

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
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Graph As", default_file_name,
                                                       "PNG Files (*.png);;All Files (*)", options=options)
            if file_path:
                if not file_path.lower().endswith('.png'):
                    file_path += '.png'
                with open(self.temp_img_path, 'rb') as temp_file:
                    with open(file_path, 'wb') as file:
                        file.write(temp_file.read())
        except FileNotFoundError as e:
            self.showErrorDialog(f"The temporary image file was not found: {e}", "File Not Found")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while downloading the graph: {e}")

    def open_preferences_dialog(self):
        """
        Open the preferences dialog and update the application settings based on user input.

        This method opens a PreferencesDialog, updates it with the current preferences, and applies the new preferences if the user accepts the changes.

        Returns:
            None
        """
        try:
            dialog = PreferencesDialog(self)
            dialog.update_preferences(self.preferences)
            if dialog.exec_() == QDialog.Accepted:
                self.preferences = dialog.get_preferences()
                self.apply_preferences()
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while opening the preferences dialog: {e}")

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
            self.showErrorDialog(f"An unexpected error occurred while applying preferences: {e}")

    def displayClimaticTrees(self):
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
            self.showErrorDialog(f"An error occurred while accessing the climatic trees: {e}", "Key Error")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

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
            self.showErrorDialog(
                f"An unexpected error occurred while displaying the selected climatic tree: {e}")

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
                    self.render_network_view(tree, label_color, edge_color, reticulation_color, layout,
                                             proportional_edge_lengths, label_internal_vertices, use_leaf_names,
                                             show_branch_length)
                else:
                    self.render_tree_view(tree, label_color, edge_color, reticulation_color, layout,
                                          proportional_edge_lengths, label_internal_vertices, use_leaf_names,
                                          show_branch_length)
        except KeyError as e:
            self.showErrorDialog(f"An error occurred while accessing the climatic tree data: {e}", "Key Error")
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred: {e}")

    def render_network_view(self, tree, label_color, edge_color, reticulation_color, layout, proportional_edge_lengths,
                            label_internal_vertices, use_leaf_names, show_branch_length):
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
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)',
                width=911,
                height=441
            )

            # Add edge annotations for branch lengths
            for annotation in edge_annotations:
                fig.add_annotation(annotation)

            temp_img_file = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
            self.temp_img_path = temp_img_file.name
            pio.write_image(fig, self.temp_img_path, format="png")
            temp_img_file.close()

            pixmap = QPixmap(self.temp_img_path)
            self.climaticTreesLabel.clear()
            self.climaticTreesLabel.setPixmap(pixmap)
            self.climaticTreesLabel.adjustSize()
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while rendering the network view: {e}")

    def render_tree_view(self, tree, label_color, edge_color, reticulation_color, layout, proportional_edge_lengths,
                         label_internal_vertices, use_leaf_names, show_branch_length):
        """
        Render the tree view of the phylogenetic tree.

        This method applies user preferences to the tree and renders it using Matplotlib.

        Args:
            tree (Phylo.BaseTree.Tree): The phylogenetic tree to render.
            label_color (str): Color for labels.
            edge_color (str): Color for edges.
            reticulation_color (str): Color for reticulation edges.
            layout (str): Layout for the tree visualization (e.g., 'horizontal').
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
                return f'{label}\n{clade.branch_length:.2f}' if show_branch_length and label else label

            Phylo.draw(tree, do_show=False, axes=ax, label_func=label_func,
                       label_colors={clade: label_color for clade in tree.find_clades()})

            ax.axis('off')  # Remove axes

            temp_img_file = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
            self.temp_img_path = temp_img_file.name
            plt.savefig(self.temp_img_path, format="png")
            temp_img_file.close()
            plt.close(fig)

            pixmap = QPixmap(self.temp_img_path)
            self.climaticTreesLabel.clear()
            self.climaticTreesLabel.setPixmap(pixmap)
            self.climaticTreesLabel.adjustSize()
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while rendering the tree view: {e}")

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
                mode='markers+text' if use_leaf_names else 'markers',
                textposition="top center",
                hoverinfo='text',
                marker=dict(
                    showscale=False,  # Disable color scale
                    colorscale='Viridis',  # Use a valid Plotly colorscale
                    size=10,
                    line_width=2,
                    color=label_color))

            for node in graph.nodes():
                x, y = pos[node]
                node_trace['x'] += (x,)
                node_trace['y'] += (y,)
                if use_leaf_names:
                    name = node.name if hasattr(node, 'name') and node.name else ''
                    node_trace['text'] += (name,)

            return node_trace
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while creating the node trace: {e}")

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
            self.showErrorDialog(f"An unexpected error occurred while getting the layout: {e}")

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
                hoverinfo='none',
                mode='lines'
            )
            edge_annotations = []

            for clade in tree.find_clades(order="level"):
                if clade.is_terminal():
                    continue
                for child in clade.clades:
                    x0, y0 = pos[clade]
                    x1, y1 = pos[child]
                    edge_trace['x'] += (x0, x1, None)
                    edge_trace['y'] += (y0, y1, None)

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
                                font=dict(color=edge_color)
                            )
                        )
            return edge_trace, edge_annotations
        except Exception as e:
            self.showErrorDialog(f"An unexpected error occurred while creating the edge trace: {e}")

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
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Graph As", default_file_name,
                                                       "PNG Files (*.png);;All Files (*)", options=options)
            if file_path:
                if not file_path.lower().endswith('.png'):
                    file_path += '.png'
                with open(self.temp_img_path, 'rb') as temp_file:
                    with open(file_path, 'wb') as file:
                        file.write(temp_file.read())
        except FileNotFoundError as e:
            self.showErrorDialog(f"The temporary image file was not found: {e}", "File Not Found")
        except Exception as e:
            self.showErrorDialog(
                f"An unexpected error occurred while downloading the climatic tree graph: {e}")


if __name__ == "__main__":
    try:
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
        sys.exit(app.exec_())
    except Exception as e:
        QtWidgets.QMessageBox.critical(None, "Application Error", f"An unexpected error occurred: {e}")
        sys.exit(1)
