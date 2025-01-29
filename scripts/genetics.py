import json
import os
import shutil
from collections import Counter, defaultdict

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import toyplot.png
import toytree
from aphylogeo.params import Params
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.ticker import MaxNLocator
from PyQt6 import QtCore, QtWidgets, uic
from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QIcon, QMovie, QPixmap
from PyQt6.QtWidgets import QApplication, QFileDialog
from Qt import loading
from utils.error_dialog import show_error_dialog
from utils.genetic_params_dialog import ParamDialog
from utils.my_dumper import update_yaml_param
from Worker import Worker


class Genetics:
    def __init__(self, main):
        self.main = main
        self.worker = None
        self.msa = []
        self.geneticTrees = []

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
            for _, value in msa_data.items():
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
            show_error_dialog(f"Key Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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
            show_error_dialog(f"Value Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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

            _, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 4), gridspec_kw={"height_ratios": [1, 8]})
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
            show_error_dialog(f"Key Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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
            show_error_dialog(f"Key Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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
            show_error_dialog(f"An unexpected error occurred: {e}")

    def update_plot(self):
        """
        Update the plot based on the current starting position and window size.

        This method reads the MSA data, standardizes the sequence lengths, and plots the alignment chart.
        The resulting plot is displayed in the specified label widget and the tab widget is updated.

        Returns:
            None
        """
        try:
            starting_position = self.main.starting_position_spinbox_2.value()
            window_size = self.main.window_size_spinbox_2.value()
            output_path = "scripts/results/sequence_alignment_plot.png"

            genetic_data = self.read_msa(self.msa)
            standardized_data = self.standardize_sequence_lengths(genetic_data)
            self.plot_alignment_chart(standardized_data, starting_position, window_size, output_path)

            pixmap = QPixmap(output_path)
            self.main.seqAlignLabel.setPixmap(pixmap)
            self.main.tabWidget.setCurrentIndex(2)

        except AttributeError as e:
            show_error_dialog(f"Attribute Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def initialize_species_list(self):
        # Load species names into the combo box
        self.main.referenceComboBox.clear()
        unique_species = set()
        for _, value in self.msa.items():
            parts = value.strip().split("\n")
            for i in range(0, len(parts), 2):
                header = parts[i].strip(">").replace("_", " ")  # Replace underscores with spaces
                unique_species.add(header)

        # Add unique species to the combo box
        for species in unique_species:
            self.main.referenceComboBox.addItem(species)

        # Optionally set the first species as the default selected item
        if self.main.referenceComboBox.count() > 0:
            self.main.referenceComboBox.setCurrentIndex(0)

        # Connect the value change signals to update the plot
        self.main.similarityWindowSizeSpinBox.valueChanged.connect(self.update_similarity_plot)
        self.main.startingPositionSimilaritySpinBox.valueChanged.connect(self.update_similarity_plot)
        self.main.referenceComboBox.currentIndexChanged.connect(self.update_similarity_plot)

        # Generate the initial plot
        self.update_similarity_plot()

    def update_similarity_plot(self):
        try:
            window_size = self.main.similarityWindowSizeSpinBox.value()
            start_pos = self.main.startingPositionSimilaritySpinBox.value()
            reference_species = self.main.referenceComboBox.currentText().replace(" ", "_")  # Convert back to original format

            sequences = defaultdict(str)
            for _, value in self.msa.items():
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
            self.main.textEditGenStats_2.setPixmap(pixmap)
            self.main.textEditGenStats_2.setFixedSize(900, 400)
            self.main.textEditGenStats_2.setScaledContents(True)

            self.main.tabWidget.setCurrentIndex(3)
        except Exception as e:
            print(f"Error updating similarity plot: {e}")

    def download_similarity_plot_chart(self):
        file_url = "scripts/results/similarity_plot.png"  # The file path
        save_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", "PNG Files (*.png);;All Files (*)")
        if not save_path:
            return  # User cancelled th
        shutil.copy(file_url, save_path)  # e save dialog

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
            options = QFileDialog.Option.ReadOnly
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

                    self.main.textEditFasta.setHtml(
                        f"<div style='background-color: #ffffff; color: #000000; padding: 10px; white-space: pre-wrap; word-wrap: break-word;'>{sequence}</div>"
                    )
                    self.main.sequenceAlignmentButtonPage1.setEnabled(True)
                    self.main.sequenceAlignmentButtonPage1.setIcon(QIcon(":inactive/sequence.svg"))
                    self.main.tabWidget.setCurrentIndex(1)
        except FileNotFoundError as e:
            show_error_dialog(f"File Not Found Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def start_alignment_analysis(self):
        """
        Perform sequence alignment and store the resulting genetic tree dictionary.

        This method calls the `callSeqAlign` method to perform sequence alignment and stores the resulting genetic tree dictionary in the `geneticTreeDict` attribute.
        """
        self.main.starting_position_spinbox_2.setEnabled(True)
        self.main.window_size_spinbox_2.setEnabled(True)
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
                item.setCheckState(Qt.CheckState.Checked)
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
            self.main.geneticTreeButtonPage1.setEnabled(True)
            self.update_plot()
            self.main.statisticsButtonPage1.setEnabled(True)
            if self.main.climaticTreeButtonPage2.isEnabled():
                self.main.resultsButton.setEnabled(True)

        def handle_error(error_message):
            loading_screen.close()
            show_error_dialog(f"An unexpected error occurred: {error_message}")

        if loading_screen := loading.Ui_LoadingDialog():
            loading_screen.setupUi(loading_screen)
            # loading_screen.close()
            # loading_screen = uic.loadUi("scripts/Qt/loading.ui")
            loading_screen.setWindowFlags(Qt.WindowType.FramelessWindowHint)  # Remove the title bar and frame

            loading_screen.setWindowModality(Qt.ApplicationModal)

            # Set the QMovie for the movieLabel
            movie = QMovie(":active/dna.gif")  # Use the resource path for the gif
            loading_screen.movieLabel.setMovie(movie)

            # Resize the movie to fit within the QLabel
            movie.setScaledSize(QtCore.QSize(100, 100))  # Set the desired size here

            # Ensure the QLabel is centered and the GIF is properly displayed
            loading_screen.movieLabel.setAlignment(QtCore.Qt.AlignCenter)
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

    def show_sequence_alignment_page(self):
        """
        Display the sequence alignment page.

        This method sets the stacked widget's current index to 1 and the tab widget's current index to 2
        to display the sequence alignment page.
        """
        try:
            self.main.stackedWidget.setCurrentIndex(1)
            self.main.tabWidget.setCurrentIndex(2)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def clear_genetic_data(self):
        """
        Clear the genetic data fields.

        This method clears the content of the text edit widgets related to FASTA sequences, sequence alignments, and genetic trees.
        It also resets the current index of the genetic statistics list to 0 and disables relevant buttons.
        """
        try:
            self.main.textEditFasta.clear()
            self.main.seqAlignLabel.clear()
            self.main.textEditGenStats_2.clear()
            self.main.sequenceAlignmentButtonPage1.setEnabled(False)
            self.main.statisticsButtonPage1.setEnabled(False)
            self.main.geneticTreeButtonPage1.setEnabled(False)
            self.main.GeneticTreeLabel.clear()
            self.main.resultsButton.setEnabled(False)
            self.geneticTreeDict = None
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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

        self.main.tabWidget.setCurrentIndex(4)
        file_path = "scripts/results/geneticTrees.json"
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
            show_error_dialog(f"The tree image file was not found: {e}", "File Not Found")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred while downloading the tree image: {e}")

    def stopWorker(self):
        if self.worker:
            self.worker.stop()

    def open_genetic_settings_window(self):
        dialog = ParamDialog()
        dialog.exec()
