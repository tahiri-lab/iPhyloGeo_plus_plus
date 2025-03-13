import os
from typing import Any, Dict

import pandas as pd
import plotly.express as px
from PyQt6.QtWidgets import QFileDialog
import plotly.io as pio

import numpy as np
from aphylogeo.params import Params
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from event_connector import blocked_signals
from Genetics.genetic_params_dialog import ParamDialog
from Genetics.genetics_plot_chart import plot_alignment_chart, read_msa, standardize_sequence_lengths
from Genetics.genetics_tree import GeneticTree
from PyQt6 import QtCore, QtWidgets
from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QMovie
from ui import loading_dialog
from utils.download_file import download_local_with_fig
from utils.error_dialog import show_error_dialog
from utils.file_caching import FileCaching
from utils.my_dumper import update_yaml_param
from worker import Worker


class Genetics:
    def __init__(self, main):
        self.main = main
        self.geneticTree = GeneticTree(main)
        self.worker = None
        self.msa = {}
        self.geneticTrees = []
        self.customConfig ={
                'modeBarButtonsToRemove': ['toImage'],
                'displaylogo' : False
            }

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

            standardized_data = standardize_sequence_lengths(self.msa)
            pixmap = plot_alignment_chart(standardized_data, starting_position, window_size, self.main.isDarkMode)

            self.main.seqAlignLabel.setPixmap(pixmap)
            self.main.tabWidget.setCurrentIndex(2)

        except AttributeError as e:
            show_error_dialog(f"Attribute Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def initialize_species_list(self):
        # Load species names into the combo box
        unique_species = set(self.msa.keys())

        with blocked_signals(self.main.referenceComboBox):
            self.main.referenceComboBox.clear()
            self.main.referenceComboBox.addItems(unique_species)
            self.main.referenceComboBox.setCurrentIndex(0)

        # Generate the initial plot
        self.update_similarity_plot()

    def update_similarity_plot(self):
        try:
            reference_species = self.main.referenceComboBox.currentText().replace(" ", "_")  # Convert back to original format

            sequences = self.msa

            max_len = max(len(seq) for seq in sequences.values())
            padded_records = []
            for header, sequence in sequences.items():
                padded_seq = sequence.ljust(max_len, "-")
                padded_records.append(SeqRecord(Seq(padded_seq), id=header))

            alignment = MultipleSeqAlignment(padded_records)
            reference_index = [record.id for record in alignment].index(reference_species.replace("_", " "))
            reference_sequence = str(alignment[reference_index].seq)
            del alignment[reference_index]
            
            similarities = []
            for record in alignment:            
                similarity = [1 if ref == res else 0 for ref, res in zip(reference_sequence, str(record.seq))]
                similarities.append(similarity)
            similarities = np.array(similarities)
            
            
            def sliding_window_avg(arr, window_size, step_size):
                return [np.mean(arr[i : i + window_size]) for i in range(0, len(arr) - window_size + 1, step_size)]

            step_size = 10  # Adjusted step size for better plotting
            windowSize = 100

            windowed_similarities = []
            for sim in similarities:
                windowed_similarities.append(sliding_window_avg(sim, windowSize, step_size))

            windowed_similarities = np.array(windowed_similarities)

            x = np.arange(0, len(reference_sequence) - windowSize + 1, step_size)

            data = {
                "Position": [],
                "Similarity": [],
                "Sequence": []
            }
            for idx, record in enumerate(alignment):
                for pos, sim_val in zip(x, windowed_similarities[idx]):
                    data["Position"].append(pos)
                    data["Similarity"].append(sim_val)
                    data["Sequence"].append(record.id)
            df = pd.DataFrame(data)

            self.fig = px.line(
                df,
                x="Position",
                y="Similarity",
                color="Sequence",
                custom_data=["Sequence"],
                labels={
                    "Position": "Position",
                    "Similarity": "Similarity",
                    "Sequence": "Sequence ID"
                },
            )
            
            self.similarity_plot_style()
            
            fig_html = self.fig.to_html(include_plotlyjs="cdn", config = self.customConfig)

            self.main.textEditGenStats_2.setHtml(fig_html)

            self.main.tabWidget.setCurrentIndex(3)
        except Exception as e:
            show_error_dialog(f"Error updating similarity plot: {e}")
            
    def similarity_plot_style(self):
        if self.main.isDarkMode:
            self.fig.update_layout(template="plotly_dark")
        else:
            self.fig.update_layout(template="plotly_white")

        self.fig.update_layout(
            hovermode="x",
            title_text="Sequence Similarity Plot",
            title_x=0.5,                            
            title_y=0.95,                          
            title_pad=dict(t=0, b=10)               
        )
        
        self.fig.update_traces(
            hovertemplate="<b>Sequence:</b> %{customdata[0]}<br>"
                            "<b>Similarity:</b> %{y:.2f}<extra></extra>")

    def download_similarity_plot_chart(self):
        download_local_with_fig(self.fig)
            
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
            options = QtWidgets.QFileDialog.Option.ReadOnly
            fullFileName, _ = QtWidgets.QFileDialog.getOpenFileName(
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
                    self.main.textEditFasta.setHtml(
                        f"<div style='padding: 10px; white-space: pre-wrap; word-wrap: break-word;'>{generate_html_from_fasta(f, self.main.isDarkMode)}</div>"
                    )
                    self.main.sequenceAlignmentButtonPage1.setEnabled(True)
                    self.main.tabWidget.setCurrentIndex(1)
        except FileNotFoundError as e:
            show_error_dialog(f"File Not Found Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def start_alignment_analysis(self):
        """
        Perform sequence alignment and store the resulting genetic tree dictionary.

        This method calls the `callSeqAlign` method to perform sequence alignment and stores the resulting genetic tree dictionary in the `geneticTrees` attribute.
        """
        self.main.starting_position_spinbox_2.setEnabled(True)
        self.main.window_size_spinbox_2.setEnabled(True)
        self.call_seq_align()

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
                QtWidgets.QApplication.processEvents()
            else:
                loading_screen.progressBar.setValue(100)
                QtWidgets.QApplication.processEvents()

        def handle_finished(result: Dict[str, Any]):
            loading_screen.close()
            msa = result["msa"]
            self.msa = read_msa(msa)
            self.geneticTrees = result["geneticTrees"]
            self.update_plot()
            self.main.geneticTreeButtonPage1.setEnabled(True)
            self.main.statisticsButtonPage1.setEnabled(True)
            if self.main.climaticTreeButtonPage2.isEnabled():
                self.main.resultsButton.setEnabled(True)

        def handle_error(error_message):
            loading_screen.close()
            show_error_dialog(f"An unexpected error occurred: {error_message}")

        if loading_screen := loading_dialog.LoadingDialog():
            loading_screen.setWindowFlags(Qt.WindowType.FramelessWindowHint)  # Remove the title bar and frame

            loading_screen.setWindowModality(Qt.WindowModality.ApplicationModal)

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

        if (result := FileCaching.get_cached_result_file(Params.reference_gene_filepath)) is not None:
            print("Using cached result")
            handle_finished(result)
            return

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
            QtWidgets.QApplication.processEvents()

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
            self.main.GeneticTreeLabel.clear()
            self.geneticTrees = None

            self.main.sequenceAlignmentButtonPage1.setEnabled(False)
            self.main.statisticsButtonPage1.setEnabled(False)
            self.main.geneticTreeButtonPage1.setEnabled(False)
            self.main.resultsButton.setEnabled(False)
            self.main.tabWidget.setCurrentIndex(0)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def stopWorker(self):
        if self.worker:
            self.worker.stop()


def open_genetic_settings_window():
    dialog = ParamDialog()
    dialog.exec()


def generate_html_from_fasta(file, isDarkMode):
    content = file.read()
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
                "T": "gray" if isDarkMode else "black",
            }
            colored_line = ""
            for char in line:
                color = nucleotide_colors.get(char, "")
                if color:
                    colored_line += f'<span style="color: {color}; font-size: 20px;">{char}</span>'
                else:
                    colored_line += char
            sequence += colored_line + " "
    return sequence
