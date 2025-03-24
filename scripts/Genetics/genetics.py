import os


from aphylogeo.params import Params
from Genetics.genetics_tree import GeneticTree
from Genetics.genetics_alignment import GeneticAlignment
from Genetics.genetics_statistics import GeneticStatistics
from PyQt6 import QtWidgets
from utils.error_dialog import show_error_dialog
from utils.my_dumper import update_yaml_param


class Genetics:
    def __init__(self, main):
        self.main = main
        self.geneticTree = GeneticTree(main)
        self.geneticAlignment = GeneticAlignment(main)
        self.geneticStatistics = GeneticStatistics(main, self.geneticAlignment)
        
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
            self.geneticAlignment.geneticTrees = None

            self.main.sequenceAlignmentButtonPage1.setEnabled(False)
            self.main.statisticsButtonPage1.setEnabled(False)
            self.main.geneticTreeButtonPage1.setEnabled(False)
            self.main.resultsButton.setEnabled(False)
            self.main.tabWidget.setCurrentIndex(0)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

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
