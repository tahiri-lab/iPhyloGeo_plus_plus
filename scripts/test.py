import sys
import os
from collections import Counter
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QWidget, QPushButton, QSpinBox
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


class Params:
    reference_gene_filepath = "F:/STUDY/projects/PFE/aPhyloGeo_plus_plus/datasets/small_seq.fasta"
    window_size = 50
    starting_position = 1


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Genetic Alignment Chart")
        self.setGeometry(100, 100, 1200, 800)

        self.textEditGenStats_2 = QLabel(self)
        self.textEditGenStats_2.setGeometry(50, 50, 1100, 600)
        self.textEditGenStats_2.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.textEditGenStats_2.setWordWrap(True)
        self.textEditGenStats_2.setFixedSize(1100, 600)

        self.starting_position_spinbox = QSpinBox(self)
        self.starting_position_spinbox.setGeometry(50, 670, 100, 30)
        self.starting_position_spinbox.setRange(1, 1000)
        self.starting_position_spinbox.setValue(Params.starting_position)

        self.window_size_spinbox = QSpinBox(self)
        self.window_size_spinbox.setGeometry(160, 670, 100, 30)
        self.window_size_spinbox.setRange(1, 1000)
        self.window_size_spinbox.setValue(Params.window_size)

        self.button = QPushButton('Load Data', self)
        self.button.clicked.connect(self.load_data_genetic)
        self.button.setGeometry(270, 670, 100, 30)

        layout = QVBoxLayout()
        layout.addWidget(self.textEditGenStats_2)
        layout.addWidget(self.starting_position_spinbox)
        layout.addWidget(self.window_size_spinbox)
        layout.addWidget(self.button)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def load_data_genetic(self):
        print(f"Loading data from: {Params.reference_gene_filepath}")
        genetic_data = self.read_fasta(Params.reference_gene_filepath)
        standardized_data = self.standardize_sequence_lengths(genetic_data)
        starting_position = self.starting_position_spinbox.value()
        window_size = self.window_size_spinbox.value()
        output_path = "./alignment_chart.png"
        self.plot_alignment_chart(standardized_data, starting_position, window_size, output_path)
        self.display_image(output_path)

    def read_fasta(self, fasta_file):
        genetic_data = {}
        with open(fasta_file, 'r') as file:
            sequence_id = None
            sequence_data = []
            for line in file:
                if line.startswith('>'):
                    if sequence_id is not None:
                        genetic_data[sequence_id] = ''.join(sequence_data)
                    sequence_id = line[1:].strip()
                    sequence_data = []
                else:
                    sequence_data.append(line.strip())
            if sequence_id is not None:
                genetic_data[sequence_id] = ''.join(sequence_data)
        return genetic_data

    def standardize_sequence_lengths(self, genetic_data):
        max_length = max(len(seq) for seq in genetic_data.values())
        standardized_data = {}
        for key, seq in genetic_data.items():
            if len(seq) < max_length:
                padded_seq = seq + '-' * (max_length - len(seq))
            else:
                padded_seq = seq[:max_length]
            standardized_data[key] = padded_seq
        return standardized_data

    def plot_alignment_chart(self, genetic_data, starting_position, window_size, output_path):
        end_position = starting_position + window_size
        truncated_data = {key: value[starting_position:end_position] for key, value in genetic_data.items()}

        alignment = MultipleSeqAlignment([
            SeqRecord(Seq(seq), id=key)
            for key, seq in truncated_data.items()
        ])

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 8, 1]})
        ax1.set_axis_off()
        ax2.set_axis_off()
        ax3.set_axis_off()

        # Calculate conservation and gaps
        conservation = self.calculate_conservation(alignment)
        gaps = self.calculate_gaps(alignment)

        # Plot conservation
        bar_width = 1.0
        ax1.bar(range(len(conservation)), conservation, color='#4CAF50', width=bar_width, align='edge')
        ax1.set_xlim(0, len(conservation))
        ax1.set_ylim(0, 1)

        # Plot gaps
        ax3.bar(range(len(gaps)), gaps, color='#F44336', width=bar_width, align='edge')
        ax3.set_ylabel('Gap')
        ax3.set_xlim(0, len(gaps))
        ax3.set_ylim(0, 1)

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
                ax2.text(j + rect_width / 2, len(seqs) - i - 1 + rect_height / 2, nucleotide, ha='center', va='center',
                         fontsize=font_size)

        consensus_seq = self.calculate_consensus(alignment)
        for j, nucleotide in enumerate(consensus_seq):
            color = colors.get(nucleotide, 'black')  # Default to black if not found
            rect = mpatches.Rectangle((j, -1), rect_width, rect_height, color=color)
            ax2.add_patch(rect)
            ax2.text(j + rect_width / 2, -1 + rect_height / 2, nucleotide, ha='center', va='center', fontsize=font_size,
                     fontweight='bold')

        # Add Consensus title
        ax2.text(-1, -1 + rect_height / 2, "Consensus", ha='right', va='center', fontsize=font_size, fontweight='bold')

        ax2.set_xlim(0, len(consensus_seq))
        ax2.set_ylim(-2, len(seqs))
        ax2.set_yticks(range(len(ids)))
        ax2.set_yticklabels(ids, fontsize=font_size)
        ax2.set_xticks(range(len(consensus_seq)))
        ax2.set_xticklabels(range(starting_position, end_position), fontsize=font_size)

        plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.05)
        plt.savefig(output_path)
        plt.close()
        print(f"Plot saved to: {output_path}")

    def calculate_conservation(self, alignment):
        conservation = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            counts = Counter(column)
            most_common = counts.most_common(1)[0][1]
            conservation.append(most_common / len(column))
        return conservation

    def calculate_gaps(self, alignment):
        gaps = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            gap_count = column.count('-')
            gaps.append(gap_count / len(column))
        return gaps

    def calculate_consensus(self, alignment):
        consensus_seq = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            most_common = Counter(column).most_common(1)[0][0]
            consensus_seq.append(most_common)
        return ''.join(consensus_seq)

    def display_image(self, image_path):
        if os.path.exists(image_path):
            print(f"Displaying image: {image_path}")
            pixmap = QPixmap(image_path)
            if pixmap.isNull():
                print("Failed to load image into QPixmap")
            else:
                print("Image loaded into QPixmap successfully")
            self.textEditGenStats_2.setPixmap(pixmap)
            self.textEditGenStats_2.repaint()  # Ensure the label is updated
        else:
            print(f"Image not found: {image_path}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWin = MainWindow()
    mainWin.show()
    sys.exit(app.exec_())
