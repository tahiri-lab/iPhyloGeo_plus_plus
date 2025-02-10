from utils.error_dialog import show_error_dialog
from utils.download_file import download_file_temporary_PLT

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter

from collections import defaultdict

def read_msa(msa_data):
    """
    Reads multiple sequence alignment (MSA) data and organizes it into a dictionary.

    Args:
        msa_data (dict): A dictionary containing MSA data where keys are identifiers and values are sequences.

    Returns:
        dict: A dictionary with sequence identifiers as keys and concatenated sequences as values.
    """
    try:       
        msa = defaultdict(str)
        for _, value in msa_data.items():
                parts = value.strip().split("\n")
                for i in range(0, len(parts), 2):
                    header = parts[i].strip(">")
                    sequence = parts[i + 1]
                    msa[header.replace("_", " ")] += sequence  # Replace underscores with spaces
                                   
        return msa

    except KeyError as e:
        show_error_dialog(f"Key Error: {e}")
    except Exception as e:
        show_error_dialog(f"An unexpected error occurred: {e}")

def standardize_sequence_lengths(genetic_data):
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

def plot_alignment_chart(genetic_data, starting_position, window_size):
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
        alignment = MultipleSeqAlignment([SeqRecord(Seq(seq), id=key) for key, seq in truncated_data.items()])

        _, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 4), gridspec_kw={"height_ratios": [1, 8]})
        ax1.set_axis_off()
        ax2.set_axis_off()

        # Calculate conservation
        conservation, _ = calculate_conservation_and_gaps(alignment)

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

        consensus_seq = calculate_consensus(alignment)
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
        return download_file_temporary_PLT("sequence_alignment_plot")
    except KeyError as e:
        show_error_dialog(f"Key Error: {e}")
    except Exception as e:
        show_error_dialog(f"An unexpected error occurred: {e}")
        
def calculate_conservation_and_gaps(alignment):
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

    return [], []

def calculate_consensus(alignment):
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
