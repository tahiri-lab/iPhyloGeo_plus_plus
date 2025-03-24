import numpy as np
import pandas as pd
import plotly.express as px

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from event_connector import blocked_signals
from utils.error_dialog import show_error_dialog
from utils.download_file import download_local_with_fig

class GeneticStatistics:
    def __init__(self, main, alignment):
        self.main = main
        self.geneticAlignment = alignment
        self.customConfig ={
                'modeBarButtonsToRemove': ['toImage', 'autoScale'],
                'displaylogo' : False
            }
        
        
    def initialize_species_list(self):
        unique_species = set(self.geneticAlignment.msa.keys())

        with blocked_signals(self.main.referenceComboBox):
            self.main.referenceComboBox.clear()
            self.main.referenceComboBox.addItems(unique_species)
            self.main.referenceComboBox.setCurrentIndex(0)

        self.update_similarity_plot()     

    def update_similarity_plot(self):
        try:
            reference_species = self.main.referenceComboBox.currentText().replace(" ", "_")  # Convert back to original format

            sequences = self.geneticAlignment.msa

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