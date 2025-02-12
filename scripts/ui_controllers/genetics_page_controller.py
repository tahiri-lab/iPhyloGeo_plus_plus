from event_connector import QtEvents, connect_decorated_methods_parent, connect_event
from Genetics.genetics import Genetics, open_genetic_settings_window


class GeneticPageController:
    def __init__(self, main_window) -> None:
        self.genetics = Genetics(main_window)
        connect_decorated_methods_parent(self, main_window)

    @connect_event("geneticTreeButtonPage1", QtEvents.clicked)
    def display_newick_trees_click(self):
        self.genetics.geneticTree.display_newick_trees()

    @connect_event(["starting_position_spinbox_2", "window_size_spinbox_2"], QtEvents.valueChanged)
    def update_plot_click(self):
        self.genetics.update_plot()

    @connect_event("statisticsButtonPage1", QtEvents.clicked)
    def initialize_species_list_click(self):
        self.genetics.initialize_species_list()

    @connect_event("downloadSimilarityButton", QtEvents.clicked)
    def download_similarity_plot_chart_click(self):
        self.genetics.download_similarity_plot_chart()

    @connect_event("fileBrowserButtonPage1", QtEvents.clicked)
    def select_fasta_file_click(self):
        self.genetics.select_fasta_file()

    @connect_event("StartSequenceAlignmentButton", QtEvents.clicked)
    def start_alignment_analysis_click(self):
        self.genetics.start_alignment_analysis()

    @connect_event("sequenceAlignmentButtonPage1", QtEvents.clicked)
    def show_sequence_alignment_page_click(self):
        self.genetics.show_sequence_alignment_page()

    @connect_event("clearButtonPage1", QtEvents.clicked)
    def clear_genetic_data_click(self):
        self.genetics.clear_genetic_data()

    @connect_event("geneticTreescomboBox", QtEvents.currentIndexChanged)
    def show_tree_click(self, index):
        self.genetics.geneticTree.show_tree(index)

    @connect_event("downloadGraphButton", QtEvents.clicked)
    def download_genetic_tree_graph_click(self):
        self.genetics.geneticTree.download_genetic_tree_graph()

    @connect_event("geneticSettingsButton", QtEvents.clicked)
    def open_genetic_settings_window_click(self):
        open_genetic_settings_window()

    @connect_event(["similarityWindowSizeSpinBox", "startingPositionSimilaritySpinBox"], QtEvents.valueChanged)
    def update_similarity_plot_valueChanged(self):
        self.genetics.update_similarity_plot()

    @connect_event(["referenceComboBox"], QtEvents.currentIndexChanged)
    def update_similarity_plot_currentIndexChanged(self):
        self.genetics.update_similarity_plot()
