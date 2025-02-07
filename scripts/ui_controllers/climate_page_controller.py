from event_connector import QtEvents, connect_decorated_methods_parent, connect_event


class ClimatePageController:
    def __init__(self, main_window) -> None:
        self.main_window = main_window
        connect_decorated_methods_parent(self, self.main_window)

    @connect_event("statisticsButtonPage2", QtEvents.clicked)
    def load_climate_statistics_click(self):
        self.main_window.climat.load_climate_statistics()

    @connect_event(["ClimaticChartSettingsAxisX", "ClimaticChartSettingsAxisY", "PlotTypesCombobox"], QtEvents.currentIndexChanged)
    def generate_climate_graph_click(self):
        self.main_window.climat.generate_climate_graph()

    @connect_event("climatePlotDownloadButton", QtEvents.clicked)
    def download_climate_plot_click(self):
        self.main_window.climat.download_climate_plot()

    @connect_event("fileBrowserButtonPage2", QtEvents.clicked)
    def load_csv_climate_file_click(self):
        self.main_window.climat.load_csv_climate_file()

    @connect_event("clearButtonPage2", QtEvents.clicked)
    def clear_climmatic_data_click(self):
        self.main_window.climat.clear_climmatic_data()

    @connect_event("climaticTreeButtonPage2", QtEvents.clicked)
    def display_climatic_trees_click(self):
        self.main_window.climat.climaticTree.display_climatic_trees()

    @connect_event("downloadGraphButton2", QtEvents.clicked)
    def download_climatic_tree_graph_click(self):
        self.main_window.climat.climaticTree.download_climatic_tree_graph()

    @connect_event("preferencesButton", QtEvents.clicked)
    def open_climatic_tree_preferences_window_click(self):
        self.main_window.climat.climaticTree.open_climatic_tree_preferences_window()

    @connect_event("climaticTreescomboBox", QtEvents.currentIndexChanged)
    def show_selected_climatic_tree_click(self, index):
        self.main_window.climat.climaticTree.show_climatic_tree(index)
