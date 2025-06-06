from event_connector import QtEvents, connect_decorated_methods_parent, connect_event

from result import Result


class ResultPageController:
    def __init__(self, main_window) -> None:
        self.result = Result(main_window)
        connect_decorated_methods_parent(self, main_window)

    @connect_event("settingsButtonPage3", QtEvents.clicked)
    def open_result_settings_window_click(self):
        self.result.open_result_settings_window()

    @connect_event("submitButtonPage3", QtEvents.clicked)
    def show_filtered_results(self):
        self.result.show_filtered_results()

    @connect_event("clearButtonPage3", QtEvents.clicked)
    def clear_result_click(self):
        self.result.clear_result()

    @connect_event(["statisticsButtonPage3"], QtEvents.clicked)
    def display_phylogeographic_trees_click(self):
        self.result.display_phylogeographic_trees()

    @connect_event(["phyloTreescomboBox", "criteriaComboBox"], QtEvents.currentIndexChanged)
    def render_tree_click(self):
        self.result.render_tree()

    @connect_event("downloadResultsPlotButton", QtEvents.clicked)
    def save_tree_graph_click(self):
        self.result.save_tree_graph()

    @connect_event("tabWidgetResult", QtEvents.currentChanged)
    def on_result_tab_changed(self, index):
        self.result.on_tab_changed(index)
