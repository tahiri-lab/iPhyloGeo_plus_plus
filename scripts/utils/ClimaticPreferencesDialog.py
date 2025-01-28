from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QLabel, QComboBox, QRadioButton, QCheckBox, QPushButton, QButtonGroup
from utils.ClimaticGraphSettings import ClimaticGraphSettings
from utils.my_dumper import update_yaml_param

try:
    ClimaticGraphSettings.load_from_file("./scripts/utils/ClimaticGraphSettings.yaml")
except FileNotFoundError:
    ClimaticGraphSettings.validate_and_set_params(ClimaticGraphSettings.PARAMETER_KEYS)
    

class ClimaticPreferencesDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        # Label Color
        self.label_color_label = QLabel("Label color:")
        self.label_color_combo = QComboBox()
        self.label_color_combo.addItems(["black", "blue", "red", "green"])
        self.label_color_combo.setCurrentText(ClimaticGraphSettings.label_color)

        # Edge Color
        self.edge_color_label = QLabel("Edge color:")
        self.edge_color_combo = QComboBox()
        self.edge_color_combo.addItems(["black", "blue", "red", "green"])
        self.edge_color_combo.setCurrentText(ClimaticGraphSettings.edge_color)

        # Reticulation Color
        self.reticulation_color_label = QLabel("Reticulation color:")
        self.reticulation_color_combo = QComboBox()
        self.reticulation_color_combo.addItems(["magenta", "cyan", "yellow", "black"])
        self.reticulation_color_combo.setCurrentText(ClimaticGraphSettings.reticulation_color)
        
        # Layout Options
        self.layout_options_label = QLabel("Layout Options:")
        self.vertical_radio = QRadioButton("Hierarchical Vertical")
        self.horizontal_radio = QRadioButton("Hierarchical Horizontal")
        self.axial_radio = QRadioButton("Axial")
        self.radial_radio = QRadioButton("Radial")
        
        match ClimaticGraphSettings.layout:
            case "vertical":
                self.vertical_radio.setChecked(True)
            case "horizontal":
                self.horizontal_radio.setChecked(True)
            case "axial":
                self.axial_radio.setChecked(True)
            case "radial":
                self.radial_radio.setChecked(True)
                

        # View Type Options
        self.view_type_label = QLabel("View Type:")
        self.network_view_radio = QRadioButton("Network View")
        self.tree_view_radio = QRadioButton("Tree View")
        match ClimaticGraphSettings.view_type:
            case "network":
                self.network_view_radio.setChecked(True)
            case "tree":
                self.tree_view_radio.setChecked(True)


        # Add view type radios to a button group to ensure only one can be selected at a time
        self.view_type_group = QButtonGroup()
        self.view_type_group.addButton(self.network_view_radio)
        self.view_type_group.addButton(self.tree_view_radio)

        # Other Options
        self.proportional_edge_lengths = QCheckBox("Proportional edge lengths")
        self.proportional_edge_lengths.setChecked(ClimaticGraphSettings.proportional_edge_lengths)
        self.label_internal_vertices = QCheckBox("Label internal vertices")
        self.label_internal_vertices.setChecked(ClimaticGraphSettings.label_internal_vertices)
        self.use_leaf_names = QCheckBox("Use leaf names")
        self.use_leaf_names.setChecked(ClimaticGraphSettings.use_leaf_names)
        self.show_branch_length = QCheckBox("Show branch lengths")
        self.show_branch_length.setChecked(ClimaticGraphSettings.show_branch_length)

        # Save Button
        self.save_button = QPushButton("Save")
        self.save_button.clicked.connect(self.submit)

        # Adding widgets to the layout
        layout.addWidget(self.label_color_label)
        layout.addWidget(self.label_color_combo)
        layout.addWidget(self.edge_color_label)
        layout.addWidget(self.edge_color_combo)
        layout.addWidget(self.reticulation_color_label)
        layout.addWidget(self.reticulation_color_combo)
        layout.addWidget(self.layout_options_label)
        layout.addWidget(self.vertical_radio)
        layout.addWidget(self.horizontal_radio)
        layout.addWidget(self.axial_radio)
        layout.addWidget(self.radial_radio)
        layout.addWidget(self.view_type_label)
        layout.addWidget(self.network_view_radio)
        layout.addWidget(self.tree_view_radio)
        layout.addWidget(self.proportional_edge_lengths)
        layout.addWidget(self.label_internal_vertices)
        layout.addWidget(self.use_leaf_names)
        layout.addWidget(self.show_branch_length)
        layout.addWidget(self.save_button)

        self.setLayout(layout)
        self.show()

    def update_preferences(self, preferences):
        self.label_color_combo.setCurrentText(preferences.get("label_color", "black"))
        self.edge_color_combo.setCurrentText(preferences.get("edge_color", "blue"))
        self.reticulation_color_combo.setCurrentText(preferences.get("reticulation_color", "red"))

        layout = preferences.get("layout", "horizontal")
        if layout == "vertical":
            self.vertical_radio.setChecked(True)
        elif layout == "horizontal":
            self.horizontal_radio.setChecked(True)
        elif layout == "axial":
            self.axial_radio.setChecked(True)
        elif layout == "radial":
            self.radial_radio.setChecked(True)

        view_type = preferences.get("view_type", "network")
        if view_type == "network":
            self.network_view_radio.setChecked(True)
        else:
            self.tree_view_radio.setChecked(True)

        self.proportional_edge_lengths.setChecked(preferences.get("proportional_edge_lengths", False))
        self.label_internal_vertices.setChecked(preferences.get("label_internal_vertices", False))
        self.use_leaf_names.setChecked(preferences.get("use_leaf_names", True))
        self.show_branch_length.setChecked(preferences.get("show_branch_length", False))

    def get_preferences(self):
        return {
            "label_color": self.label_color_combo.currentText(),
            "edge_color": self.edge_color_combo.currentText(),
            "reticulation_color": self.reticulation_color_combo.currentText(),
            "layout": "vertical"
            if self.vertical_radio.isChecked()
            else "horizontal"
            if self.horizontal_radio.isChecked()
            else "axial"
            if self.axial_radio.isChecked()
            else "radial",
            "view_type": "network" if self.network_view_radio.isChecked() else "tree",
            "proportional_edge_lengths": self.proportional_edge_lengths.isChecked(),
            "label_internal_vertices": self.label_internal_vertices.isChecked(),
            "use_leaf_names": self.use_leaf_names.isChecked(),
            "show_branch_length": self.show_branch_length.isChecked(),
        }
        
    def submit(self):
        
        for property_name, new_value in self.get_preferences().items():
            update_yaml_param(ClimaticGraphSettings, "scripts/utils/ClimaticGraphSettings.yaml", property_name, new_value)
        self.accept()  # Close the dialog and indicate success


if __name__ == "__main__":
    import sys
    
    app = QtWidgets.QApplication(sys.argv)
    ui = ClimaticPreferencesDialog()
    sys.exit(app.exec())
