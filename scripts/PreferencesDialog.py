from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QComboBox, QRadioButton, QCheckBox, QLineEdit, QPushButton

class PreferencesDialog(QDialog):
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

        # Edge Color
        self.edge_color_label = QLabel("Edge color:")
        self.edge_color_combo = QComboBox()
        self.edge_color_combo.addItems(["black", "blue", "red", "green"])

        # Reticulation Color
        self.reticulation_color_label = QLabel("Reticulation color:")
        self.reticulation_color_combo = QComboBox()
        self.reticulation_color_combo.addItems(["magenta", "cyan", "yellow", "black"])

        # Layout Options
        self.layout_options_label = QLabel("Layout Options:")
        self.vertical_radio = QRadioButton("Hierarchical Vertical")
        self.horizontal_radio = QRadioButton("Hierarchical Horizontal")
        self.horizontal_radio.setChecked(True)
        self.axial_radio = QRadioButton("Axial")
        self.radial_radio = QRadioButton("Radial")

        # Other Options
        self.proportional_edge_lengths = QCheckBox("Proportional edge lengths")
        self.label_internal_vertices = QCheckBox("Label internal vertices")
        self.use_leaf_names = QCheckBox("Use leaf names")
        self.use_leaf_names.setChecked(True)
        self.root_with_leaf_node = QCheckBox("Root with leaf/node")
        self.root_leaf_node_input = QLineEdit("0")

        # Save Button
        self.save_button = QPushButton("Save")
        self.save_button.clicked.connect(self.accept)

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
        layout.addWidget(self.proportional_edge_lengths)
        layout.addWidget(self.label_internal_vertices)
        layout.addWidget(self.use_leaf_names)
        layout.addWidget(self.root_with_leaf_node)
        layout.addWidget(self.root_leaf_node_input)
        layout.addWidget(self.save_button)

        self.setLayout(layout)

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

        self.proportional_edge_lengths.setChecked(preferences.get("proportional_edge_lengths", False))
        self.label_internal_vertices.setChecked(preferences.get("label_internal_vertices", False))
        self.use_leaf_names.setChecked(preferences.get("use_leaf_names", True))
        self.root_with_leaf_node.setChecked(preferences.get("root_with_leaf_node", False))
        self.root_leaf_node_input.setText(preferences.get("root_leaf_node", "0"))

    def get_preferences(self):
        return {
            "label_color": self.label_color_combo.currentText(),
            "edge_color": self.edge_color_combo.currentText(),
            "reticulation_color": self.reticulation_color_combo.currentText(),
            "layout": "vertical" if self.vertical_radio.isChecked() else
                      "horizontal" if self.horizontal_radio.isChecked() else
                      "axial" if self.axial_radio.isChecked() else "radial",
            "proportional_edge_lengths": self.proportional_edge_lengths.isChecked(),
            "label_internal_vertices": self.label_internal_vertices.isChecked(),
            "use_leaf_names": self.use_leaf_names.isChecked(),
            "root_with_leaf_node": self.root_with_leaf_node.isChecked(),
            "root_leaf_node": self.root_leaf_node_input.text(),
        }

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ui = PreferencesDialog()
    ui.show()
    sys.exit(app.exec_())
