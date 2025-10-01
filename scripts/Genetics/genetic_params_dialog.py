import sys

from aphylogeo.params import Params
from PyQt6.QtGui import QIntValidator
from PyQt6.QtWidgets import QApplication, QComboBox, QDialog, QFormLayout, QLabel, QLineEdit, QPushButton
from utils.file_caching import FileCaching
from utils.my_dumper import update_yaml_param

try:
    Params.load_from_file("params.yaml")
except FileNotFoundError:
    Params.validate_and_set_params(Params.PARAMETER_KEYS)


class ParamDialog(QDialog):
    def __init__(self):
        super().__init__()

        self.params = None  # Initialize the parameters to None
        self.initUI()

    def initUI(self):
        layout = QFormLayout()

        # Bootstrap threshold
        self.bootstrap_threshold_label = QLabel("Bootstrap Threshold:")
        self.bootstrap_threshold_input = QLineEdit(str(Params.bootstrap_threshold))
        self.bootstrap_threshold_input.setValidator(QIntValidator(0, 10000))  # Adjust range as necessary
        layout.addRow(self.bootstrap_threshold_label, self.bootstrap_threshold_input)

        # Window size
        self.window_size_label = QLabel("Window Size:")
        self.window_size_input = QLineEdit(str(Params.window_size))
        self.window_size_input.setValidator(QIntValidator(1, 10000))  # Adjust range as necessary
        layout.addRow(self.window_size_label, self.window_size_input)

        # Step size
        self.step_size_label = QLabel("Step Size:")
        self.step_size_input = QLineEdit(str(Params.step_size))
        self.step_size_input.setValidator(QIntValidator(1, 10000))  # Adjust range as necessary
        layout.addRow(self.step_size_label, self.step_size_input)

        # Bootstrap amount
        self.bootstrap_amount_label = QLabel("Bootstrap Amount:")
        self.bootstrap_amount_input = QLineEdit(str(Params.bootstrap_amount))
        self.bootstrap_amount_input.setValidator(QIntValidator(1, 10000))  # Adjust range as necessary
        layout.addRow(self.bootstrap_amount_label, self.bootstrap_amount_input)

        # Alignment method
        self.alignment_method_label = QLabel("Alignment Method:")
        self.alignment_method_input = QComboBox()
        self.alignment_method_input.addItems(["pairwiseAligner", "MUSCLE", "CLUSTALW", "MAFFT"])
        self.alignment_method_input.setCurrentIndex(int(Params.alignment_method) - 1)  # Default to MUSCLE
        layout.addRow(self.alignment_method_label, self.alignment_method_input)

        # Fit method
        self.fit_method_label = QLabel("Fit Method:")
        self.fit_method_input = QComboBox()
        self.fit_method_input.addItems([
            "Wider Fit by elongating with Gap (starAlignment)",
            "Narrow-fit prevent elongation with gap when possible",
        ])
        self.fit_method_input.setCurrentIndex(int(Params.fit_method) - 1)  # Default to Wider Fit
        layout.addRow(self.fit_method_label, self.fit_method_input)

        # Tree type
        self.tree_type_label = QLabel("Tree Type:")
        self.tree_type_input = QComboBox()
        self.tree_type_input.addItems(["BioPython consensus tree", "FastTree application"])
        self.tree_type_input.setCurrentIndex(int(Params.tree_type) - 1)  # Default to FastTree application
        layout.addRow(self.tree_type_label, self.tree_type_input)

        # Rate similarity
        self.rate_similarity_label = QLabel("Rate Similarity:")
        self.rate_similarity_input = QLineEdit(str(Params.rate_similarity))
        self.rate_similarity_input.setValidator(QIntValidator(0, 100))  # Rate similarity as percentage
        layout.addRow(self.rate_similarity_label, self.rate_similarity_input)

        # Method similarity
        self.method_similarity_label = QLabel("Method Similarity:")
        self.method_similarity_input = QComboBox()
        self.method_similarity_input.addItems([
            "Hamming distance",
            "Levenshtein distance",
            "Damerau-Levenshtein distance",
            "Jaro similarity",
            "Jaro-Winkler similarity",
            "Smith–Waterman similarity",
            "Jaccard similarity",
            "Sørensen-Dice similarity",
        ])
        self.method_similarity_input.setCurrentIndex(int(Params.method_similarity) - 1)  # Default to Hamming distance
        layout.addRow(self.method_similarity_label, self.method_similarity_input)

        # Submit button
        self.submit_button = QPushButton("Save Settings")
        self.submit_button.clicked.connect(self.submit)
        layout.addRow(self.submit_button)

        self.setLayout(layout)
        self.setWindowTitle("Parameter Dialog")

        # Set style
        self.setStyleSheet(
            """
            QWidget {
                font-family: Arial, sans-serif;
            }
            QLabel {
                font-size: 14px;
            }
            QLineEdit {
                border: 1px solid #ccc;
                border-radius: 5px;
                padding: 5px;
                font-size: 14px;
            }
            QComboBox {
                border: 1px solid #ccc;
                border-radius: 5px;
                padding: 5px;
                font-size: 14px;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 10px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """
        )

        self.resize(400, 300)
        self.show()

    def submit(self):
        self.params = {
            "bootstrap_threshold": int(self.bootstrap_threshold_input.text()),
            "window_size": int(self.window_size_input.text()),
            "step_size": int(self.step_size_input.text()),
            "bootstrap_amount": int(self.bootstrap_amount_input.text()),
            "alignment_method": str(self.alignment_method_input.currentIndex() + 1),
            "fit_method": str(self.fit_method_input.currentIndex() + 1),
            "tree_type": str(self.tree_type_input.currentIndex() + 1),
            "rate_similarity": int(self.rate_similarity_input.text()),
            "method_similarity": str(self.method_similarity_input.currentIndex() + 1),
        }
        for property_name, new_value in self.params.items():
            update_yaml_param(Params, "scripts/utils/params.yaml", property_name, new_value)
        FileCaching.clear_cache()
        self.accept()  # Close the dialog and indicate success


if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = ParamDialog()
    sys.exit(app.exec())
