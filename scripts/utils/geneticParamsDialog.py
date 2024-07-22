import sys
from PyQt5.QtWidgets import QApplication, QWidget, QFormLayout, QLabel, QLineEdit, QComboBox, QPushButton

class ParamDialog(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        layout = QFormLayout()

        # Bootstrap threshold
        self.bootstrap_threshold_label = QLabel('Bootstrap Threshold:')
        self.bootstrap_threshold_input = QLineEdit('10')
        layout.addRow(self.bootstrap_threshold_label, self.bootstrap_threshold_input)

        # Window size
        self.window_size_label = QLabel('Window Size:')
        self.window_size_input = QLineEdit('200')
        layout.addRow(self.window_size_label, self.window_size_input)

        # Step size
        self.step_size_label = QLabel('Step Size:')
        self.step_size_input = QLineEdit('100')
        layout.addRow(self.step_size_label, self.step_size_input)

        # Bootstrap amount
        self.bootstrap_amount_label = QLabel('Bootstrap Amount:')
        self.bootstrap_amount_input = QLineEdit('100')
        layout.addRow(self.bootstrap_amount_label, self.bootstrap_amount_input)

        # Alignment method
        self.alignment_method_label = QLabel('Alignment Method:')
        self.alignment_method_input = QComboBox()
        self.alignment_method_input.addItems(['pairwiseAligner', 'MUSCLE', 'CLUSTALW', 'MAFFT'])
        self.alignment_method_input.setCurrentIndex(1) # Default to MUSCLE
        layout.addRow(self.alignment_method_label, self.alignment_method_input)

        # Fit method
        self.fit_method_label = QLabel('Fit Method:')
        self.fit_method_input = QComboBox()
        self.fit_method_input.addItems(['Wider Fit by elongating with Gap (starAlignment)', 'Narrow-fit prevent elongation with gap when possible'])
        self.fit_method_input.setCurrentIndex(0) # Default to Wider Fit
        layout.addRow(self.fit_method_label, self.fit_method_input)

        # Tree type
        self.tree_type_label = QLabel('Tree Type:')
        self.tree_type_input = QComboBox()
        self.tree_type_input.addItems(['BioPython consensus tree', 'FastTree application'])
        self.tree_type_input.setCurrentIndex(1) # Default to FastTree application
        layout.addRow(self.tree_type_label, self.tree_type_input)

        # Rate similarity
        self.rate_similarity_label = QLabel('Rate Similarity:')
        self.rate_similarity_input = QLineEdit('90')
        layout.addRow(self.rate_similarity_label, self.rate_similarity_input)

        # Method similarity
        self.method_similarity_label = QLabel('Method Similarity:')
        self.method_similarity_input = QComboBox()
        self.method_similarity_input.addItems([
            'Hamming distance',
            'Levenshtein distance',
            'Damerau-Levenshtein distance',
            'Jaro similarity',
            'Jaro-Winkler similarity',
            'Smith–Waterman similarity',
            'Jaccard similarity',
            'Sørensen-Dice similarity'
        ])
        self.method_similarity_input.setCurrentIndex(0) # Default to Hamming distance
        layout.addRow(self.method_similarity_label, self.method_similarity_input)

        # Submit button
        self.submit_button = QPushButton('Submit')
        self.submit_button.clicked.connect(self.submit)
        layout.addRow(self.submit_button)

        self.setLayout(layout)
        self.setWindowTitle('Parameter Dialog')

        # Set style
        self.setStyleSheet("""
            QWidget {
                background-color: #f0f0f0;
                font-family: Arial, sans-serif;
            }
            QLabel {
                color: #333;
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
        """)

        self.resize(400, 300)
        self.show()

    def submit(self):
        params = {
            'bootstrap_threshold': int(self.bootstrap_threshold_input.text()),
            'window_size': int(self.window_size_input.text()),
            'step_size': int(self.step_size_input.text()),
            'bootstrap_amount': int(self.bootstrap_amount_input.text()),
            'alignment_method': str(self.alignment_method_input.currentIndex() + 1),
            'fit_method': str(self.fit_method_input.currentIndex() + 1),
            'tree_type': str(self.tree_type_input.currentIndex() + 1),
            'rate_similarity': int(self.rate_similarity_input.text()),
            'method_similarity': str(self.method_similarity_input.currentIndex() + 1),
        }
        print(params)  # Here you can handle the parameters as needed
        self.close()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = ParamDialog()
    sys.exit(app.exec_())
