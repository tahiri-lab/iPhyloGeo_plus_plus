import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel, QScrollArea, QTextBrowser
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt

class UiHowToUse(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Help - aPhyloGeo')
        self.setGeometry(100, 100, 800, 600)
        self.initUI()

    def initUI(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        layout.setContentsMargins(20, 20, 20, 20)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)

        content_widget = QWidget()
        scroll_area.setWidget(content_widget)
        content_layout = QVBoxLayout(content_widget)
        content_layout.setContentsMargins(20, 20, 20, 20)

        title_label = QLabel('Help - aPhyloGeo')
        title_label.setFont(QFont('Helvetica', 28, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet('color: #2C3E50; margin-bottom: 20px;')
        content_layout.addWidget(title_label)

        sections = [
            {
                'title': 'Genetic Data',
                'description': 'To get the genetic information contained in a file.',
                'instructions': [
                    'Click on "Genetic Data" then on "File Browser" and select the file you want to access (.fasta extension recommended).',
                    'Once the file is selected, you will see a colorful representation of the genetic information contained in the file.'
                ]
            },
            {
                'title': 'Sequence Alignment',
                'description': 'To create the sequence alignment of the obtained data, click on "Sequence Alignment".',
                'instructions': [
                    'Get the "Genetic Data" of your file. This step may take some time and resources.',
                    'Once the alignment is completed, you will see a line for each variant whose RNA sequence has been aligned.'
                ]
            },
            {
                'title': 'Climatic Data',
                'description': 'To get the climatic data extracted from a spreadsheet.',
                'instructions': [
                    'Click on "Climatic Data" then on "File Browser" and select the file you want to access (.csv extension recommended).',
                    'Once the file is selected, you will see an array of the values retrieved. If the last 2 columns of the file are \'LAT\' and \'LONG\', a map will be rendered.'
                ]
            },
            {
                'title': 'Statistics',
                'description': 'To get the climatic bar graph of the data extracted.',
                'instructions': [
                    'Get the "Climatic Data" of your file.',
                    'Once the data is extracted, click on "Statistics".'
                ]
            },
            {
                'title': 'Results',
                'description': 'To get the results extracted from your file.',
                'instructions': [
                    'Click on "Results" then on "Submit" to filter the data according to a specific threshold.',
                    'Click on "Settings" to choose the filter you want to apply. The results are automatically saved in a file called "output.csv".'
                ]
            }
        ]

        for section in sections:
            section_widget = QWidget()
            section_layout = QVBoxLayout(section_widget)
            section_layout.setContentsMargins(0, 0, 0, 20)

            section_title = QLabel(section['title'])
            section_title.setFont(QFont('Helvetica', 22, QFont.Bold))
            section_title.setStyleSheet('color: #2C3E50; margin-bottom: 10px;')
            section_title.setWordWrap(True)
            section_layout.addWidget(section_title)

            section_description = QLabel(section['description'])
            section_description.setFont(QFont('Helvetica', 16))
            section_description.setStyleSheet('color: #7F8C8D; margin-bottom: 10px;')
            section_description.setWordWrap(True)
            section_layout.addWidget(section_description)

            for instruction in section['instructions']:
                instruction_label = QLabel(f'- {instruction}')
                instruction_label.setFont(QFont('Helvetica', 14))
                instruction_label.setStyleSheet('color: #2C3E50; margin-left: 20px;')
                instruction_label.setWordWrap(True)
                section_layout.addWidget(instruction_label)

            content_layout.addWidget(section_widget)

        additional_info = QTextBrowser()
        additional_info.setOpenExternalLinks(True)
        additional_info.setHtml(
            'You will obtain more information here: <a href="https://github.com/tahiri-lab/aPhyloGeo_plus_plus/blob/main/README.md">README.md</a>'
        )
        additional_info.setStyleSheet('margin-top: 20px; color: #2980B9;')
        content_layout.addWidget(additional_info)

        content_layout.addStretch()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = UiHowToUse()
    window.show()
    sys.exit(app.exec_())
