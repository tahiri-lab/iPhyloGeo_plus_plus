from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QListWidget, QListWidgetItem, QSizePolicy, QScrollArea
from PyQt5.QtGui import QFont, QColor
from PyQt5.QtCore import Qt
import sys
class UiHowToUse(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Instructions')
        self.setGeometry(100, 100, 600, 800)

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget(scroll_area)
        scroll_area.setWidget(scroll_content)

        layout = QVBoxLayout(scroll_content)

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
            section_layout = QVBoxLayout()
            section_label = QLabel(section['title'])
            section_label.setFont(QFont('Arial', 16, QFont.Bold))
            section_label.setStyleSheet('color: #333; border-bottom: 2px solid #4CAF50; padding-bottom: 5px; text-transform: uppercase;')

            description_label = QLabel(section['description'])
            description_label.setFont(QFont('Arial', 12, QFont.StyleItalic))
            description_label.setStyleSheet('color: #666;')

            section_layout.addWidget(section_label)
            section_layout.addWidget(description_label)

            instruction_list = QListWidget()
            for instruction in section['instructions']:
                item = QListWidgetItem(f"➔ {instruction}")
                item.setFont(QFont('Arial', 12))
                item.setBackground(QColor('#e0f7fa'))
                item.setTextAlignment(Qt.AlignLeft)
                instruction_list.addItem(item)

            instruction_list.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            section_layout.addWidget(instruction_list)
            layout.addLayout(section_layout)

        footer_label = QLabel('You will obtain more information here: <a href="https://github.com/tahiri-lab/aPhyloGeo_plus_plus/blob/main/README.md">README.md</a>')
        footer_label.setOpenExternalLinks(True)
        footer_label.setAlignment(Qt.AlignCenter)
        footer_label.setStyleSheet('font-size: 0.9em; color: #777;')
        layout.addWidget(footer_label)

        main_layout = QVBoxLayout(self)
        main_layout.addWidget(scroll_area)
        self.setLayout(main_layout)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWin = UiHowToUse()
    mainWin.show()
    sys.exit(app.exec_())