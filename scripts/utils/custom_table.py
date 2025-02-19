import re
from decimal import Decimal

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QHeaderView, QTableWidget, QTableWidgetItem


def create_sleek_table(df, hardCodeSize=False):
    num_rows, num_columns = df.shape
    table_widget = QTableWidget(num_rows, num_columns)
    table_widget.setStyleSheet(
        """
        QTableWidget {
            border: 1px solid #ddd;
        }
        QHeaderView::section {
            background-color: #4CAF50;
            color: white;
            font-weight: bold;
            border: none;
            padding: 5px;
        }
        QTableWidget::item {
            border: none;
            padding: 5px;
            font-weight: bold;
        }
    """
    )
    table_widget.setAlternatingRowColors(True)

    if horizontal_header := table_widget.horizontalHeader():
        horizontal_header.setStretchLastSection(True)
        horizontal_header.setVisible(True)
        horizontal_header.setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)
        if hardCodeSize:
            horizontal_header.setDefaultSectionSize(150)
        else:
            horizontal_header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)

    if vertical_header := table_widget.verticalHeader():
        vertical_header.setVisible(False)

    for col in range(num_columns):
        item = QTableWidgetItem(transform_string_to_fit_table(df.columns[col]))
        item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        table_widget.setHorizontalHeaderItem(col, item)

    for row in range(num_rows):
        for col in range(num_columns):
            value = transform_string_to_fit_table(str(df.iloc[row, col]))
            item = QTableWidgetItem(value)
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            table_widget.setItem(row, col, item)

    return table_widget


def transform_string_to_fit_table(value: str):
    if value.__contains__(".") and value.split(".")[1].__len__() >= 3 and re.search("^[0-9\\-]*\\.[0-9]*", value) is not None:
        value = str(round(Decimal(value), 3))
    elif len(value) > 10 and value.__contains__(" "):
        value = value.replace(" ", "\n", 1)
    return value
