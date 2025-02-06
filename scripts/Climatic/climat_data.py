import io
import re
import folium
from decimal import Decimal

from PyQt6.QtWidgets import QTableWidget, QTableWidgetItem
from PyQt6.QtCore import Qt

from utils.error_dialog import show_error_dialog

def create_sleek_table(df):
    num_rows, num_columns = df.shape
    table_widget = QTableWidget(num_rows, num_columns)
    table_widget.setStyleSheet(
        """
        QTableWidget {
            background-color: #f2f2f2;
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
        }
    """
    )
    if horizontal_header := table_widget.horizontalHeader():
        horizontal_header.setStretchLastSection(True)
        horizontal_header.setVisible(True)
        horizontal_header.setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)
        horizontal_header.setDefaultSectionSize(150)

    if vertical_header := table_widget.verticalHeader():
        vertical_header.setVisible(False)

    for col in range(num_columns):
        item = QTableWidgetItem(df.columns[col])
        item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        table_widget.setHorizontalHeaderItem(col, item)

    for row in range(num_rows):
        for col in range(num_columns):
            value = str(df.iloc[row, col])
            if re.search("^[0-9\\-]*\\.[0-9]*", value) is not None:
                value = str(round(Decimal(value), 3))
            item = QTableWidgetItem(value)
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            table_widget.setItem(row, col, item)

    return table_widget
            
        
def get_folium_data(lat, long):
    """
    Create and display a folium map with given latitude and longitude.

    This method generates a map centered on the mean latitude and longitude of the provided coordinates.
    It places markers on the map for each coordinate pair and displays the map in a QWebEngineView.

    Args:
        lat (list): List of latitudes.
        long (list): List of longitudes.
    """
    try:
        if not lat or not long or len(lat) != len(long):
            raise ValueError("Latitude and longitude lists must be non-empty and of the same length.")

        lat = [float(x) for x in lat]  # Convert all elements in lat to float
        long = [float(x) for x in long]  # Convert all elements in long to float

        mean_lat = sum(lat) / len(lat)
        mean_long = sum(long) / len(long)

        m = folium.Map(location=[mean_lat, mean_long], zoom_start=4, tiles="OpenStreetMap")  # Adjusted zoom level for better visibility

        # Add markers to the map
        for latitude, longitude in zip(lat, long):
            folium.Marker([latitude, longitude]).add_to(m)

        # Adjust the map to fit all markers
        m.fit_bounds([[min(lat), min(long)], [max(lat), max(long)]])

        data = io.BytesIO()
        m.save(data, close_file=False)

        return data.getvalue().decode()

    except ValueError as e:
        show_error_dialog(f"Value Error: {e}")
    except Exception as e:
        show_error_dialog(f"An unexpected error occurred: {e}")