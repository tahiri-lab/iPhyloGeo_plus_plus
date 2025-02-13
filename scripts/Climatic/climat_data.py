import io
import folium

from utils.error_dialog import show_error_dialog
    
        
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