from pandas import DataFrame
import plotly.express as px

from event_connector import blocked_signals
from utils.download_file import download_local_with_fig

class ClimaticStatistics:
    def __init__(self, main):
        self.main = main
        self.data = DataFrame
        self.customConfig ={
                'modeBarButtonsToRemove': ['toImage', 'autoScale', 'select2d', 'lasso2d'],
                'displaylogo' : False
            }
        
    def load_climate_statistics(self):
        """
        Load climate data from a CSV file, update the UI elements with the column names, and switch to the appropriate tab.

        This method reads the climate data from a CSV file (excluding the first column), updates combo boxes
        for X and Y axis settings with the column names, and switches to the climate data tab.

        Returns:
            None
        """        
        columns = self.data.columns.tolist()
        self.first_column_name = self.data.columns[0]

        if columns:
            columns.pop(0)

        with blocked_signals(self.main.ClimaticChartSettingsAxisX, self.main.ClimaticChartSettingsAxisY):
            self.main.ClimaticChartSettingsAxisX.clear()
            self.main.ClimaticChartSettingsAxisX.addItems(columns)
            self.main.ClimaticChartSettingsAxisY.clear()
            self.main.ClimaticChartSettingsAxisY.addItems(columns)
            self.main.ClimaticChartSettingsAxisY.setCurrentIndex(1)
            self.main.tabWidget2.setCurrentIndex(2)

        self.generate_climate_graph()

    def generate_climate_graph(self):
        """
        Generate and display a graph based on the selected X and Y axis data and the chosen plot type.

        This method reads the selected data columns and plot type from the UI, generates the corresponding graph,
        and displays it in the specified QLabel widget.

        Returns:
            None
        """
        x_data = self.main.ClimaticChartSettingsAxisX.currentText()
        y_data = self.main.ClimaticChartSettingsAxisY.currentText()
        plot_type = self.main.PlotTypesCombobox.currentText()
        
        self.main.ClimaticChartSettingsAxisX.setEnabled(True)
        self.main.ClimaticChartSettingsAxisY.setEnabled(True)
        
        match plot_type:
            case "Scatter Plot":
                self.generate_scatter_plot(x_data, y_data)
            case "Line Plot":
                self.generate_line_plot(x_data, y_data)
            case "Bar Graph":
                self.generate_bar_graph(x_data, y_data)
            case "Violin Plot":
                self.generate_violin_plot(x_data, y_data)
            case "Pie Plot":
                self.generate_pie_plot(x_data, y_data)
            case "Correlation":
                self.generate_correlation_plot()
                self.main.ClimaticChartSettingsAxisX.setEnabled(False)
                self.main.ClimaticChartSettingsAxisY.setEnabled(False)
                  
        self.main.climatGraphView.setHtml(self.fig.to_html(include_plotlyjs='cdn', config= self.customConfig))

    def generate_bar_graph(self, x_data, y_data):
        self.fig = px.bar(
        data_frame = self.data,
        x=x_data,
        y=y_data,
        hover_name= self.first_column_name
        )

    def generate_scatter_plot(self, x_data, y_data):
        self.fig = px.scatter(
            data_frame = self.data,
            x = x_data,
            y = y_data,
            hover_name = self.first_column_name
        )

    def generate_line_plot(self, x_data, y_data):
        self.fig = px.line(
            data_frame = self.data,
            x = x_data,
            y = y_data,
            hover_data=[self.first_column_name],
            markers=True
        )

    def generate_pie_plot(self, x_data, y_data):
        self.fig = px.pie(
        data_frame = self.data,
        names=self.data[self.first_column_name],
        values=y_data,
        hover_data=[x_data]
        )

    def generate_violin_plot(self, x_data, y_data):   
        self.fig = px.violin(
            data_frame = self.data,
            x = x_data,
            y = y_data,
            hover_name = self.first_column_name
        )
           
    def generate_correlation_plot(self):
        matrix_correlation = self.data.iloc[:, 1:].corr(method="spearman").round(2)
        self.fig = px.imshow(
            matrix_correlation,
            aspect="auto"
        )
         
    def download_climate_plot(self):
        """
        Download the generated plot.

        This method allows the user to download the currently displayed plot.

        Returns:
            None
        """
        download_local_with_fig(self.fig)