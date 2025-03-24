from PyQt6 import QtCore, QtWidgets
from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QMovie

from ui import loading_dialog

from typing import Any, Dict
from aphylogeo.params import Params

from utils.file_caching import FileCaching
from utils.error_dialog import show_error_dialog

from worker import Worker

from Genetics.genetic_params_dialog import ParamDialog
from Genetics.genetics_alignment_plot import plot_alignment_chart, read_msa, standardize_sequence_lengths

class GeneticAlignment:
    def __init__(self, main):
        self.main = main
        self.msa = {}
        self.geneticTrees = []
        self.worker = None
        
    def show_sequence_alignment_page(self):
        """
        Display the sequence alignment page.

        This method sets the stacked widget's current index to 1 and the tab widget's current index to 2
        to display the sequence alignment page.
        """
        
        self.main.stackedWidget.setCurrentIndex(1)
        self.main.tabWidget.setCurrentIndex(2)
        
    def update_plot(self):
        """
        Update the plot based on the current starting position and window size.

        This method reads the MSA data, standardizes the sequence lengths, and plots the alignment chart.
        The resulting plot is displayed in the specified label widget and the tab widget is updated.

        Returns:
            None
        """
        try:
            starting_position = self.main.starting_position_spinbox_2.value()
            window_size = self.main.window_size_spinbox_2.value()

            standardized_data = standardize_sequence_lengths(self.msa)
            pixmap = plot_alignment_chart(standardized_data, starting_position, window_size, self.main.isDarkMode)

            self.main.seqAlignLabel.setPixmap(pixmap)
            self.main.tabWidget.setCurrentIndex(2)

        except AttributeError as e:
            show_error_dialog(f"Attribute Error: {e}")
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")
            
            
    def start_alignment_analysis(self):
        """
        Perform sequence alignment and store the resulting genetic tree dictionary.

        This method calls the `callSeqAlign` method to perform sequence alignment and stores the resulting genetic tree dictionary in the `geneticTrees` attribute.
        """
        self.main.starting_position_spinbox_2.setEnabled(True)
        self.main.window_size_spinbox_2.setEnabled(True)
        self.call_seq_align()
        
        
    def call_seq_align(self):
        """
        Execute the sequence alignment pipeline and display progress using a worker thread.

        This method performs the following steps:
        1. Loads sequences from the reference gene file.
        2. Aligns the loaded sequences.
        3. Generates genetic trees based on the aligned sequences.
        4. Prepares and displays the results in the UI.
        5. Saves the alignment and genetic tree results to JSON files.

        Returns:
            None
        """

        def update_progress(step):
            if step < loading_screen.checkListWidget.count():
                item = loading_screen.checkListWidget.item(step)
                item.setCheckState(Qt.CheckState.Checked)
                progress_value = int((step + 1) * (100 / loading_screen.checkListWidget.count()))
                loading_screen.progressBar.setValue(progress_value)
                QtWidgets.QApplication.processEvents()
            else:
                loading_screen.progressBar.setValue(100)
                QtWidgets.QApplication.processEvents()

        def handle_finished(result: Dict[str, Any]):
            loading_screen.close()
            msa = result["msa"]
            self.msa = read_msa(msa)
            self.geneticTrees = result["geneticTrees"]
            self.update_plot()
            self.main.geneticTreeButtonPage1.setEnabled(True)
            self.main.statisticsButtonPage1.setEnabled(True)
            if self.main.climaticTreeButtonPage2.isEnabled():
                self.main.resultsButton.setEnabled(True)

        def handle_error(error_message):
            loading_screen.close()
            show_error_dialog(f"An unexpected error occurred: {error_message}")

        if loading_screen := loading_dialog.LoadingDialog():
            loading_screen.setWindowFlags(Qt.WindowType.FramelessWindowHint)  # Remove the title bar and frame

            loading_screen.setWindowModality(Qt.WindowModality.ApplicationModal)

            # Set the QMovie for the movieLabel
            movie = QMovie(":active/dna.gif")  # Use the resource path for the gif
            loading_screen.movieLabel.setMovie(movie)

            # Resize the movie to fit within the QLabel
            movie.setScaledSize(QtCore.QSize(100, 100))  # Set the desired size here

            # Ensure the QLabel is centered and the GIF is properly displayed
            loading_screen.movieLabel.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            movie.start()

            # Show the loading screen
            loading_screen.show()
        QtWidgets.QApplication.processEvents()

        if (result := FileCaching.get_cached_result_file(Params.reference_gene_filepath)) is not None:
            print("Using cached result")
            handle_finished(result)
            return

        self.workerThread = QThread()
        self.worker = Worker(Params.reference_gene_filepath)
        self.worker.moveToThread(self.workerThread)

        self.worker.progress.connect(update_progress)
        self.worker.finished.connect(handle_finished)
        self.worker.error.connect(handle_error)

        self.workerThread.started.connect(self.worker.run)
        self.worker.finished.connect(self.workerThread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.workerThread.finished.connect(self.workerThread.deleteLater)

        self.workerThread.start()

        # Use a loop to wait until the thread finishes and the result is set
        while self.workerThread.isRunning():
            QtWidgets.QApplication.processEvents()
            
    def stopWorker(self):
        if self.worker:
            self.worker.stop()
            
            
def open_genetic_settings_window():
    dialog = ParamDialog()
    dialog.exec()