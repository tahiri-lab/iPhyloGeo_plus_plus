from PyQt6.QtGui import QIcon
from utils.error_dialog import show_error_dialog
from utils.help import HelpDialog


class Navigation:
    def __init__(self, main):
        self.main = main

    def show_home_section(self):
        """
        Display the home page of the application.

        This method sets the icons for the climatic data and genetic data buttons to their inactive states
        and displays the home page by setting the stacked widget's current index to 0.
        """
        try:
            self.main.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.main.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.main.homeButton.setIcon(QIcon(":active/home.png"))
            self.main.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.main.stackedWidget.setCurrentIndex(0)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def show_genetic_section(self):
        """
        Display the genetic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the genetic data page
        by setting the stacked widget's current index to 1, and sets the tab widget's current index to 0.
        """
        try:
            self.main.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.main.geneticDataButton.setIcon(QIcon(":active/genetic.svg"))
            self.main.homeButton.setIcon(QIcon(":other/home.svg"))
            self.main.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.main.stackedWidget.setCurrentIndex(1)
            self.main.tabWidget.setCurrentIndex(0)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def show_climate_section(self):
        """
        Display the climatic data page of the application.

        This method sets the icons for the climatic data and genetic data buttons, displays the climatic data page
        by setting the stacked widget's current index to 2, and sets the tab widget's current index to 0.
        """
        try:
            self.main.climaticDataButton.setIcon(QIcon(":active/climaticData.png"))
            self.main.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.main.homeButton.setIcon(QIcon(":other/home.svg"))
            self.main.resultsButton.setIcon(QIcon(":inactive/result.svg"))
            self.main.stackedWidget.setCurrentIndex(2)
            self.main.tabWidget2.setCurrentIndex(0)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def show_results_section(self):
        """
        Display the results page of the application.

        This method sets the stacked widget's current index to 3 to display the results page.
        """
        try:
            self.main.climaticDataButton.setIcon(QIcon(":inactive/climaticData.svg"))
            self.main.geneticDataButton.setIcon(QIcon(":inactive/genetic.svg"))
            self.main.homeButton.setIcon(QIcon(":other/home.svg"))
            self.main.resultsButton.setIcon(QIcon(":active/result.svg"))
            self.main.stackedWidget.setCurrentIndex(3)
        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    def open_help_window(self):
        """
        Initialize and display the 'How to Use' window.

        This method creates a new QMainWindow instance, sets up its UI using the UiHowToUse class, and displays the window.
        """

        dialog = HelpDialog()
        dialog.exec()
