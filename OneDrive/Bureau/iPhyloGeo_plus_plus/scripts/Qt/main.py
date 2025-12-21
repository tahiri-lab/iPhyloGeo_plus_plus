import sys
import shutil

import qtmodern.styles
import qtmodern.windows
from aphylogeo.params import Params, os
from event_connector import QtEvents, connect_decorated_methods, connect_event
from navigation import Navigation
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import Qt, QThread, QSize
from PyQt6.QtGui import QColor, QIcon
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import QGraphicsDropShadowEffect, QVBoxLayout
from Qt import main_ui
from ui_controllers import ClimatePageController, GeneticPageController, ResultPageController
from ui_helpers import create_shadow_effect, get_button_style, style_buttons
from utils import resources_rc  # noqa: F401  # Import the compiled resource module for resolving image resource path
from utils.error_dialog import show_error_dialog

try:
    Params.load_from_file("params.yaml")
except FileNotFoundError:
    Params.validate_and_set_params(Params.PARAMETER_KEYS)



class UiMainWindow(main_ui.Ui_MainWindow, QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.navigation = Navigation(self)
        self.setupUi(self)
        self.setup_ui()
        self.geneticsPage = GeneticPageController(self)
        self.climatePage = ClimatePageController(self)
        self.resultPage = ResultPageController(self)
        connect_decorated_methods(self)

    def setup_ui(self):
        """
        Setup the UI components and initialize the main window.

        This method connects various UI buttons to their corresponding event handlers,
        sets up styles and effects for UI elements, and initializes the state of the application.
        """
        try:
            self.maplayout = QVBoxLayout(self.graphicsViewClimData)
            self.mapView = QWebEngineView(self.graphicsViewClimData)
            self.maplayout.addWidget(self.mapView)
            self.graphicsViewClimData.setLayout(self.maplayout)

            self.climatTableLayout = QVBoxLayout(self.textEditClimData)
            self.resultTableLayout = QVBoxLayout(self.textEditResults)

            self.setObjectName("MainWindow")
            self.window_size_spinbox_2.setRange(1, 1000)
            self.starting_position_spinbox_2.setRange(1, 1000)
            self.isDarkMode = False  # Keep track of the state
            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
            self.stackedWidget.setCurrentIndex(0)

            self.buttons = [
                self.geneticDataButton,
                self.climaticDataButton,
                self.helpButton,
                self.homeButton,
                self.resultsButton,
            ]

            self.buttons_vertical = [
                self.fileBrowserButtonPage1,
                self.sequenceAlignmentButtonPage1,
                self.clearButtonPage1,
                self.statisticsButtonPage1,
                self.geneticTreeButtonPage1,
                self.fileBrowserButtonPage2,
                self.clearButtonPage2,
                self.climaticTreeButtonPage2,
                self.statisticsButtonPage2,
                self.settingsButtonPage3,
                self.submitButtonPage3,
                self.statisticsButtonPage3,
                self.clearButtonPage3,
            ]

            # Define cursor and stylesheet for all buttons
            for button in self.buttons:
                button.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
                button.setStyleSheet(
                    """
                    QPushButton {
                        padding: 10px 20px;
                        font-weight: bold;
                        background-color: #DEDDDA;
                        border-radius: 20px;
                    }
                    QPushButton:hover {
                        background-color: #B7B7B6;
                    }
                    QPushButton:pressed {
                        background-color: #DEDDDA;
                    }
                """
                )
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 140))
                shadow_effect.setOffset(3, 3)
                button.setGraphicsEffect(shadow_effect)

            for buttonV in self.buttons_vertical:
                buttonV.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
                buttonV.setStyleSheet(
                    """
                    QPushButton {
                        border-radius: 14px;
                        background-color: #EEEEEE;
                        padding: 10px 20px;
                        font-weight: bold;
                    }
                    QPushButton:hover {
                        background-color: #D7D7D7;
                    }
                    QPushButton:pressed {
                        background-color: #EEEEEE;
                    }
                """
                )
                shadow_effect = QGraphicsDropShadowEffect()
                shadow_effect.setBlurRadius(10)
                shadow_effect.setColor(QColor(0, 0, 0, 110))
                shadow_effect.setOffset(3, 3)
                buttonV.setGraphicsEffect(shadow_effect)

            self.darkModeButton.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
            self.darkModeButton.setStyleSheet(
                """
                QPushButton {
                    padding: 10px 20px;
                    font-weight: bold;
                    background-color: #DEDDDA;
                    border-radius: 14px;
                }
                QPushButton:hover {
                    background-color: #B7B7B6;
                }
                QPushButton:pressed {
                    background-color: #DEDDDA;
                }
            """
            )
            shadow_effect = QGraphicsDropShadowEffect()
            shadow_effect.setBlurRadius(10)
            shadow_effect.setColor(QColor(0, 0, 0, 140))
            shadow_effect.setOffset(3, 3)
            self.darkModeButton.setGraphicsEffect(shadow_effect)

            QtCore.QMetaObject.connectSlotsByName(self)

        except AttributeError as e:
            show_error_dialog(f"An error occurred while setting up the UI: {e}", "Attribute Error")
        except Exception as e:
            show_error_dialog(
                f"An unexpected error occurred: {e}",
                "Unexpected Error",
            )

    @connect_event("homeButton", QtEvents.clicked)
    def show_home_section_click(self):
        self.navigation.show_home_section()

    @connect_event("geneticDataButton", QtEvents.clicked)
    def show_genetic_section_click(self):
        self.navigation.show_genetic_section()

    @connect_event("climaticDataButton", QtEvents.clicked)
    def show_climate_section_click(self):
        self.navigation.show_climate_section()

    @connect_event("resultsButton", QtEvents.clicked)
    def show_results_section_click(self):
        self.navigation.show_results_section()

    @connect_event("helpButton", QtEvents.clicked)
    def open_help_window_click(self):
        self.navigation.open_help_window()

    ################################

    def stop_thread(self):
        self.geneticsPage.genetics.geneticAlignment.stop_worker()
        currThread = self.thread()
        if currThread is QThread:
            if currThread.isRunning():
                currThread.stop()
                currThread.quit()
                currThread.wait()

    def closeEvent(self, event):
        self.stop_thread()
        event.accept()

    @connect_event("darkModeButton", QtEvents.clicked)
    def toggle_dark_mode(self):
        """
        Toggle the application's dark mode setting.

        This method switches the application's theme between dark mode and light mode. It updates the isDarkMode attribute,
        applies the corresponding style, and changes the icon of the darkModeButton.

        Attributes:
            isDarkMode (bool): A flag indicating whether dark mode is currently enabled.

        Actions:
            If dark mode is enabled, apply the dark style and set the darkModeButton icon to the 'light' icon.
            If dark mode is disabled, apply the light style and set the darkModeButton icon to the 'dark' icon.
        """
        try:
            self.isDarkMode = not self.isDarkMode

            style_buttons(self.buttons, self.isDarkMode)
            style_buttons(self.buttons_vertical, self.isDarkMode)

            if self.isDarkMode:
                qtmodern.styles.dark(app)
                self.top_frame.setStyleSheet("background-color: #646464;")
                self.darkModeButton.setIcon(QIcon(":other/light.png"))  # Set the 'light' icon for dark mode
            else:
                qtmodern.styles.light(app)
                self.top_frame.setStyleSheet("background-color: rgb(222, 221, 218);")
                self.darkModeButton.setIcon(QIcon(":other/dark.png"))  # Set the 'dark' icon

            # Common settings for both modes
            self.darkModeButton.setCursor(Qt.CursorShape.PointingHandCursor)
            self.darkModeButton.setStyleSheet(get_button_style(self.isDarkMode))
            self.darkModeButton.setGraphicsEffect(create_shadow_effect(10, 140))

        except AttributeError as e:
            show_error_dialog(f"An error occurred while setting attributes: {e}", "Attribute Error")

        except Exception as e:
            show_error_dialog(f"An unexpected error occurred: {e}")

    ################################################


if __name__ == "__main__":
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--disable-gpu"
    app = QtWidgets.QApplication(sys.argv)

    qtmodern.styles.light(app)

    window = UiMainWindow()

    mw = qtmodern.windows.ModernWindow(window)

    if os.path.exists("scripts/utils/params.yaml") is False:
        shutil.copy("scripts/utils/params_default.yaml", "scripts/utils/params.yaml")

    if primary_screen := app.primaryScreen():
        screen_geometry = primary_screen.availableGeometry()
        center_point = screen_geometry.center()
        x = center_point.x() - mw.width() // 2
        y = center_point.y() - mw.height() // 2

        mw.move(x, y)
        mw.show()

    sys.exit(app.exec())
