import sys
import shutil
import os

os.environ["QT_API"] = "pyqt6"

#import qtmodern
#import qtmodern.styles
#import qtmodern.windows
from aphylogeo.params import Params, os
from event_connector import QtEvents, connect_decorated_methods, connect_event
from navigation import Navigation
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QColor, QIcon
#from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import QGraphicsDropShadowEffect, QVBoxLayout
from Qt import main_ui
from ui_controllers import ClimatePageController, GeneticPageController, ResultPageController
from ui_helpers import create_shadow_effect, get_button_style, style_buttons
from utils import resources_rc  # noqa: F401  # Import the compiled resource module for resolving image resource path
from utils.error_dialog import show_error_dialog

print(QtCore.__file__)