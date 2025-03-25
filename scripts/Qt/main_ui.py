# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'main.ui'
##
## Created by: Qt User Interface Compiler version 6.8.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import (QApplication, QComboBox, QDoubleSpinBox, QFrame,
    QGraphicsView, QGridLayout, QHBoxLayout, QLabel,
    QMainWindow, QPushButton, QSizePolicy, QSpacerItem,
    QSpinBox, QStackedWidget, QStatusBar, QTabWidget,
    QTextBrowser, QTextEdit, QVBoxLayout, QWidget)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1200, 710)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMinimumSize(QSize(1200, 710))
        MainWindow.setMaximumSize(QSize(10000, 10000))
        icon = QIcon()
        icon.addFile(u"../../img/other/sherbrooke.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.verticalLayout = QVBoxLayout(self.centralwidget)
        self.verticalLayout.setSpacing(10)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.top_frame = QFrame(self.centralwidget)
        self.top_frame.setObjectName(u"top_frame")
        self.top_frame.setEnabled(True)
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Maximum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.top_frame.sizePolicy().hasHeightForWidth())
        self.top_frame.setSizePolicy(sizePolicy1)
        palette = QPalette()
        brush = QBrush(QColor(0, 0, 0, 255))
        brush.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.WindowText, brush)
        brush1 = QBrush(QColor(222, 221, 218, 255))
        brush1.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Button, brush1)
        brush2 = QBrush(QColor(255, 255, 255, 255))
        brush2.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Light, brush2)
        brush3 = QBrush(QColor(250, 250, 249, 255))
        brush3.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Midlight, brush3)
        brush4 = QBrush(QColor(123, 122, 122, 255))
        brush4.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Dark, brush4)
        brush5 = QBrush(QColor(164, 163, 163, 255))
        brush5.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Mid, brush5)
        palette.setBrush(QPalette.Active, QPalette.Text, brush)
        palette.setBrush(QPalette.Active, QPalette.BrightText, brush2)
        palette.setBrush(QPalette.Active, QPalette.ButtonText, brush)
        palette.setBrush(QPalette.Active, QPalette.Base, brush1)
        palette.setBrush(QPalette.Active, QPalette.Window, brush1)
        palette.setBrush(QPalette.Active, QPalette.Shadow, brush)
        palette.setBrush(QPalette.Active, QPalette.AlternateBase, brush3)
        brush6 = QBrush(QColor(255, 255, 220, 255))
        brush6.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.ToolTipBase, brush6)
        palette.setBrush(QPalette.Active, QPalette.ToolTipText, brush)
        brush7 = QBrush(QColor(0, 0, 0, 128))
        brush7.setStyle(Qt.NoBrush)
#if QT_VERSION >= QT_VERSION_CHECK(5, 12, 0)
        palette.setBrush(QPalette.Active, QPalette.PlaceholderText, brush7)
#endif
        palette.setBrush(QPalette.Inactive, QPalette.WindowText, brush)
        palette.setBrush(QPalette.Inactive, QPalette.Button, brush1)
        palette.setBrush(QPalette.Inactive, QPalette.Light, brush2)
        palette.setBrush(QPalette.Inactive, QPalette.Midlight, brush3)
        palette.setBrush(QPalette.Inactive, QPalette.Dark, brush4)
        palette.setBrush(QPalette.Inactive, QPalette.Mid, brush5)
        palette.setBrush(QPalette.Inactive, QPalette.Text, brush)
        palette.setBrush(QPalette.Inactive, QPalette.BrightText, brush2)
        palette.setBrush(QPalette.Inactive, QPalette.ButtonText, brush)
        palette.setBrush(QPalette.Inactive, QPalette.Base, brush1)
        palette.setBrush(QPalette.Inactive, QPalette.Window, brush1)
        palette.setBrush(QPalette.Inactive, QPalette.Shadow, brush)
        palette.setBrush(QPalette.Inactive, QPalette.AlternateBase, brush3)
        palette.setBrush(QPalette.Inactive, QPalette.ToolTipBase, brush6)
        palette.setBrush(QPalette.Inactive, QPalette.ToolTipText, brush)
        brush8 = QBrush(QColor(0, 0, 0, 128))
        brush8.setStyle(Qt.NoBrush)
#if QT_VERSION >= QT_VERSION_CHECK(5, 12, 0)
        palette.setBrush(QPalette.Inactive, QPalette.PlaceholderText, brush8)
#endif
        palette.setBrush(QPalette.Disabled, QPalette.WindowText, brush4)
        palette.setBrush(QPalette.Disabled, QPalette.Button, brush1)
        palette.setBrush(QPalette.Disabled, QPalette.Light, brush2)
        palette.setBrush(QPalette.Disabled, QPalette.Midlight, brush3)
        palette.setBrush(QPalette.Disabled, QPalette.Dark, brush4)
        palette.setBrush(QPalette.Disabled, QPalette.Mid, brush5)
        palette.setBrush(QPalette.Disabled, QPalette.Text, brush4)
        palette.setBrush(QPalette.Disabled, QPalette.BrightText, brush2)
        palette.setBrush(QPalette.Disabled, QPalette.ButtonText, brush4)
        palette.setBrush(QPalette.Disabled, QPalette.Base, brush1)
        palette.setBrush(QPalette.Disabled, QPalette.Window, brush1)
        palette.setBrush(QPalette.Disabled, QPalette.Shadow, brush)
        brush9 = QBrush(QColor(246, 245, 244, 255))
        brush9.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Disabled, QPalette.AlternateBase, brush9)
        palette.setBrush(QPalette.Disabled, QPalette.ToolTipBase, brush6)
        palette.setBrush(QPalette.Disabled, QPalette.ToolTipText, brush)
        brush10 = QBrush(QColor(0, 0, 0, 128))
        brush10.setStyle(Qt.NoBrush)
#if QT_VERSION >= QT_VERSION_CHECK(5, 12, 0)
        palette.setBrush(QPalette.Disabled, QPalette.PlaceholderText, brush10)
#endif
        self.top_frame.setPalette(palette)
        self.top_frame.setStyleSheet(u"background-color: rgb(222, 221, 218);")
        self.top_frame.setFrameShape(QFrame.Shape.NoFrame)
        self.top_frame.setFrameShadow(QFrame.Shadow.Sunken)
        self.horizontalLayout = QHBoxLayout(self.top_frame)
        self.horizontalLayout.setSpacing(30)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(9, 15, 9, 15)
        self.appLogo = QLabel(self.top_frame)
        self.appLogo.setObjectName(u"appLogo")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Ignored)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.appLogo.sizePolicy().hasHeightForWidth())
        self.appLogo.setSizePolicy(sizePolicy2)
        self.appLogo.setMaximumSize(QSize(350, 10000))
        self.appLogo.setPixmap(QPixmap(u"../../img/other/Logo_decoupe.png"))
        self.appLogo.setScaledContents(True)
        self.appLogo.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.appLogo.setMargin(0)

        self.horizontalLayout.addWidget(self.appLogo)

        self.homeButton = QPushButton(self.top_frame)
        self.homeButton.setObjectName(u"homeButton")
        sizePolicy3 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.homeButton.sizePolicy().hasHeightForWidth())
        self.homeButton.setSizePolicy(sizePolicy3)
        font = QFont()
        font.setKerning(True)
        self.homeButton.setFont(font)
        icon1 = QIcon()
        icon1.addFile(u":/active/home.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.homeButton.setIcon(icon1)
        self.homeButton.setIconSize(QSize(75, 75))
        self.homeButton.setCheckable(False)
        self.homeButton.setFlat(False)

        self.horizontalLayout.addWidget(self.homeButton)

        self.geneticDataButton = QPushButton(self.top_frame)
        self.geneticDataButton.setObjectName(u"geneticDataButton")
        self.geneticDataButton.setEnabled(True)
        sizePolicy3.setHeightForWidth(self.geneticDataButton.sizePolicy().hasHeightForWidth())
        self.geneticDataButton.setSizePolicy(sizePolicy3)
        self.geneticDataButton.setFont(font)
        icon2 = QIcon()
        icon2.addFile(u":/inactive/genetic.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.geneticDataButton.setIcon(icon2)
        self.geneticDataButton.setIconSize(QSize(75, 75))
        self.geneticDataButton.setCheckable(False)
        self.geneticDataButton.setFlat(False)

        self.horizontalLayout.addWidget(self.geneticDataButton)

        self.climaticDataButton = QPushButton(self.top_frame)
        self.climaticDataButton.setObjectName(u"climaticDataButton")
        self.climaticDataButton.setEnabled(True)
        sizePolicy3.setHeightForWidth(self.climaticDataButton.sizePolicy().hasHeightForWidth())
        self.climaticDataButton.setSizePolicy(sizePolicy3)
        self.climaticDataButton.setFont(font)
        icon3 = QIcon()
        icon3.addFile(u":/inactive/climaticData.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.climaticDataButton.setIcon(icon3)
        self.climaticDataButton.setIconSize(QSize(75, 75))
        self.climaticDataButton.setCheckable(False)
        self.climaticDataButton.setFlat(False)

        self.horizontalLayout.addWidget(self.climaticDataButton)

        self.resultsButton = QPushButton(self.top_frame)
        self.resultsButton.setObjectName(u"resultsButton")
        self.resultsButton.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.resultsButton.sizePolicy().hasHeightForWidth())
        self.resultsButton.setSizePolicy(sizePolicy3)
        icon4 = QIcon()
        icon4.addFile(u":/inactive/result.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        icon4.addFile(u":/active/result.svg", QSize(), QIcon.Mode.Active, QIcon.State.On)
        icon4.addFile(u"\n"
"              :/active/result.svg", QSize(), QIcon.Mode.Selected, QIcon.State.On)
        self.resultsButton.setIcon(icon4)
        self.resultsButton.setIconSize(QSize(75, 75))
        self.resultsButton.setFlat(False)

        self.horizontalLayout.addWidget(self.resultsButton)

        self.helpButton = QPushButton(self.top_frame)
        self.helpButton.setObjectName(u"helpButton")
        sizePolicy3.setHeightForWidth(self.helpButton.sizePolicy().hasHeightForWidth())
        self.helpButton.setSizePolicy(sizePolicy3)
        icon5 = QIcon()
        icon5.addFile(u":/other/help.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.helpButton.setIcon(icon5)
        self.helpButton.setIconSize(QSize(75, 75))
        self.helpButton.setFlat(False)

        self.horizontalLayout.addWidget(self.helpButton)

        self.darkModeButton = QPushButton(self.top_frame)
        self.darkModeButton.setObjectName(u"darkModeButton")
        sizePolicy3.setHeightForWidth(self.darkModeButton.sizePolicy().hasHeightForWidth())
        self.darkModeButton.setSizePolicy(sizePolicy3)
        self.darkModeButton.setStyleSheet(u"QPushButton {\n"
"              color: white; /* White text */\n"
"              border: none;\n"
"              border-radius: 5px; /* Rounded corners */\n"
"              padding: 10px 20px; /* Button padding */\n"
"              font-weight: bold; /* Slightly bold text */\n"
"              box-shadow 0.3s ease, transform 0.3s ease; /* Transitions for smooth effects */\n"
"              }\n"
"\n"
"              QPushButton:hover {\n"
"              background-color: #444; /* Slightly lighter on hover */\n"
"              border-color: #777; /* Even lighter border on hover */\n"
"              box-shadow: 0 0 10px 2px rgba(255, 255, 255, 0.3); /* Subtle glow effect */\n"
"              transform: scale(1.05); /* Slightly increase size on hover */\n"
"              }\n"
"\n"
"              /* Optional: Color Pulsing Animation */\n"
"              @keyframes pulse {\n"
"\n"
"              0% { background-color: #444; }\n"
"\n"
"              100% { background-color: #555; }\n"
"              }\n"
"\n"
"        "
                        "      QPushButton:hover {\n"
"              /* ... other hover styles ... */\n"
"              animation: pulse 1.5s infinite alternate; /* Adjust speed and behavior if needed */\n"
"              }\n"
"            ")
        icon6 = QIcon()
        icon6.addFile(u":/other/dark.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.darkModeButton.setIcon(icon6)
        self.darkModeButton.setIconSize(QSize(40, 40))

        self.horizontalLayout.addWidget(self.darkModeButton)

        self.horizontalLayout.setStretch(0, 1)
        self.horizontalLayout.setStretch(1, 2)
        self.horizontalLayout.setStretch(2, 2)
        self.horizontalLayout.setStretch(3, 2)
        self.horizontalLayout.setStretch(4, 2)
        self.horizontalLayout.setStretch(5, 2)
        self.horizontalLayout.setStretch(6, 1)

        self.verticalLayout.addWidget(self.top_frame)

        self.stackedWidget = QStackedWidget(self.centralwidget)
        self.stackedWidget.setObjectName(u"stackedWidget")
        self.stackedWidget.setEnabled(True)
        sizePolicy.setHeightForWidth(self.stackedWidget.sizePolicy().hasHeightForWidth())
        self.stackedWidget.setSizePolicy(sizePolicy)
        self.HomePage = QWidget()
        self.HomePage.setObjectName(u"HomePage")
        self.verticalLayout_4 = QVBoxLayout(self.HomePage)
        self.verticalLayout_4.setSpacing(0)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.HomeTab = QFrame(self.HomePage)
        self.HomeTab.setObjectName(u"HomeTab")
        self.HomeTab.setFrameShape(QFrame.Shape.NoFrame)
        self.HomeTab.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_5 = QVBoxLayout(self.HomeTab)
        self.verticalLayout_5.setSpacing(0)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.HomeText = QTextEdit(self.HomeTab)
        self.HomeText.setObjectName(u"HomeText")
        sizePolicy.setHeightForWidth(self.HomeText.sizePolicy().hasHeightForWidth())
        self.HomeText.setSizePolicy(sizePolicy)
        font1 = QFont()
        font1.setPointSize(8)
        self.HomeText.setFont(font1)
        self.HomeText.setReadOnly(True)

        self.verticalLayout_5.addWidget(self.HomeText)


        self.verticalLayout_4.addWidget(self.HomeTab)

        self.stackedWidget.addWidget(self.HomePage)
        self.GeneticDataPage = QWidget()
        self.GeneticDataPage.setObjectName(u"GeneticDataPage")
        self.horizontalLayout_5 = QHBoxLayout(self.GeneticDataPage)
        self.horizontalLayout_5.setSpacing(10)
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.GeneticDataSideButtons = QFrame(self.GeneticDataPage)
        self.GeneticDataSideButtons.setObjectName(u"GeneticDataSideButtons")
        sizePolicy4 = QSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Minimum)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.GeneticDataSideButtons.sizePolicy().hasHeightForWidth())
        self.GeneticDataSideButtons.setSizePolicy(sizePolicy4)
        self.GeneticDataSideButtons.setMinimumSize(QSize(200, 0))
        self.GeneticDataSideButtons.setFrameShape(QFrame.Shape.NoFrame)
        self.GeneticDataSideButtons.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_6 = QVBoxLayout(self.GeneticDataSideButtons)
        self.verticalLayout_6.setSpacing(30)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.verticalLayout_6.setContentsMargins(5, 0, 5, 10)
        self.fileBrowserButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.fileBrowserButtonPage1.setObjectName(u"fileBrowserButtonPage1")
        sizePolicy5 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.fileBrowserButtonPage1.sizePolicy().hasHeightForWidth())
        self.fileBrowserButtonPage1.setSizePolicy(sizePolicy5)
        self.fileBrowserButtonPage1.setMinimumSize(QSize(0, 80))
        self.fileBrowserButtonPage1.setFont(font)
        icon7 = QIcon()
        icon7.addFile(u":/other/Browse.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.fileBrowserButtonPage1.setIcon(icon7)
        self.fileBrowserButtonPage1.setIconSize(QSize(60, 90))
        self.fileBrowserButtonPage1.setCheckable(False)
        self.fileBrowserButtonPage1.setFlat(False)

        self.verticalLayout_6.addWidget(self.fileBrowserButtonPage1)

        self.sequenceAlignmentButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.sequenceAlignmentButtonPage1.setObjectName(u"sequenceAlignmentButtonPage1")
        self.sequenceAlignmentButtonPage1.setEnabled(False)
        sizePolicy5.setHeightForWidth(self.sequenceAlignmentButtonPage1.sizePolicy().hasHeightForWidth())
        self.sequenceAlignmentButtonPage1.setSizePolicy(sizePolicy5)
        self.sequenceAlignmentButtonPage1.setMinimumSize(QSize(0, 80))
        icon8 = QIcon()
        icon8.addFile(u":/inactive/sequence.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.sequenceAlignmentButtonPage1.setIcon(icon8)
        self.sequenceAlignmentButtonPage1.setIconSize(QSize(60, 70))
        self.sequenceAlignmentButtonPage1.setFlat(False)

        self.verticalLayout_6.addWidget(self.sequenceAlignmentButtonPage1)

        self.statisticsButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.statisticsButtonPage1.setObjectName(u"statisticsButtonPage1")
        self.statisticsButtonPage1.setEnabled(False)
        sizePolicy5.setHeightForWidth(self.statisticsButtonPage1.sizePolicy().hasHeightForWidth())
        self.statisticsButtonPage1.setSizePolicy(sizePolicy5)
        self.statisticsButtonPage1.setMinimumSize(QSize(0, 80))
        icon9 = QIcon()
        icon9.addFile(u":/inactive/statistics.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.statisticsButtonPage1.setIcon(icon9)
        self.statisticsButtonPage1.setIconSize(QSize(50, 50))
        self.statisticsButtonPage1.setFlat(False)

        self.verticalLayout_6.addWidget(self.statisticsButtonPage1)

        self.geneticTreeButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.geneticTreeButtonPage1.setObjectName(u"geneticTreeButtonPage1")
        self.geneticTreeButtonPage1.setEnabled(False)
        sizePolicy5.setHeightForWidth(self.geneticTreeButtonPage1.sizePolicy().hasHeightForWidth())
        self.geneticTreeButtonPage1.setSizePolicy(sizePolicy5)
        self.geneticTreeButtonPage1.setMinimumSize(QSize(0, 80))
        icon10 = QIcon()
        icon10.addFile(u":/inactive/tree.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.geneticTreeButtonPage1.setIcon(icon10)
        self.geneticTreeButtonPage1.setIconSize(QSize(100, 100))
        self.geneticTreeButtonPage1.setFlat(False)

        self.verticalLayout_6.addWidget(self.geneticTreeButtonPage1)

        self.clearButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.clearButtonPage1.setObjectName(u"clearButtonPage1")
        self.clearButtonPage1.setEnabled(True)
        sizePolicy5.setHeightForWidth(self.clearButtonPage1.sizePolicy().hasHeightForWidth())
        self.clearButtonPage1.setSizePolicy(sizePolicy5)
        self.clearButtonPage1.setMinimumSize(QSize(0, 80))
        self.clearButtonPage1.setFont(font)
        icon11 = QIcon()
        icon11.addFile(u":/other/erase.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.clearButtonPage1.setIcon(icon11)
        self.clearButtonPage1.setIconSize(QSize(60, 90))
        self.clearButtonPage1.setCheckable(False)
        self.clearButtonPage1.setFlat(False)

        self.verticalLayout_6.addWidget(self.clearButtonPage1)

        self.verticalLayout_6.setStretch(0, 1)
        self.verticalLayout_6.setStretch(1, 1)
        self.verticalLayout_6.setStretch(2, 1)
        self.verticalLayout_6.setStretch(3, 1)
        self.verticalLayout_6.setStretch(4, 1)

        self.horizontalLayout_5.addWidget(self.GeneticDataSideButtons)

        self.GeneticDataTabs = QFrame(self.GeneticDataPage)
        self.GeneticDataTabs.setObjectName(u"GeneticDataTabs")
        sizePolicy.setHeightForWidth(self.GeneticDataTabs.sizePolicy().hasHeightForWidth())
        self.GeneticDataTabs.setSizePolicy(sizePolicy)
        self.GeneticDataTabs.setFrameShape(QFrame.Shape.NoFrame)
        self.GeneticDataTabs.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_7 = QVBoxLayout(self.GeneticDataTabs)
        self.verticalLayout_7.setSpacing(0)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.verticalLayout_7.setContentsMargins(0, 0, 6, 0)
        self.tabWidget = QTabWidget(self.GeneticDataTabs)
        self.tabWidget.setObjectName(u"tabWidget")
        self.GenTab1 = QWidget()
        self.GenTab1.setObjectName(u"GenTab1")
        self.verticalLayout_8 = QVBoxLayout(self.GenTab1)
        self.verticalLayout_8.setSpacing(0)
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.verticalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.textEditGenStart = QTextEdit(self.GenTab1)
        self.textEditGenStart.setObjectName(u"textEditGenStart")
        sizePolicy.setHeightForWidth(self.textEditGenStart.sizePolicy().hasHeightForWidth())
        self.textEditGenStart.setSizePolicy(sizePolicy)
        self.textEditGenStart.setFont(font1)
        self.textEditGenStart.setFrameShape(QFrame.Shape.StyledPanel)
        self.textEditGenStart.setReadOnly(True)

        self.verticalLayout_8.addWidget(self.textEditGenStart)

        self.tabWidget.addTab(self.GenTab1, "")
        self.GenTab2 = QWidget()
        self.GenTab2.setObjectName(u"GenTab2")
        self.verticalLayout_9 = QVBoxLayout(self.GenTab2)
        self.verticalLayout_9.setSpacing(0)
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.verticalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.textEditFasta = QTextEdit(self.GenTab2)
        self.textEditFasta.setObjectName(u"textEditFasta")

        self.verticalLayout_9.addWidget(self.textEditFasta)

        self.tabWidget.addTab(self.GenTab2, "")
        self.GenTab3 = QWidget()
        self.GenTab3.setObjectName(u"GenTab3")
        self.verticalLayout_19 = QVBoxLayout(self.GenTab3)
        self.verticalLayout_19.setSpacing(6)
        self.verticalLayout_19.setObjectName(u"verticalLayout_19")
        self.verticalLayout_19.setContentsMargins(0, 0, 0, 0)
        self.horizontalWidget_2 = QWidget(self.GenTab3)
        self.horizontalWidget_2.setObjectName(u"horizontalWidget_2")
        sizePolicy.setHeightForWidth(self.horizontalWidget_2.sizePolicy().hasHeightForWidth())
        self.horizontalWidget_2.setSizePolicy(sizePolicy)
        self.horizontalWidget_2.setMaximumSize(QSize(16777215, 250))
        self.horizontalLayout_6 = QHBoxLayout(self.horizontalWidget_2)
        self.horizontalLayout_6.setSpacing(30)
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.horizontalLayout_6.setContentsMargins(20, 20, 10, 20)
        self.geneticSettingsButton = QPushButton(self.horizontalWidget_2)
        self.geneticSettingsButton.setObjectName(u"geneticSettingsButton")
        sizePolicy6 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)
        sizePolicy6.setHeightForWidth(self.geneticSettingsButton.sizePolicy().hasHeightForWidth())
        self.geneticSettingsButton.setSizePolicy(sizePolicy6)
        self.geneticSettingsButton.setMaximumSize(QSize(16777215, 70))
        self.geneticSettingsButton.setStyleSheet(u"QPushButton {\n"
"                      /* Basic styling */\n"
"                      background-color: #3498db; /* Button color */\n"
"                      color: white; /* Text color */\n"
"                      font-size: 14px; /* Font size */\n"
"                      font-weight: bold; /* Bold text */\n"
"                      border-radius: 8px; /* Rounded corners */\n"
"                      padding: 10px 20px; /* Padding */\n"
"                      border: 2px solid #2980b9; /* Border */\n"
"                      }\n"
"\n"
"                      QPushButton:hover {\n"
"                      /* Hover state styling */\n"
"                      background-color: #2980b9; /* Darker blue on hover */\n"
"                      border-color: #1c5980; /* Darker border on hover */\n"
"                      }\n"
"\n"
"                      QPushButton:pressed {\n"
"                      /* Pressed state styling */\n"
"                      background-color: #1c5980; /* Even darker blue on press */\n"
"         "
                        "             border-color: #145374; /* Darker border on press */\n"
"                      color: #ecf0f1; /* Lighter text on press */\n"
"                      }")

        self.horizontalLayout_6.addWidget(self.geneticSettingsButton)

        self.horizontalSpacer_5 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_6.addItem(self.horizontalSpacer_5)

        self.StartSequenceAlignmentButton = QPushButton(self.horizontalWidget_2)
        self.StartSequenceAlignmentButton.setObjectName(u"StartSequenceAlignmentButton")
        sizePolicy6.setHeightForWidth(self.StartSequenceAlignmentButton.sizePolicy().hasHeightForWidth())
        self.StartSequenceAlignmentButton.setSizePolicy(sizePolicy6)
        self.StartSequenceAlignmentButton.setMaximumSize(QSize(16777215, 70))
        font2 = QFont()
        font2.setBold(True)
        font2.setKerning(True)
        self.StartSequenceAlignmentButton.setFont(font2)
        self.StartSequenceAlignmentButton.setStyleSheet(u"QPushButton {\n"
"                      /* Basic styling */\n"
"                      background-color: #008331; /* Button color */\n"
"                      color: white; /* Text color */\n"
"                      font-size: 14px; /* Font size */\n"
"                      font-weight: bold; /* Bold text */\n"
"                      border-radius: 8px; /* Rounded corners */\n"
"                      padding: 10px 20px; /* Padding */\n"
"                      border: 2px solid #2980b9; /* Border */\n"
"                      }\n"
"\n"
"                      QPushButton:hover {\n"
"                      /* Hover state styling */\n"
"                      background-color: #7a9244; /* Darker blue on hover */\n"
"                      border-color: #1c5980; /* Darker border on hover */\n"
"                      }\n"
"\n"
"                      QPushButton:pressed {\n"
"                      /* Pressed state styling */\n"
"                      background-color: #1c5980; /* Even darker blue on press */\n"
"         "
                        "             border-color: #145374; /* Darker border on press */\n"
"                      color: #ecf0f1; /* Lighter text on press */\n"
"                      }")
        self.StartSequenceAlignmentButton.setIconSize(QSize(60, 90))
        self.StartSequenceAlignmentButton.setCheckable(False)
        self.StartSequenceAlignmentButton.setFlat(False)

        self.horizontalLayout_6.addWidget(self.StartSequenceAlignmentButton)

        self.horizontalSpacer_6 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_6.addItem(self.horizontalSpacer_6)

        self.verticalLayout_17 = QVBoxLayout()
        self.verticalLayout_17.setSpacing(0)
        self.verticalLayout_17.setObjectName(u"verticalLayout_17")
        self.verticalSpacer_12 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_17.addItem(self.verticalSpacer_12)

        self.label_3 = QLabel(self.horizontalWidget_2)
        self.label_3.setObjectName(u"label_3")
        sizePolicy7 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy7)

        self.verticalLayout_17.addWidget(self.label_3)

        self.window_size_spinbox_2 = QSpinBox(self.horizontalWidget_2)
        self.window_size_spinbox_2.setObjectName(u"window_size_spinbox_2")
        self.window_size_spinbox_2.setEnabled(False)
        self.window_size_spinbox_2.setValue(35)

        self.verticalLayout_17.addWidget(self.window_size_spinbox_2)

        self.verticalSpacer_13 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_17.addItem(self.verticalSpacer_13)


        self.horizontalLayout_6.addLayout(self.verticalLayout_17)

        self.verticalLayout_18 = QVBoxLayout()
        self.verticalLayout_18.setSpacing(0)
        self.verticalLayout_18.setObjectName(u"verticalLayout_18")
        self.verticalSpacer_14 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_18.addItem(self.verticalSpacer_14)

        self.label_4 = QLabel(self.horizontalWidget_2)
        self.label_4.setObjectName(u"label_4")
        sizePolicy7.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy7)

        self.verticalLayout_18.addWidget(self.label_4)

        self.starting_position_spinbox_2 = QSpinBox(self.horizontalWidget_2)
        self.starting_position_spinbox_2.setObjectName(u"starting_position_spinbox_2")
        self.starting_position_spinbox_2.setEnabled(False)
        self.starting_position_spinbox_2.setValue(1)

        self.verticalLayout_18.addWidget(self.starting_position_spinbox_2)

        self.verticalSpacer_15 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_18.addItem(self.verticalSpacer_15)


        self.horizontalLayout_6.addLayout(self.verticalLayout_18)

        self.horizontalLayout_6.setStretch(0, 1)
        self.horizontalLayout_6.setStretch(1, 2)
        self.horizontalLayout_6.setStretch(2, 1)
        self.horizontalLayout_6.setStretch(3, 1)
        self.horizontalLayout_6.setStretch(4, 1)
        self.horizontalLayout_6.setStretch(5, 1)

        self.verticalLayout_19.addWidget(self.horizontalWidget_2)

        self.seqAlignLabel = QLabel(self.GenTab3)
        self.seqAlignLabel.setObjectName(u"seqAlignLabel")
        self.seqAlignLabel.setEnabled(True)
        sizePolicy.setHeightForWidth(self.seqAlignLabel.sizePolicy().hasHeightForWidth())
        self.seqAlignLabel.setSizePolicy(sizePolicy)
        self.seqAlignLabel.setMinimumSize(QSize(900, 400))
        self.seqAlignLabel.setMaximumSize(QSize(10000, 10000))
        self.seqAlignLabel.setScaledContents(True)

        self.verticalLayout_19.addWidget(self.seqAlignLabel)

        self.verticalLayout_19.setStretch(0, 1)
        self.verticalLayout_19.setStretch(1, 4)
        self.tabWidget.addTab(self.GenTab3, "")
        self.GenTab4 = QWidget()
        self.GenTab4.setObjectName(u"GenTab4")
        self.verticalLayout_10 = QVBoxLayout(self.GenTab4)
        self.verticalLayout_10.setSpacing(0)
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.verticalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.genStatsHeader = QFrame(self.GenTab4)
        self.genStatsHeader.setObjectName(u"genStatsHeader")
        sizePolicy.setHeightForWidth(self.genStatsHeader.sizePolicy().hasHeightForWidth())
        self.genStatsHeader.setSizePolicy(sizePolicy)
        self.genStatsHeader.setMaximumSize(QSize(10000, 250))
        self.horizontalLayout_7 = QHBoxLayout(self.genStatsHeader)
        self.horizontalLayout_7.setSpacing(20)
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.GenStatsTitle = QLabel(self.genStatsHeader)
        self.GenStatsTitle.setObjectName(u"GenStatsTitle")
        font3 = QFont()
        font3.setPointSize(16)
        self.GenStatsTitle.setFont(font3)
        self.GenStatsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.horizontalLayout_7.addWidget(self.GenStatsTitle)

        self.verticalLayout_22 = QVBoxLayout()
        self.verticalLayout_22.setSpacing(0)
        self.verticalLayout_22.setObjectName(u"verticalLayout_22")
        self.verticalSpacer_10 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_22.addItem(self.verticalSpacer_10)

        self.label_7 = QLabel(self.genStatsHeader)
        self.label_7.setObjectName(u"label_7")
        sizePolicy7.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy7)

        self.verticalLayout_22.addWidget(self.label_7)

        self.referenceComboBox = QComboBox(self.genStatsHeader)
        self.referenceComboBox.setObjectName(u"referenceComboBox")

        self.verticalLayout_22.addWidget(self.referenceComboBox)

        self.verticalSpacer_11 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_22.addItem(self.verticalSpacer_11)


        self.horizontalLayout_7.addLayout(self.verticalLayout_22)

        self.downloadSimilarityButton = QPushButton(self.genStatsHeader)
        self.downloadSimilarityButton.setObjectName(u"downloadSimilarityButton")
        sizePolicy6.setHeightForWidth(self.downloadSimilarityButton.sizePolicy().hasHeightForWidth())
        self.downloadSimilarityButton.setSizePolicy(sizePolicy6)
        self.downloadSimilarityButton.setMaximumSize(QSize(16777215, 50))

        self.horizontalLayout_7.addWidget(self.downloadSimilarityButton)

        self.horizontalLayout_7.setStretch(0, 3)
        self.horizontalLayout_7.setStretch(1, 1)
        self.horizontalLayout_7.setStretch(2, 1)

        self.verticalLayout_10.addWidget(self.genStatsHeader)

        self.textEditGenStats_2 = QWebEngineView(self.GenTab4)
        self.textEditGenStats_2.setObjectName(u"textEditGenStats_2")
        self.textEditGenStats_2.setUrl(QUrl(u"about:blank"))

        self.verticalLayout_10.addWidget(self.textEditGenStats_2)

        self.verticalLayout_10.setStretch(0, 1)
        self.verticalLayout_10.setStretch(1, 4)
        self.tabWidget.addTab(self.GenTab4, "")
        self.GenTab5 = QWidget()
        self.GenTab5.setObjectName(u"GenTab5")
        sizePolicy.setHeightForWidth(self.GenTab5.sizePolicy().hasHeightForWidth())
        self.GenTab5.setSizePolicy(sizePolicy)
        self.verticalLayout_24 = QVBoxLayout(self.GenTab5)
        self.verticalLayout_24.setSpacing(0)
        self.verticalLayout_24.setObjectName(u"verticalLayout_24")
        self.verticalLayout_24.setContentsMargins(0, 0, 0, 0)
        self.horizontalWidget_4 = QWidget(self.GenTab5)
        self.horizontalWidget_4.setObjectName(u"horizontalWidget_4")
        self.horizontalWidget_4.setMaximumSize(QSize(16777215, 250))
        self.horizontalLayout_9 = QHBoxLayout(self.horizontalWidget_4)
        self.horizontalLayout_9.setSpacing(20)
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.horizontalLayout_9.setContentsMargins(-1, -1, 10, -1)
        self.horizontalSpacer_7 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_9.addItem(self.horizontalSpacer_7)

        self.geneticTreescomboBox = QComboBox(self.horizontalWidget_4)
        self.geneticTreescomboBox.setObjectName(u"geneticTreescomboBox")

        self.horizontalLayout_9.addWidget(self.geneticTreescomboBox)

        self.downloadGraphButton = QPushButton(self.horizontalWidget_4)
        self.downloadGraphButton.setObjectName(u"downloadGraphButton")
        sizePolicy6.setHeightForWidth(self.downloadGraphButton.sizePolicy().hasHeightForWidth())
        self.downloadGraphButton.setSizePolicy(sizePolicy6)
        self.downloadGraphButton.setMinimumSize(QSize(0, 35))
        self.downloadGraphButton.setMaximumSize(QSize(16777215, 50))

        self.horizontalLayout_9.addWidget(self.downloadGraphButton)

        self.horizontalLayout_9.setStretch(0, 3)
        self.horizontalLayout_9.setStretch(1, 1)
        self.horizontalLayout_9.setStretch(2, 1)

        self.verticalLayout_24.addWidget(self.horizontalWidget_4)

        self.GeneticTreeLabel = QLabel(self.GenTab5)
        self.GeneticTreeLabel.setObjectName(u"GeneticTreeLabel")
        self.GeneticTreeLabel.setEnabled(True)
        self.GeneticTreeLabel.setMinimumSize(QSize(0, 0))
        self.GeneticTreeLabel.setMaximumSize(QSize(1000000, 100000))
        self.GeneticTreeLabel.setScaledContents(True)

        self.verticalLayout_24.addWidget(self.GeneticTreeLabel)

        self.verticalLayout_24.setStretch(0, 1)
        self.verticalLayout_24.setStretch(1, 4)
        self.tabWidget.addTab(self.GenTab5, "")

        self.verticalLayout_7.addWidget(self.tabWidget)


        self.horizontalLayout_5.addWidget(self.GeneticDataTabs)

        self.horizontalLayout_5.setStretch(0, 1)
        self.horizontalLayout_5.setStretch(1, 5)
        self.stackedWidget.addWidget(self.GeneticDataPage)
        self.ClimaticDataPage = QWidget()
        self.ClimaticDataPage.setObjectName(u"ClimaticDataPage")
        self.horizontalLayout_2 = QHBoxLayout(self.ClimaticDataPage)
        self.horizontalLayout_2.setSpacing(10)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.ClimaticDataSideButtons = QFrame(self.ClimaticDataPage)
        self.ClimaticDataSideButtons.setObjectName(u"ClimaticDataSideButtons")
        sizePolicy4.setHeightForWidth(self.ClimaticDataSideButtons.sizePolicy().hasHeightForWidth())
        self.ClimaticDataSideButtons.setSizePolicy(sizePolicy4)
        self.ClimaticDataSideButtons.setMinimumSize(QSize(200, 0))
        self.ClimaticDataSideButtons.setFrameShape(QFrame.Shape.NoFrame)
        self.ClimaticDataSideButtons.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_2 = QVBoxLayout(self.ClimaticDataSideButtons)
        self.verticalLayout_2.setSpacing(50)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.verticalLayout_2.setContentsMargins(5, 0, 5, 10)
        self.fileBrowserButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.fileBrowserButtonPage2.setObjectName(u"fileBrowserButtonPage2")
        sizePolicy5.setHeightForWidth(self.fileBrowserButtonPage2.sizePolicy().hasHeightForWidth())
        self.fileBrowserButtonPage2.setSizePolicy(sizePolicy5)
        self.fileBrowserButtonPage2.setMinimumSize(QSize(0, 95))
        self.fileBrowserButtonPage2.setFont(font)
        self.fileBrowserButtonPage2.setIcon(icon7)
        self.fileBrowserButtonPage2.setIconSize(QSize(60, 90))
        self.fileBrowserButtonPage2.setCheckable(False)
        self.fileBrowserButtonPage2.setFlat(False)

        self.verticalLayout_2.addWidget(self.fileBrowserButtonPage2)

        self.statisticsButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.statisticsButtonPage2.setObjectName(u"statisticsButtonPage2")
        self.statisticsButtonPage2.setEnabled(False)
        sizePolicy5.setHeightForWidth(self.statisticsButtonPage2.sizePolicy().hasHeightForWidth())
        self.statisticsButtonPage2.setSizePolicy(sizePolicy5)
        self.statisticsButtonPage2.setMinimumSize(QSize(0, 95))
        self.statisticsButtonPage2.setIcon(icon9)
        self.statisticsButtonPage2.setIconSize(QSize(50, 50))
        self.statisticsButtonPage2.setFlat(False)

        self.verticalLayout_2.addWidget(self.statisticsButtonPage2)

        self.climaticTreeButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.climaticTreeButtonPage2.setObjectName(u"climaticTreeButtonPage2")
        self.climaticTreeButtonPage2.setEnabled(False)
        sizePolicy5.setHeightForWidth(self.climaticTreeButtonPage2.sizePolicy().hasHeightForWidth())
        self.climaticTreeButtonPage2.setSizePolicy(sizePolicy5)
        self.climaticTreeButtonPage2.setMinimumSize(QSize(0, 95))
        icon12 = QIcon()
        icon12.addFile(u":/inactive/climatic.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.climaticTreeButtonPage2.setIcon(icon12)
        self.climaticTreeButtonPage2.setIconSize(QSize(70, 60))
        self.climaticTreeButtonPage2.setFlat(False)

        self.verticalLayout_2.addWidget(self.climaticTreeButtonPage2)

        self.clearButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.clearButtonPage2.setObjectName(u"clearButtonPage2")
        sizePolicy5.setHeightForWidth(self.clearButtonPage2.sizePolicy().hasHeightForWidth())
        self.clearButtonPage2.setSizePolicy(sizePolicy5)
        self.clearButtonPage2.setMinimumSize(QSize(0, 95))
        self.clearButtonPage2.setFont(font)
        self.clearButtonPage2.setIcon(icon11)
        self.clearButtonPage2.setIconSize(QSize(60, 90))
        self.clearButtonPage2.setCheckable(False)
        self.clearButtonPage2.setFlat(False)

        self.verticalLayout_2.addWidget(self.clearButtonPage2)

        self.verticalLayout_2.setStretch(0, 1)
        self.verticalLayout_2.setStretch(1, 1)
        self.verticalLayout_2.setStretch(2, 1)
        self.verticalLayout_2.setStretch(3, 1)

        self.horizontalLayout_2.addWidget(self.ClimaticDataSideButtons)

        self.ClimaticDataTabs = QFrame(self.ClimaticDataPage)
        self.ClimaticDataTabs.setObjectName(u"ClimaticDataTabs")
        sizePolicy.setHeightForWidth(self.ClimaticDataTabs.sizePolicy().hasHeightForWidth())
        self.ClimaticDataTabs.setSizePolicy(sizePolicy)
        self.ClimaticDataTabs.setFrameShape(QFrame.Shape.NoFrame)
        self.ClimaticDataTabs.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_3 = QVBoxLayout(self.ClimaticDataTabs)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.tabWidget2 = QTabWidget(self.ClimaticDataTabs)
        self.tabWidget2.setObjectName(u"tabWidget2")
        self.ClimTab1 = QWidget()
        self.ClimTab1.setObjectName(u"ClimTab1")
        self.verticalLayout_21 = QVBoxLayout(self.ClimTab1)
        self.verticalLayout_21.setSpacing(0)
        self.verticalLayout_21.setObjectName(u"verticalLayout_21")
        self.verticalLayout_21.setContentsMargins(0, 0, 0, 0)
        self.climate_data_settings = QWidget(self.ClimTab1)
        self.climate_data_settings.setObjectName(u"climate_data_settings")
        self.climate_data_settings.setAutoFillBackground(True)
        self.horizontalLayout_14 = QHBoxLayout(self.climate_data_settings)
        self.horizontalLayout_14.setObjectName(u"horizontalLayout_14")
        self.horizontalLayout_14.setContentsMargins(10, 0, 0, 0)
        self.label_2 = QLabel(self.climate_data_settings)
        self.label_2.setObjectName(u"label_2")

        self.horizontalLayout_14.addWidget(self.label_2)

        self.max_correlation_climat = QDoubleSpinBox(self.climate_data_settings)
        self.max_correlation_climat.setObjectName(u"max_correlation_climat")
        self.max_correlation_climat.setDecimals(4)
        self.max_correlation_climat.setMaximum(1.000000000000000)
        self.max_correlation_climat.setValue(1.000000000000000)

        self.horizontalLayout_14.addWidget(self.max_correlation_climat)

        self.horizontalSpacer_11 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_14.addItem(self.horizontalSpacer_11)

        self.label = QLabel(self.climate_data_settings)
        self.label.setObjectName(u"label")

        self.horizontalLayout_14.addWidget(self.label)

        self.min_variance_climat = QDoubleSpinBox(self.climate_data_settings)
        self.min_variance_climat.setObjectName(u"min_variance_climat")
        self.min_variance_climat.setDecimals(4)
        self.min_variance_climat.setMaximum(100000.000000000000000)

        self.horizontalLayout_14.addWidget(self.min_variance_climat)

        self.horizontalSpacer_10 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_14.addItem(self.horizontalSpacer_10)

        self.horizontalLayout_14.setStretch(0, 1)
        self.horizontalLayout_14.setStretch(1, 1)
        self.horizontalLayout_14.setStretch(2, 1)
        self.horizontalLayout_14.setStretch(3, 1)
        self.horizontalLayout_14.setStretch(4, 1)
        self.horizontalLayout_14.setStretch(5, 4)

        self.verticalLayout_21.addWidget(self.climate_data_settings)

        self.textEditClimStart = QTextBrowser(self.ClimTab1)
        self.textEditClimStart.setObjectName(u"textEditClimStart")
        self.textEditClimStart.setAutoFillBackground(True)
        self.textEditClimStart.setFrameShape(QFrame.Shape.StyledPanel)
        self.textEditClimStart.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.textEditClimStart.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        self.verticalLayout_21.addWidget(self.textEditClimStart)

        self.verticalLayout_21.setStretch(0, 1)
        self.verticalLayout_21.setStretch(1, 6)
        self.tabWidget2.addTab(self.ClimTab1, "")
        self.ClimTab2 = QWidget()
        self.ClimTab2.setObjectName(u"ClimTab2")
        self.verticalLayout_13 = QVBoxLayout(self.ClimTab2)
        self.verticalLayout_13.setSpacing(0)
        self.verticalLayout_13.setObjectName(u"verticalLayout_13")
        self.verticalLayout_13.setContentsMargins(0, 0, 0, 0)
        self.textEditClimData = QTextBrowser(self.ClimTab2)
        self.textEditClimData.setObjectName(u"textEditClimData")

        self.verticalLayout_13.addWidget(self.textEditClimData)

        self.graphicsViewClimData = QGraphicsView(self.ClimTab2)
        self.graphicsViewClimData.setObjectName(u"graphicsViewClimData")

        self.verticalLayout_13.addWidget(self.graphicsViewClimData)

        self.tabWidget2.addTab(self.ClimTab2, "")
        self.ClimTab4 = QWidget()
        self.ClimTab4.setObjectName(u"ClimTab4")
        self.horizontalLayout_10 = QHBoxLayout(self.ClimTab4)
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.horizontalLayout_10.setContentsMargins(0, 0, 0, 0)
        self.frameClimStats = QFrame(self.ClimTab4)
        self.frameClimStats.setObjectName(u"frameClimStats")
        self.frameClimStats.setFrameShape(QFrame.Shape.NoFrame)
        self.frameClimStats.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_11 = QHBoxLayout(self.frameClimStats)
        self.horizontalLayout_11.setSpacing(4)
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.horizontalLayout_11.setContentsMargins(10, 0, 5, 5)
        self.verticalWidget = QWidget(self.frameClimStats)
        self.verticalWidget.setObjectName(u"verticalWidget")
        sizePolicy8 = QSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Expanding)
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(0)
        sizePolicy8.setHeightForWidth(self.verticalWidget.sizePolicy().hasHeightForWidth())
        self.verticalWidget.setSizePolicy(sizePolicy8)
        self.verticalWidget.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.verticalWidget.setAutoFillBackground(False)
        self.verticalLayout_15 = QVBoxLayout(self.verticalWidget)
        self.verticalLayout_15.setObjectName(u"verticalLayout_15")
        self.StatisticsTitle = QLabel(self.verticalWidget)
        self.StatisticsTitle.setObjectName(u"StatisticsTitle")
        sizePolicy7.setHeightForWidth(self.StatisticsTitle.sizePolicy().hasHeightForWidth())
        self.StatisticsTitle.setSizePolicy(sizePolicy7)
        self.StatisticsTitle.setFont(font3)
        self.StatisticsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.verticalLayout_15.addWidget(self.StatisticsTitle)

        self.verticalSpacer_3 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_15.addItem(self.verticalSpacer_3)

        self.ClimaticChartSettingsTitle = QLabel(self.verticalWidget)
        self.ClimaticChartSettingsTitle.setObjectName(u"ClimaticChartSettingsTitle")
        sizePolicy7.setHeightForWidth(self.ClimaticChartSettingsTitle.sizePolicy().hasHeightForWidth())
        self.ClimaticChartSettingsTitle.setSizePolicy(sizePolicy7)
        self.ClimaticChartSettingsTitle.setFont(font3)
        self.ClimaticChartSettingsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.verticalLayout_15.addWidget(self.ClimaticChartSettingsTitle)

        self.ClimaticChartSettingsTextAxisX = QLabel(self.verticalWidget)
        self.ClimaticChartSettingsTextAxisX.setObjectName(u"ClimaticChartSettingsTextAxisX")
        sizePolicy7.setHeightForWidth(self.ClimaticChartSettingsTextAxisX.sizePolicy().hasHeightForWidth())
        self.ClimaticChartSettingsTextAxisX.setSizePolicy(sizePolicy7)
        self.ClimaticChartSettingsTextAxisX.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.verticalLayout_15.addWidget(self.ClimaticChartSettingsTextAxisX)

        self.ClimaticChartSettingsAxisX = QComboBox(self.verticalWidget)
        self.ClimaticChartSettingsAxisX.setObjectName(u"ClimaticChartSettingsAxisX")
        sizePolicy7.setHeightForWidth(self.ClimaticChartSettingsAxisX.sizePolicy().hasHeightForWidth())
        self.ClimaticChartSettingsAxisX.setSizePolicy(sizePolicy7)
        self.ClimaticChartSettingsAxisX.setMinimumSize(QSize(0, 40))
        self.ClimaticChartSettingsAxisX.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.ClimaticChartSettingsAxisX.setAutoFillBackground(False)
        self.ClimaticChartSettingsAxisX.setEditable(False)
        self.ClimaticChartSettingsAxisX.setInsertPolicy(QComboBox.InsertPolicy.InsertAtBottom)
        self.ClimaticChartSettingsAxisX.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.ClimaticChartSettingsAxisX.setFrame(True)

        self.verticalLayout_15.addWidget(self.ClimaticChartSettingsAxisX)

        self.verticalSpacer_2 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_15.addItem(self.verticalSpacer_2)

        self.ClimaticChartSettingsTextAxisY = QLabel(self.verticalWidget)
        self.ClimaticChartSettingsTextAxisY.setObjectName(u"ClimaticChartSettingsTextAxisY")
        sizePolicy7.setHeightForWidth(self.ClimaticChartSettingsTextAxisY.sizePolicy().hasHeightForWidth())
        self.ClimaticChartSettingsTextAxisY.setSizePolicy(sizePolicy7)
        self.ClimaticChartSettingsTextAxisY.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.verticalLayout_15.addWidget(self.ClimaticChartSettingsTextAxisY)

        self.ClimaticChartSettingsAxisY = QComboBox(self.verticalWidget)
        self.ClimaticChartSettingsAxisY.setObjectName(u"ClimaticChartSettingsAxisY")
        sizePolicy7.setHeightForWidth(self.ClimaticChartSettingsAxisY.sizePolicy().hasHeightForWidth())
        self.ClimaticChartSettingsAxisY.setSizePolicy(sizePolicy7)
        self.ClimaticChartSettingsAxisY.setMinimumSize(QSize(0, 40))
        self.ClimaticChartSettingsAxisY.setEditable(False)
        self.ClimaticChartSettingsAxisY.setFrame(True)

        self.verticalLayout_15.addWidget(self.ClimaticChartSettingsAxisY)

        self.verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_15.addItem(self.verticalSpacer)

        self.PlotsTypesComboBox = QLabel(self.verticalWidget)
        self.PlotsTypesComboBox.setObjectName(u"PlotsTypesComboBox")
        sizePolicy7.setHeightForWidth(self.PlotsTypesComboBox.sizePolicy().hasHeightForWidth())
        self.PlotsTypesComboBox.setSizePolicy(sizePolicy7)
        self.PlotsTypesComboBox.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.verticalLayout_15.addWidget(self.PlotsTypesComboBox)

        self.PlotTypesCombobox = QComboBox(self.verticalWidget)
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.setObjectName(u"PlotTypesCombobox")
        sizePolicy7.setHeightForWidth(self.PlotTypesCombobox.sizePolicy().hasHeightForWidth())
        self.PlotTypesCombobox.setSizePolicy(sizePolicy7)
        self.PlotTypesCombobox.setMinimumSize(QSize(0, 40))
        self.PlotTypesCombobox.setEditable(False)
        self.PlotTypesCombobox.setFrame(True)

        self.verticalLayout_15.addWidget(self.PlotTypesCombobox)

        self.verticalSpacer_4 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_15.addItem(self.verticalSpacer_4)

        self.horizontalLayout_8 = QHBoxLayout()
        self.horizontalLayout_8.setSpacing(0)
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.horizontalSpacer_2 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_8.addItem(self.horizontalSpacer_2)

        self.climatePlotDownloadButton = QPushButton(self.verticalWidget)
        self.climatePlotDownloadButton.setObjectName(u"climatePlotDownloadButton")
        sizePolicy9 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Maximum)
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.climatePlotDownloadButton.sizePolicy().hasHeightForWidth())
        self.climatePlotDownloadButton.setSizePolicy(sizePolicy9)
        self.climatePlotDownloadButton.setMaximumSize(QSize(10000, 16777215))

        self.horizontalLayout_8.addWidget(self.climatePlotDownloadButton)

        self.horizontalSpacer_4 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_8.addItem(self.horizontalSpacer_4)

        self.horizontalLayout_8.setStretch(0, 1)
        self.horizontalLayout_8.setStretch(1, 1)
        self.horizontalLayout_8.setStretch(2, 1)

        self.verticalLayout_15.addLayout(self.horizontalLayout_8)

        self.verticalSpacer_5 = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout_15.addItem(self.verticalSpacer_5)


        self.horizontalLayout_11.addWidget(self.verticalWidget)

        self.climatGraphView = QWebEngineView(self.frameClimStats)
        self.climatGraphView.setObjectName(u"climatGraphView")
        self.climatGraphView.setUrl(QUrl(u"about:blank"))

        self.horizontalLayout_11.addWidget(self.climatGraphView)

        self.horizontalLayout_11.setStretch(0, 2)
        self.horizontalLayout_11.setStretch(1, 4)

        self.horizontalLayout_10.addWidget(self.frameClimStats)

        self.tabWidget2.addTab(self.ClimTab4, "")
        self.ClimTab3 = QWidget()
        self.ClimTab3.setObjectName(u"ClimTab3")
        self.verticalLayout_14 = QVBoxLayout(self.ClimTab3)
        self.verticalLayout_14.setObjectName(u"verticalLayout_14")
        self.verticalLayout_14.setContentsMargins(0, 0, 0, 0)
        self.frameClimTab3 = QFrame(self.ClimTab3)
        self.frameClimTab3.setObjectName(u"frameClimTab3")
        self.frameClimTab3.setFrameShape(QFrame.Shape.NoFrame)
        self.frameClimTab3.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_16 = QVBoxLayout(self.frameClimTab3)
        self.verticalLayout_16.setSpacing(10)
        self.verticalLayout_16.setObjectName(u"verticalLayout_16")
        self.verticalLayout_16.setContentsMargins(10, 15, 10, 5)
        self.horizontalWidget = QWidget(self.frameClimTab3)
        self.horizontalWidget.setObjectName(u"horizontalWidget")
        sizePolicy3.setHeightForWidth(self.horizontalWidget.sizePolicy().hasHeightForWidth())
        self.horizontalWidget.setSizePolicy(sizePolicy3)
        self.horizontalLayout_12 = QHBoxLayout(self.horizontalWidget)
        self.horizontalLayout_12.setObjectName(u"horizontalLayout_12")
        self.horizontalLayout_12.setContentsMargins(0, 0, 0, 0)
        self.preferencesButton = QPushButton(self.horizontalWidget)
        self.preferencesButton.setObjectName(u"preferencesButton")
        sizePolicy3.setHeightForWidth(self.preferencesButton.sizePolicy().hasHeightForWidth())
        self.preferencesButton.setSizePolicy(sizePolicy3)
        self.preferencesButton.setMinimumSize(QSize(0, 50))

        self.horizontalLayout_12.addWidget(self.preferencesButton)

        self.horizontalSpacer_3 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_12.addItem(self.horizontalSpacer_3)

        self.downloadGraphButton2 = QPushButton(self.horizontalWidget)
        self.downloadGraphButton2.setObjectName(u"downloadGraphButton2")
        sizePolicy3.setHeightForWidth(self.downloadGraphButton2.sizePolicy().hasHeightForWidth())
        self.downloadGraphButton2.setSizePolicy(sizePolicy3)
        self.downloadGraphButton2.setMinimumSize(QSize(0, 50))

        self.horizontalLayout_12.addWidget(self.downloadGraphButton2)

        self.climaticTreescomboBox = QComboBox(self.horizontalWidget)
        self.climaticTreescomboBox.setObjectName(u"climaticTreescomboBox")
        sizePolicy3.setHeightForWidth(self.climaticTreescomboBox.sizePolicy().hasHeightForWidth())
        self.climaticTreescomboBox.setSizePolicy(sizePolicy3)

        self.horizontalLayout_12.addWidget(self.climaticTreescomboBox)

        self.horizontalLayout_12.setStretch(0, 1)
        self.horizontalLayout_12.setStretch(1, 2)
        self.horizontalLayout_12.setStretch(2, 1)
        self.horizontalLayout_12.setStretch(3, 1)

        self.verticalLayout_16.addWidget(self.horizontalWidget)

        self.climaticTreesLabel = QLabel(self.frameClimTab3)
        self.climaticTreesLabel.setObjectName(u"climaticTreesLabel")
        sizePolicy.setHeightForWidth(self.climaticTreesLabel.sizePolicy().hasHeightForWidth())
        self.climaticTreesLabel.setSizePolicy(sizePolicy)
        self.climaticTreesLabel.setMinimumSize(QSize(0, 0))
        self.climaticTreesLabel.setMaximumSize(QSize(1000000, 1000000))
        self.climaticTreesLabel.setScaledContents(True)

        self.verticalLayout_16.addWidget(self.climaticTreesLabel)

        self.verticalLayout_16.setStretch(0, 1)
        self.verticalLayout_16.setStretch(1, 4)

        self.verticalLayout_14.addWidget(self.frameClimTab3)

        self.tabWidget2.addTab(self.ClimTab3, "")

        self.verticalLayout_3.addWidget(self.tabWidget2)

        self.verticalLayout_3.setStretch(0, 1)

        self.horizontalLayout_2.addWidget(self.ClimaticDataTabs)

        self.horizontalLayout_2.setStretch(0, 1)
        self.horizontalLayout_2.setStretch(1, 5)
        self.stackedWidget.addWidget(self.ClimaticDataPage)
        self.ResultsPage = QWidget()
        self.ResultsPage.setObjectName(u"ResultsPage")
        self.horizontalLayout_4 = QHBoxLayout(self.ResultsPage)
        self.horizontalLayout_4.setSpacing(10)
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.ResultsSideButtons = QFrame(self.ResultsPage)
        self.ResultsSideButtons.setObjectName(u"ResultsSideButtons")
        sizePolicy4.setHeightForWidth(self.ResultsSideButtons.sizePolicy().hasHeightForWidth())
        self.ResultsSideButtons.setSizePolicy(sizePolicy4)
        self.ResultsSideButtons.setMinimumSize(QSize(200, 0))
        self.ResultsSideButtons.setMaximumSize(QSize(16777215, 16777215))
        self.ResultsSideButtons.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsSideButtons.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_11 = QVBoxLayout(self.ResultsSideButtons)
        self.verticalLayout_11.setSpacing(50)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.verticalLayout_11.setContentsMargins(5, 0, 5, 10)
        self.settingsButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.settingsButtonPage3.setObjectName(u"settingsButtonPage3")
        sizePolicy5.setHeightForWidth(self.settingsButtonPage3.sizePolicy().hasHeightForWidth())
        self.settingsButtonPage3.setSizePolicy(sizePolicy5)
        self.settingsButtonPage3.setMinimumSize(QSize(0, 95))
        icon13 = QIcon()
        icon13.addFile(u":/inactive/settings.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.settingsButtonPage3.setIcon(icon13)
        self.settingsButtonPage3.setIconSize(QSize(50, 50))
        self.settingsButtonPage3.setFlat(False)

        self.verticalLayout_11.addWidget(self.settingsButtonPage3)

        self.submitButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.submitButtonPage3.setObjectName(u"submitButtonPage3")
        sizePolicy5.setHeightForWidth(self.submitButtonPage3.sizePolicy().hasHeightForWidth())
        self.submitButtonPage3.setSizePolicy(sizePolicy5)
        self.submitButtonPage3.setMinimumSize(QSize(0, 95))
        icon14 = QIcon()
        icon14.addFile(u":/inactive/submit.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.submitButtonPage3.setIcon(icon14)
        self.submitButtonPage3.setIconSize(QSize(50, 50))
        self.submitButtonPage3.setFlat(False)

        self.verticalLayout_11.addWidget(self.submitButtonPage3)

        self.statisticsButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.statisticsButtonPage3.setObjectName(u"statisticsButtonPage3")
        sizePolicy5.setHeightForWidth(self.statisticsButtonPage3.sizePolicy().hasHeightForWidth())
        self.statisticsButtonPage3.setSizePolicy(sizePolicy5)
        self.statisticsButtonPage3.setMinimumSize(QSize(0, 95))
        self.statisticsButtonPage3.setIcon(icon9)
        self.statisticsButtonPage3.setIconSize(QSize(50, 50))
        self.statisticsButtonPage3.setFlat(False)

        self.verticalLayout_11.addWidget(self.statisticsButtonPage3)

        self.clearButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.clearButtonPage3.setObjectName(u"clearButtonPage3")
        sizePolicy5.setHeightForWidth(self.clearButtonPage3.sizePolicy().hasHeightForWidth())
        self.clearButtonPage3.setSizePolicy(sizePolicy5)
        self.clearButtonPage3.setMinimumSize(QSize(0, 95))
        self.clearButtonPage3.setFont(font)
        self.clearButtonPage3.setIcon(icon11)
        self.clearButtonPage3.setIconSize(QSize(60, 90))
        self.clearButtonPage3.setCheckable(False)
        self.clearButtonPage3.setFlat(False)

        self.verticalLayout_11.addWidget(self.clearButtonPage3)

        self.verticalLayout_11.setStretch(0, 1)
        self.verticalLayout_11.setStretch(1, 1)
        self.verticalLayout_11.setStretch(2, 1)
        self.verticalLayout_11.setStretch(3, 1)

        self.horizontalLayout_4.addWidget(self.ResultsSideButtons)

        self.ResultsTab = QFrame(self.ResultsPage)
        self.ResultsTab.setObjectName(u"ResultsTab")
        sizePolicy.setHeightForWidth(self.ResultsTab.sizePolicy().hasHeightForWidth())
        self.ResultsTab.setSizePolicy(sizePolicy)
        self.ResultsTab.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsTab.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout_12 = QVBoxLayout(self.ResultsTab)
        self.verticalLayout_12.setObjectName(u"verticalLayout_12")
        self.verticalLayout_12.setContentsMargins(0, 0, 6, 0)
        self.ResultsTitle = QLabel(self.ResultsTab)
        self.ResultsTitle.setObjectName(u"ResultsTitle")
        self.ResultsTitle.setFont(font3)
        self.ResultsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.verticalLayout_12.addWidget(self.ResultsTitle)

        self.tabWidgetResult = QTabWidget(self.ResultsTab)
        self.tabWidgetResult.setObjectName(u"tabWidgetResult")
        self.resultTab = QWidget()
        self.resultTab.setObjectName(u"resultTab")
        sizePolicy.setHeightForWidth(self.resultTab.sizePolicy().hasHeightForWidth())
        self.resultTab.setSizePolicy(sizePolicy)
        self.verticalLayout_23 = QVBoxLayout(self.resultTab)
        self.verticalLayout_23.setSpacing(0)
        self.verticalLayout_23.setObjectName(u"verticalLayout_23")
        self.verticalLayout_23.setContentsMargins(0, 0, 0, 0)
        self.textEditResults = QTextBrowser(self.resultTab)
        self.textEditResults.setObjectName(u"textEditResults")
        self.textEditResults.setMinimumSize(QSize(0, 0))

        self.verticalLayout_23.addWidget(self.textEditResults)

        self.tabWidgetResult.addTab(self.resultTab, "")
        self.resultStatTab = QWidget()
        self.resultStatTab.setObjectName(u"resultStatTab")
        self.verticalLayout_25 = QVBoxLayout(self.resultStatTab)
        self.verticalLayout_25.setSpacing(0)
        self.verticalLayout_25.setObjectName(u"verticalLayout_25")
        self.verticalLayout_25.setContentsMargins(0, 0, 0, 10)
        self.frameResultsStats = QFrame(self.resultStatTab)
        self.frameResultsStats.setObjectName(u"frameResultsStats")
        sizePolicy.setHeightForWidth(self.frameResultsStats.sizePolicy().hasHeightForWidth())
        self.frameResultsStats.setSizePolicy(sizePolicy)
        self.frameResultsStats.setFrameShape(QFrame.Shape.NoFrame)
        self.frameResultsStats.setFrameShadow(QFrame.Shadow.Raised)
        self.horizontalLayout_13 = QHBoxLayout(self.frameResultsStats)
        self.horizontalLayout_13.setSpacing(0)
        self.horizontalLayout_13.setObjectName(u"horizontalLayout_13")
        self.horizontalLayout_13.setContentsMargins(0, 0, 0, 0)
        self.horizontalSpacer_9 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_13.addItem(self.horizontalSpacer_9)

        self.ResultsStatsTitle = QLabel(self.frameResultsStats)
        self.ResultsStatsTitle.setObjectName(u"ResultsStatsTitle")
        self.ResultsStatsTitle.setFont(font3)
        self.ResultsStatsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.horizontalLayout_13.addWidget(self.ResultsStatsTitle)

        self.horizontalSpacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_13.addItem(self.horizontalSpacer)

        self.downloadResultsPlotButton = QPushButton(self.frameResultsStats)
        self.downloadResultsPlotButton.setObjectName(u"downloadResultsPlotButton")
        self.downloadResultsPlotButton.setMinimumSize(QSize(120, 40))

        self.horizontalLayout_13.addWidget(self.downloadResultsPlotButton)

        self.horizontalSpacer_8 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_13.addItem(self.horizontalSpacer_8)

        self.ResultsStatsFilters = QFrame(self.frameResultsStats)
        self.ResultsStatsFilters.setObjectName(u"ResultsStatsFilters")
        self.ResultsStatsFilters.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsStatsFilters.setFrameShadow(QFrame.Shadow.Raised)
        self.gridLayout = QGridLayout(self.ResultsStatsFilters)
        self.gridLayout.setSpacing(0)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout.setContentsMargins(0, 0, 9, 0)
        self.ResultsStatsListTitle = QLabel(self.ResultsStatsFilters)
        self.ResultsStatsListTitle.setObjectName(u"ResultsStatsListTitle")
        self.ResultsStatsListTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.gridLayout.addWidget(self.ResultsStatsListTitle, 0, 0, 1, 1)

        self.criteriaComboBox = QComboBox(self.ResultsStatsFilters)
        self.criteriaComboBox.setObjectName(u"criteriaComboBox")
        sizePolicy10 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        sizePolicy10.setHorizontalStretch(0)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.criteriaComboBox.sizePolicy().hasHeightForWidth())
        self.criteriaComboBox.setSizePolicy(sizePolicy10)
        self.criteriaComboBox.setFrame(True)

        self.gridLayout.addWidget(self.criteriaComboBox, 0, 1, 1, 1)

        self.ResultsStatsListTitle_2 = QLabel(self.ResultsStatsFilters)
        self.ResultsStatsListTitle_2.setObjectName(u"ResultsStatsListTitle_2")
        self.ResultsStatsListTitle_2.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.gridLayout.addWidget(self.ResultsStatsListTitle_2, 1, 0, 1, 1)

        self.phyloTreescomboBox = QComboBox(self.ResultsStatsFilters)
        self.phyloTreescomboBox.setObjectName(u"phyloTreescomboBox")
        sizePolicy10.setHeightForWidth(self.phyloTreescomboBox.sizePolicy().hasHeightForWidth())
        self.phyloTreescomboBox.setSizePolicy(sizePolicy10)
        self.phyloTreescomboBox.setFrame(True)

        self.gridLayout.addWidget(self.phyloTreescomboBox, 1, 1, 1, 1)


        self.horizontalLayout_13.addWidget(self.ResultsStatsFilters)

        self.horizontalLayout_13.setStretch(0, 2)
        self.horizontalLayout_13.setStretch(1, 2)
        self.horizontalLayout_13.setStretch(2, 1)
        self.horizontalLayout_13.setStretch(3, 1)
        self.horizontalLayout_13.setStretch(4, 1)
        self.horizontalLayout_13.setStretch(5, 4)

        self.verticalLayout_25.addWidget(self.frameResultsStats)

        self.PhyloTreeLabel = QLabel(self.resultStatTab)
        self.PhyloTreeLabel.setObjectName(u"PhyloTreeLabel")
        sizePolicy.setHeightForWidth(self.PhyloTreeLabel.sizePolicy().hasHeightForWidth())
        self.PhyloTreeLabel.setSizePolicy(sizePolicy)
        self.PhyloTreeLabel.setMinimumSize(QSize(0, 0))
        self.PhyloTreeLabel.setMaximumSize(QSize(10000, 10000))
        self.PhyloTreeLabel.setScaledContents(True)

        self.verticalLayout_25.addWidget(self.PhyloTreeLabel)

        self.verticalLayout_25.setStretch(0, 1)
        self.verticalLayout_25.setStretch(1, 4)
        self.tabWidgetResult.addTab(self.resultStatTab, "")
        self.resultMapTab = QWidget()
        self.resultMapTab.setObjectName(u"resultMapTab")
        self.verticalLayout_resultMapTab = QVBoxLayout(self.resultMapTab)
        self.verticalLayout_resultMapTab.setSpacing(0)
        self.verticalLayout_resultMapTab.setObjectName(u"verticalLayout_resultMapTab")
        self.verticalLayout_resultMapTab.setContentsMargins(0, 0, 0, 0)
        self.mapImage = QLabel(self.resultMapTab)
        self.mapImage.setObjectName(u"mapImage")
        sizePolicy.setHeightForWidth(self.mapImage.sizePolicy().hasHeightForWidth())
        self.mapImage.setSizePolicy(sizePolicy)
        self.mapImage.setMinimumSize(QSize(0, 0))
        self.mapImage.setMaximumSize(QSize(10000, 10000))
        self.mapImage.setScaledContents(True)

        self.verticalLayout_resultMapTab.addWidget(self.mapImage)

        self.tabWidgetResult.addTab(self.resultMapTab, "")

        self.verticalLayout_12.addWidget(self.tabWidgetResult)

        self.verticalLayout_12.setStretch(0, 1)
        self.verticalLayout_12.setStretch(1, 10)

        self.horizontalLayout_4.addWidget(self.ResultsTab)

        self.horizontalLayout_4.setStretch(0, 1)
        self.horizontalLayout_4.setStretch(1, 5)
        self.stackedWidget.addWidget(self.ResultsPage)

        self.verticalLayout.addWidget(self.stackedWidget)

        self.verticalLayout.setStretch(0, 1)
        self.verticalLayout.setStretch(1, 6)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)

        self.helpButton.setDefault(False)
        self.stackedWidget.setCurrentIndex(2)
        self.tabWidget.setCurrentIndex(3)
        self.tabWidget2.setCurrentIndex(0)
        self.tabWidgetResult.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"iPhyloGeo", None))
        self.appLogo.setText("")
#if QT_CONFIG(tooltip)
        self.homeButton.setToolTip(QCoreApplication.translate("MainWindow", u"Home", None))
#endif // QT_CONFIG(tooltip)
        self.homeButton.setText("")
#if QT_CONFIG(shortcut)
        self.homeButton.setShortcut(QCoreApplication.translate("MainWindow", u"Down", None))
#endif // QT_CONFIG(shortcut)
#if QT_CONFIG(tooltip)
        self.geneticDataButton.setToolTip(QCoreApplication.translate("MainWindow", u"Genetic Data", None))
#endif // QT_CONFIG(tooltip)
        self.geneticDataButton.setText("")
#if QT_CONFIG(shortcut)
        self.geneticDataButton.setShortcut(QCoreApplication.translate("MainWindow", u"Down", None))
#endif // QT_CONFIG(shortcut)
#if QT_CONFIG(tooltip)
        self.climaticDataButton.setToolTip(QCoreApplication.translate("MainWindow", u"Climatic Data", None))
#endif // QT_CONFIG(tooltip)
        self.climaticDataButton.setText("")
        self.resultsButton.setText("")
#if QT_CONFIG(shortcut)
        self.resultsButton.setShortcut(QCoreApplication.translate("MainWindow", u"Ctrl+S", None))
#endif // QT_CONFIG(shortcut)
#if QT_CONFIG(tooltip)
        self.helpButton.setToolTip(QCoreApplication.translate("MainWindow", u"How to use the application", None))
#endif // QT_CONFIG(tooltip)
        self.helpButton.setText("")
        self.darkModeButton.setText("")
        self.HomeText.setMarkdown(QCoreApplication.translate("MainWindow", u"**Welcome to iPhyloGeo**\n"
"\n"
"Thank you for downloading our software.\n"
"\n"
"Here is your guide to using ***iPhyloGeo:***\n"
"\n"
"You can get help by clicking the Help button in the top right corner. \n"
"\n"
"To toggle Dark mode, click the moon-shaped button in the top right corner.\n"
"\n"
"You can process two types of data extractions: **\"Genetic Data\"** and **\"Climatic\n"
"Data\"**.\n"
"\n"
"You can return to the home menu via the **\"Home\"** button.\n"
"\n"
"Once you've selected the method you want to apply to your data, a side menu\n"
"with new buttons will appear, \n"
"\n"
"allowing you to choose and manipulate the files.\n"
"\n"
"Results will be outputted in a CSV file that you can rename and place wherever\n"
"you want.\n"
"\n"
"For more specific information, click the Help button to get details about each\n"
"algorithm.  \n"
"\n"
"**\u00a9 2024 iPhyloGeo. All rights reserved.**\n"
"\n"
"", None))
        self.HomeText.setHtml(QCoreApplication.translate("MainWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:'Segoe UI'; font-size:8pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"center\" style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:20pt; font-weight:700;\">Welcome to iPhyloGeo</span></p>\n"
"<p align=\"center\" style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Thank you for downloading our software.</span></p>\n"
"<p style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-righ"
                        "t:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:16pt;\">Here is your guide to using </span><span style=\" font-size:16pt; font-weight:700; font-style:italic;\">iPhyloGeo:</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">You can get help by clicking the Help button in the top right corner. </span></p>\n"
"<p style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">To toggle Dark mode, click the moon-shaped button in the top right corner.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\""
                        " margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">You can process two types of data extractions: </span><span style=\" font-size:14pt; font-weight:700; color:#0000ff;\">&quot;Genetic Data&quot;</span><span style=\" font-size:14pt;\"> and </span><span style=\" font-size:14pt; font-weight:700; color:#0000ff;\">&quot;Climatic Data&quot;</span><span style=\" font-size:14pt;\">.</span></p>\n"
"<p style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">You can return to the home menu via the </span><span style=\" font-size:14pt; font-weight:700; color:#0000ff;\">&quot;Home&quot;</span><span style=\" font-size:14pt;\"> button.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:5px; margi"
                        "n-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">Once you've selected the method you want to apply to your data, a side menu with new buttons will appear, </span></p>\n"
"<p style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">allowing you to choose and manipulate the files.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">Results will be outputted in a CSV file that you can rename and place wherever you want.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\""
                        "><br /></p>\n"
"<p style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt;\">For more specific information, click the Help button to get details about each algorithm.  </span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p align=\"center\" style=\" margin-top:5px; margin-bottom:5px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt; font-weight:700; color:#535353;\">\u00a9 2024 iPhyloGeo. All rights reserved.</span></p></body></html>", None))
        self.fileBrowserButtonPage1.setText(QCoreApplication.translate("MainWindow", u" File Browser", None))
        self.sequenceAlignmentButtonPage1.setText(QCoreApplication.translate("MainWindow", u" Alignment", None))
        self.statisticsButtonPage1.setText(QCoreApplication.translate("MainWindow", u" Statistics", None))
        self.geneticTreeButtonPage1.setText(QCoreApplication.translate("MainWindow", u"Genetic Tree", None))
        self.clearButtonPage1.setText(QCoreApplication.translate("MainWindow", u" Clear", None))
        self.textEditGenStart.setHtml(QCoreApplication.translate("MainWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><title>Get Started</title><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:'Segoe UI'; font-size:8pt; font-weight:400; font-style:normal;\">\n"
"<h1 align=\"center\" style=\" margin-top:40px; margin-bottom:20px; margin-left:37px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:36pt; font-weight:700; color:#0056b3;\">Get Started!</span></h1>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:37px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">"
                        "For the genetic data extraction, follow these steps to obtain the final results:</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">1. Select a </span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-weight:600; color:#0056b3;\">.fasta </span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">file using the File Browser button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">2. Proceed to the sequence alignment using the Sequence button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','s"
                        "ans-serif'; font-size:16pt;\">3. Create a genetic tree using the Genetic Tree button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">4. Export your results by clicking the Results button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:37px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-style:italic;\">For more specific information, click the </span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-weight:600; font-style:italic;\">Help</span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-style:italic;\"> button to get all the details you need.</span></p>\n"
"<p align=\"center\" style=\" margin-top:40px; margin-bottom:0px; margin-left:37px; margin-right:37px; -qt-block"
                        "-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:11pt; color:#888888;\">\u00a9 2024 iPhyloGeo. All rights reserved.</span></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab1), QCoreApplication.translate("MainWindow", u"Get Started !", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab2), QCoreApplication.translate("MainWindow", u"Fasta File", None))
        self.geneticSettingsButton.setText(QCoreApplication.translate("MainWindow", u"Settings", None))
        self.StartSequenceAlignmentButton.setText(QCoreApplication.translate("MainWindow", u"Start", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"Window size", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"Starting position", None))
        self.seqAlignLabel.setText("")
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab3), QCoreApplication.translate("MainWindow", u"Sequence Alignment", None))
        self.GenStatsTitle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\"\n"
"                        font-size:24pt; text-decoration: underline;\">Alignment\n"
"                        Chart</span></p></body></html>", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"Reference", None))
        self.downloadSimilarityButton.setText(QCoreApplication.translate("MainWindow", u"Download", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab4), QCoreApplication.translate("MainWindow", u"Species Stats", None))
        self.downloadGraphButton.setText(QCoreApplication.translate("MainWindow", u"Download", None))
        self.GeneticTreeLabel.setText("")
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab5), QCoreApplication.translate("MainWindow", u"Genetic Tree", None))
        self.fileBrowserButtonPage2.setText(QCoreApplication.translate("MainWindow", u" File Browser", None))
        self.statisticsButtonPage2.setText(QCoreApplication.translate("MainWindow", u" Statistics", None))
        self.climaticTreeButtonPage2.setText(QCoreApplication.translate("MainWindow", u" Climatic Tree", None))
        self.clearButtonPage2.setText(QCoreApplication.translate("MainWindow", u" Clear", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"Max Correlation", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"Min Variance", None))
        self.textEditClimStart.setHtml(QCoreApplication.translate("MainWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><title>Get Started</title><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:'Segoe UI'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<h1 align=\"center\" style=\" margin-top:40px; margin-bottom:20px; margin-left:37px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:36pt; font-weight:600; color:#0056b3;\">Get Started!</span></h1>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:37px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">"
                        "For the climatic data extraction, follow these steps to obtain the final results:</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">1. Select a </span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-weight:600; color:#0056b3;\">.csv</span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\"> file using the File Browser button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">2. Proceed to the climatic tree creation using the Climatic Tree button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'A"
                        "rial','sans-serif'; font-size:16pt;\">3. Generate the statistics using the Statistics button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:57px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt;\">4. Export and visualize your results by clicking the Results button.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:20px; margin-left:37px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-style:italic;\">For more specific information, click the </span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-weight:600; font-style:italic;\">Help</span><span style=\" font-family:'Arial','sans-serif'; font-size:16pt; font-style:italic;\"> button to get all the details you need.</span><span style=\" font-family:'Arial','sans-serif'; font-size:11pt;\"> </span><span style=\" font-fam"
                        "ily:'MS Shell Dlg 2'; font-size:11pt;\">      </span><span style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt;\">                </span></p>\n"
"<p align=\"center\" style=\" margin-top:40px; margin-bottom:0px; margin-left:37px; margin-right:37px; -qt-block-indent:0; text-indent:0px; line-height:160%;\"><span style=\" font-family:'Arial','sans-serif'; font-size:11pt; color:#888888;\">\u00a9 2024 iPhyloGeo. All rights reserved.</span></p></body></html>", None))
        self.tabWidget2.setTabText(self.tabWidget2.indexOf(self.ClimTab1), QCoreApplication.translate("MainWindow", u"Get Started !", None))
        self.tabWidget2.setTabText(self.tabWidget2.indexOf(self.ClimTab2), QCoreApplication.translate("MainWindow", u"Climatic Data", None))
        self.StatisticsTitle.setText(QCoreApplication.translate("MainWindow", u"\n"
"                        <html><head/><body><p><span style=\"\n"
"                        font-size:24pt; text-decoration:\n"
"                        underline;\">Statistics</span></p></body></html>", None))
#if QT_CONFIG(whatsthis)
        self.ClimaticChartSettingsTitle.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.ClimaticChartSettingsTitle.setText(QCoreApplication.translate("MainWindow", u"Generate your Graph", None))
#if QT_CONFIG(whatsthis)
        self.ClimaticChartSettingsTextAxisX.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.ClimaticChartSettingsTextAxisX.setText(QCoreApplication.translate("MainWindow", u"Insert X axis data", None))
#if QT_CONFIG(whatsthis)
        self.ClimaticChartSettingsTextAxisY.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.ClimaticChartSettingsTextAxisY.setText(QCoreApplication.translate("MainWindow", u"Insert Y axis data", None))
#if QT_CONFIG(whatsthis)
        self.PlotsTypesComboBox.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.PlotsTypesComboBox.setText(QCoreApplication.translate("MainWindow", u"Choose Plot Type", None))
        self.PlotTypesCombobox.setItemText(0, QCoreApplication.translate("MainWindow", u"Scatter Plot", None))
        self.PlotTypesCombobox.setItemText(1, QCoreApplication.translate("MainWindow", u"Line Plot", None))
        self.PlotTypesCombobox.setItemText(2, QCoreApplication.translate("MainWindow", u"Bar Graph", None))
        self.PlotTypesCombobox.setItemText(3, QCoreApplication.translate("MainWindow", u"Violin Plot", None))
        self.PlotTypesCombobox.setItemText(4, QCoreApplication.translate("MainWindow", u"Pie Plot", None))
        self.PlotTypesCombobox.setItemText(5, QCoreApplication.translate("MainWindow", u"Correlation", None))

        self.climatePlotDownloadButton.setText(QCoreApplication.translate("MainWindow", u"Download", None))
        self.tabWidget2.setTabText(self.tabWidget2.indexOf(self.ClimTab4), QCoreApplication.translate("MainWindow", u"Statistics", None))
        self.preferencesButton.setText(QCoreApplication.translate("MainWindow", u"Preferences", None))
        self.downloadGraphButton2.setText(QCoreApplication.translate("MainWindow", u"Download Graph", None))
        self.climaticTreesLabel.setText("")
        self.tabWidget2.setTabText(self.tabWidget2.indexOf(self.ClimTab3), QCoreApplication.translate("MainWindow", u"Climatic Tree", None))
        self.settingsButtonPage3.setText(QCoreApplication.translate("MainWindow", u" Settings", None))
        self.submitButtonPage3.setText(QCoreApplication.translate("MainWindow", u"Save", None))
        self.statisticsButtonPage3.setText(QCoreApplication.translate("MainWindow", u" Statistics", None))
        self.clearButtonPage3.setText(QCoreApplication.translate("MainWindow", u" Clear", None))
        self.ResultsTitle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\"\n"
"                  font-size:20pt;\n"
"                  font-weight:600;\">Results</span></p></body></html>", None))
        self.tabWidgetResult.setTabText(self.tabWidgetResult.indexOf(self.resultTab), QCoreApplication.translate("MainWindow", u"Results", None))
        self.ResultsStatsTitle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\"\n"
"                    font-size:24pt; text-decoration:\n"
"                    underline;\">Statistics</span></p></body></html>", None))
        self.downloadResultsPlotButton.setText(QCoreApplication.translate("MainWindow", u"Download", None))
#if QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle.setText(QCoreApplication.translate("MainWindow", u"Condition", None))
#if QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle_2.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle_2.setText(QCoreApplication.translate("MainWindow", u"trees", None))
        self.PhyloTreeLabel.setText("")
        self.tabWidgetResult.setTabText(self.tabWidgetResult.indexOf(self.resultStatTab), QCoreApplication.translate("MainWindow", u"Statistics", None))
        self.mapImage.setText("")
        self.tabWidgetResult.setTabText(self.tabWidgetResult.indexOf(self.resultMapTab), QCoreApplication.translate("MainWindow", u"Map", None))
    # retranslateUi

