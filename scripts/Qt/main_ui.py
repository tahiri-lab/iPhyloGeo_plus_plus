# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'main.ui'
##
## Created by: Qt User Interface Compiler version 6.8.1
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
from PySide6.QtWidgets import (QApplication, QComboBox, QFrame, QGraphicsView,
    QLabel, QMainWindow, QPushButton, QSizePolicy,
    QSpinBox, QStackedWidget, QStatusBar, QTabWidget,
    QTextBrowser, QTextEdit, QWidget)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1143, 670)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMinimumSize(QSize(1143, 670))
        MainWindow.setMaximumSize(QSize(1143, 670))
        icon = QIcon()
        icon.addFile(u"../../../../../../hazem/.designer/img/sherbrooke.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.top_frame = QFrame(self.centralwidget)
        self.top_frame.setObjectName(u"top_frame")
        self.top_frame.setGeometry(QRect(0, -1, 1151, 101))
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
        self.helpButton = QPushButton(self.top_frame)
        self.helpButton.setObjectName(u"helpButton")
        self.helpButton.setGeometry(QRect(970, 10, 101, 81))
        icon1 = QIcon()
        icon1.addFile(u":/other/help.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.helpButton.setIcon(icon1)
        self.helpButton.setIconSize(QSize(80, 80))
        self.helpButton.setFlat(False)
        self.climaticDataButton = QPushButton(self.top_frame)
        self.climaticDataButton.setObjectName(u"climaticDataButton")
        self.climaticDataButton.setEnabled(True)
        self.climaticDataButton.setGeometry(QRect(630, 10, 101, 81))
        font = QFont()
        font.setKerning(True)
        self.climaticDataButton.setFont(font)
        icon2 = QIcon()
        icon2.addFile(u":/inactive/climaticData.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.climaticDataButton.setIcon(icon2)
        self.climaticDataButton.setIconSize(QSize(75, 75))
        self.climaticDataButton.setCheckable(False)
        self.climaticDataButton.setFlat(False)
        self.appLogo = QLabel(self.top_frame)
        self.appLogo.setObjectName(u"appLogo")
        self.appLogo.setGeometry(QRect(0, 10, 271, 81))
        self.appLogo.setPixmap(QPixmap(u"../../img/other/Logo_decoupe.png"))
        self.appLogo.setScaledContents(True)
        self.appLogo.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.geneticDataButton = QPushButton(self.top_frame)
        self.geneticDataButton.setObjectName(u"geneticDataButton")
        self.geneticDataButton.setEnabled(True)
        self.geneticDataButton.setGeometry(QRect(460, 10, 101, 81))
        self.geneticDataButton.setFont(font)
        icon3 = QIcon()
        icon3.addFile(u":/inactive/genetic.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.geneticDataButton.setIcon(icon3)
        self.geneticDataButton.setIconSize(QSize(75, 75))
        self.geneticDataButton.setCheckable(False)
        self.geneticDataButton.setFlat(False)
        self.darkModeButton = QPushButton(self.top_frame)
        self.darkModeButton.setObjectName(u"darkModeButton")
        self.darkModeButton.setGeometry(QRect(1090, 10, 41, 41))
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
        icon4 = QIcon()
        icon4.addFile(u":/other/dark.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.darkModeButton.setIcon(icon4)
        self.darkModeButton.setIconSize(QSize(40, 40))
        self.homeButton = QPushButton(self.top_frame)
        self.homeButton.setObjectName(u"homeButton")
        self.homeButton.setGeometry(QRect(290, 10, 101, 81))
        self.homeButton.setFont(font)
        icon5 = QIcon()
        icon5.addFile(u":/active/home.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.homeButton.setIcon(icon5)
        self.homeButton.setIconSize(QSize(70, 70))
        self.homeButton.setCheckable(False)
        self.homeButton.setFlat(False)
        self.resultsButton = QPushButton(self.top_frame)
        self.resultsButton.setObjectName(u"resultsButton")
        self.resultsButton.setEnabled(False)
        self.resultsButton.setGeometry(QRect(800, 10, 101, 81))
        icon6 = QIcon()
        icon6.addFile(u":/inactive/result.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        icon6.addFile(u":/active/result.svg", QSize(), QIcon.Mode.Active, QIcon.State.On)
        icon6.addFile(u"\n"
"              :/active/result.svg", QSize(), QIcon.Mode.Selected, QIcon.State.On)
        self.resultsButton.setIcon(icon6)
        self.resultsButton.setIconSize(QSize(60, 60))
        self.resultsButton.setFlat(False)
        self.stackedWidget = QStackedWidget(self.centralwidget)
        self.stackedWidget.setObjectName(u"stackedWidget")
        self.stackedWidget.setEnabled(True)
        self.stackedWidget.setGeometry(QRect(0, 99, 1151, 581))
        sizePolicy.setHeightForWidth(self.stackedWidget.sizePolicy().hasHeightForWidth())
        self.stackedWidget.setSizePolicy(sizePolicy)
        self.HomePage = QWidget()
        self.HomePage.setObjectName(u"HomePage")
        self.HomeTab = QFrame(self.HomePage)
        self.HomeTab.setObjectName(u"HomeTab")
        self.HomeTab.setGeometry(QRect(0, 0, 1141, 551))
        self.HomeTab.setFrameShape(QFrame.Shape.NoFrame)
        self.HomeTab.setFrameShadow(QFrame.Shadow.Raised)
        self.HomeText = QTextEdit(self.HomeTab)
        self.HomeText.setObjectName(u"HomeText")
        self.HomeText.setGeometry(QRect(10, 30, 1121, 531))
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.HomeText.sizePolicy().hasHeightForWidth())
        self.HomeText.setSizePolicy(sizePolicy1)
        font1 = QFont()
        font1.setPointSize(8)
        self.HomeText.setFont(font1)
        self.HomeText.setReadOnly(True)
        self.stackedWidget.addWidget(self.HomePage)
        self.GeneticDataPage = QWidget()
        self.GeneticDataPage.setObjectName(u"GeneticDataPage")
        self.GeneticDataSideButtons = QFrame(self.GeneticDataPage)
        self.GeneticDataSideButtons.setObjectName(u"GeneticDataSideButtons")
        self.GeneticDataSideButtons.setGeometry(QRect(10, 0, 191, 541))
        self.GeneticDataSideButtons.setFrameShape(QFrame.Shape.NoFrame)
        self.GeneticDataSideButtons.setFrameShadow(QFrame.Shadow.Raised)
        self.sequenceAlignmentButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.sequenceAlignmentButtonPage1.setObjectName(u"sequenceAlignmentButtonPage1")
        self.sequenceAlignmentButtonPage1.setEnabled(False)
        self.sequenceAlignmentButtonPage1.setGeometry(QRect(0, 120, 181, 81))
        icon7 = QIcon()
        icon7.addFile(u":/inactive/sequence.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.sequenceAlignmentButtonPage1.setIcon(icon7)
        self.sequenceAlignmentButtonPage1.setIconSize(QSize(60, 70))
        self.sequenceAlignmentButtonPage1.setFlat(False)
        self.statisticsButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.statisticsButtonPage1.setObjectName(u"statisticsButtonPage1")
        self.statisticsButtonPage1.setEnabled(False)
        self.statisticsButtonPage1.setGeometry(QRect(0, 230, 181, 81))
        icon8 = QIcon()
        icon8.addFile(u":/inactive/statistics.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.statisticsButtonPage1.setIcon(icon8)
        self.statisticsButtonPage1.setIconSize(QSize(50, 50))
        self.statisticsButtonPage1.setFlat(False)
        self.fileBrowserButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.fileBrowserButtonPage1.setObjectName(u"fileBrowserButtonPage1")
        self.fileBrowserButtonPage1.setGeometry(QRect(0, 10, 181, 81))
        self.fileBrowserButtonPage1.setFont(font)
        icon9 = QIcon()
        icon9.addFile(u":/other/Browse.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.fileBrowserButtonPage1.setIcon(icon9)
        self.fileBrowserButtonPage1.setIconSize(QSize(60, 90))
        self.fileBrowserButtonPage1.setCheckable(False)
        self.fileBrowserButtonPage1.setFlat(False)
        self.clearButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.clearButtonPage1.setObjectName(u"clearButtonPage1")
        self.clearButtonPage1.setEnabled(True)
        self.clearButtonPage1.setGeometry(QRect(0, 450, 181, 81))
        self.clearButtonPage1.setFont(font)
        icon10 = QIcon()
        icon10.addFile(u":/other/erase.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.clearButtonPage1.setIcon(icon10)
        self.clearButtonPage1.setIconSize(QSize(60, 90))
        self.clearButtonPage1.setCheckable(False)
        self.clearButtonPage1.setFlat(False)
        self.geneticTreeButtonPage1 = QPushButton(self.GeneticDataSideButtons)
        self.geneticTreeButtonPage1.setObjectName(u"geneticTreeButtonPage1")
        self.geneticTreeButtonPage1.setEnabled(False)
        self.geneticTreeButtonPage1.setGeometry(QRect(0, 340, 181, 81))
        icon11 = QIcon()
        icon11.addFile(u":/inactive/tree.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.geneticTreeButtonPage1.setIcon(icon11)
        self.geneticTreeButtonPage1.setIconSize(QSize(100, 100))
        self.geneticTreeButtonPage1.setFlat(False)
        self.GeneticDataTabs = QFrame(self.GeneticDataPage)
        self.GeneticDataTabs.setObjectName(u"GeneticDataTabs")
        self.GeneticDataTabs.setGeometry(QRect(200, 0, 941, 541))
        sizePolicy1.setHeightForWidth(self.GeneticDataTabs.sizePolicy().hasHeightForWidth())
        self.GeneticDataTabs.setSizePolicy(sizePolicy1)
        self.GeneticDataTabs.setFrameShape(QFrame.Shape.NoFrame)
        self.GeneticDataTabs.setFrameShadow(QFrame.Shadow.Raised)
        self.tabWidget = QTabWidget(self.GeneticDataTabs)
        self.tabWidget.setObjectName(u"tabWidget")
        self.tabWidget.setGeometry(QRect(0, 0, 921, 541))
        self.GenTab1 = QWidget()
        self.GenTab1.setObjectName(u"GenTab1")
        self.textEditGenStart = QTextEdit(self.GenTab1)
        self.textEditGenStart.setObjectName(u"textEditGenStart")
        self.textEditGenStart.setGeometry(QRect(0, 0, 926, 511))
        sizePolicy1.setHeightForWidth(self.textEditGenStart.sizePolicy().hasHeightForWidth())
        self.textEditGenStart.setSizePolicy(sizePolicy1)
        self.textEditGenStart.setFont(font1)
        self.textEditGenStart.setReadOnly(True)
        self.tabWidget.addTab(self.GenTab1, "")
        self.GenTab2 = QWidget()
        self.GenTab2.setObjectName(u"GenTab2")
        self.textEditFasta = QTextEdit(self.GenTab2)
        self.textEditFasta.setObjectName(u"textEditFasta")
        self.textEditFasta.setGeometry(QRect(0, 0, 926, 511))
        self.tabWidget.addTab(self.GenTab2, "")
        self.GenTab3 = QWidget()
        self.GenTab3.setObjectName(u"GenTab3")
        self.seqAlignLabel = QLabel(self.GenTab3)
        self.seqAlignLabel.setObjectName(u"seqAlignLabel")
        self.seqAlignLabel.setEnabled(True)
        self.seqAlignLabel.setGeometry(QRect(10, 90, 900, 400))
        self.seqAlignLabel.setMinimumSize(QSize(900, 400))
        self.seqAlignLabel.setMaximumSize(QSize(900, 400))
        self.GenStatsFilters = QFrame(self.GenTab3)
        self.GenStatsFilters.setObjectName(u"GenStatsFilters")
        self.GenStatsFilters.setGeometry(QRect(440, 10, 491, 61))
        self.GenStatsFilters.setFrameShape(QFrame.Shape.NoFrame)
        self.GenStatsFilters.setFrameShadow(QFrame.Shadow.Raised)
        self.starting_position_spinbox_2 = QSpinBox(self.GenStatsFilters)
        self.starting_position_spinbox_2.setObjectName(u"starting_position_spinbox_2")
        self.starting_position_spinbox_2.setEnabled(False)
        self.starting_position_spinbox_2.setGeometry(QRect(350, 30, 111, 22))
        self.starting_position_spinbox_2.setValue(1)
        self.window_size_spinbox_2 = QSpinBox(self.GenStatsFilters)
        self.window_size_spinbox_2.setObjectName(u"window_size_spinbox_2")
        self.window_size_spinbox_2.setEnabled(False)
        self.window_size_spinbox_2.setGeometry(QRect(180, 30, 111, 22))
        self.window_size_spinbox_2.setValue(35)
        self.label_3 = QLabel(self.GenStatsFilters)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setGeometry(QRect(190, 10, 81, 16))
        self.label_4 = QLabel(self.GenStatsFilters)
        self.label_4.setObjectName(u"label_4")
        self.label_4.setGeometry(QRect(360, 10, 121, 16))
        self.StartSequenceAlignmentButton = QPushButton(self.GenTab3)
        self.StartSequenceAlignmentButton.setObjectName(u"StartSequenceAlignmentButton")
        self.StartSequenceAlignmentButton.setGeometry(QRect(380, 20, 111, 51))
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
        self.geneticSettingsButton = QPushButton(self.GenTab3)
        self.geneticSettingsButton.setObjectName(u"geneticSettingsButton")
        self.geneticSettingsButton.setGeometry(QRect(30, 17, 111, 51))
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
        self.tabWidget.addTab(self.GenTab3, "")
        self.GenTab4 = QWidget()
        self.GenTab4.setObjectName(u"GenTab4")
        self.GenStats = QFrame(self.GenTab4)
        self.GenStats.setObjectName(u"GenStats")
        self.GenStats.setGeometry(QRect(0, 0, 921, 511))
        self.GenStats.setFrameShape(QFrame.Shape.NoFrame)
        self.GenStats.setFrameShadow(QFrame.Shadow.Raised)
        self.GenStatsFilters1 = QFrame(self.GenStats)
        self.GenStatsFilters1.setObjectName(u"GenStatsFilters1")
        self.GenStatsFilters1.setGeometry(QRect(440, 20, 481, 61))
        self.GenStatsFilters1.setFrameShape(QFrame.Shape.NoFrame)
        self.GenStatsFilters1.setFrameShadow(QFrame.Shadow.Raised)
        self.downloadSimilarityButton = QPushButton(self.GenStatsFilters1)
        self.downloadSimilarityButton.setObjectName(u"downloadSimilarityButton")
        self.downloadSimilarityButton.setGeometry(QRect(360, 10, 93, 28))
        self.startingPositionSimilaritySpinBox = QSpinBox(self.GenStatsFilters1)
        self.startingPositionSimilaritySpinBox.setObjectName(u"startingPositionSimilaritySpinBox")
        self.startingPositionSimilaritySpinBox.setGeometry(QRect(70, 20, 91, 22))
        self.startingPositionSimilaritySpinBox.setMaximum(9999999)
        self.label_6 = QLabel(self.GenStatsFilters1)
        self.label_6.setObjectName(u"label_6")
        self.label_6.setGeometry(QRect(70, 0, 111, 16))
        self.referenceComboBox = QComboBox(self.GenStatsFilters1)
        self.referenceComboBox.setObjectName(u"referenceComboBox")
        self.referenceComboBox.setGeometry(QRect(210, 20, 73, 22))
        self.label_7 = QLabel(self.GenStatsFilters1)
        self.label_7.setObjectName(u"label_7")
        self.label_7.setGeometry(QRect(210, 0, 111, 16))
        self.GenStatsTitle = QLabel(self.GenStats)
        self.GenStatsTitle.setObjectName(u"GenStatsTitle")
        self.GenStatsTitle.setGeometry(QRect(-100, 0, 541, 81))
        font3 = QFont()
        font3.setPointSize(16)
        self.GenStatsTitle.setFont(font3)
        self.GenStatsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.textEditGenStats_2 = QLabel(self.GenStats)
        self.textEditGenStats_2.setObjectName(u"textEditGenStats_2")
        self.textEditGenStats_2.setGeometry(QRect(10, 100, 900, 400))
        self.textEditGenStats_2.setMinimumSize(QSize(900, 400))
        self.textEditGenStats_2.setMaximumSize(QSize(800, 400))
        self.similarityWindowSizeSpinBox = QSpinBox(self.GenStats)
        self.similarityWindowSizeSpinBox.setObjectName(u"similarityWindowSizeSpinBox")
        self.similarityWindowSizeSpinBox.setGeometry(QRect(370, 40, 91, 22))
        self.similarityWindowSizeSpinBox.setMinimum(1)
        self.similarityWindowSizeSpinBox.setMaximum(1000)
        self.similarityWindowSizeSpinBox.setValue(100)
        self.label_5 = QLabel(self.GenStats)
        self.label_5.setObjectName(u"label_5")
        self.label_5.setGeometry(QRect(380, 20, 81, 16))
        self.tabWidget.addTab(self.GenTab4, "")
        self.GenTab5 = QWidget()
        self.GenTab5.setObjectName(u"GenTab5")
        self.GeneticTreeLabel = QLabel(self.GenTab5)
        self.GeneticTreeLabel.setObjectName(u"GeneticTreeLabel")
        self.GeneticTreeLabel.setEnabled(True)
        self.GeneticTreeLabel.setGeometry(QRect(0, 60, 921, 450))
        self.GeneticTreeLabel.setMinimumSize(QSize(921, 450))
        self.GeneticTreeLabel.setMaximumSize(QSize(921, 450))
        self.geneticTreescomboBox = QComboBox(self.GenTab5)
        self.geneticTreescomboBox.setObjectName(u"geneticTreescomboBox")
        self.geneticTreescomboBox.setGeometry(QRect(740, 10, 131, 22))
        self.downloadGraphButton = QPushButton(self.GenTab5)
        self.downloadGraphButton.setObjectName(u"downloadGraphButton")
        self.downloadGraphButton.setGeometry(QRect(610, 0, 121, 41))
        self.tabWidget.addTab(self.GenTab5, "")
        self.stackedWidget.addWidget(self.GeneticDataPage)
        self.ClimaticDataPage = QWidget()
        self.ClimaticDataPage.setObjectName(u"ClimaticDataPage")
        self.ClimaticDataSideButtons = QFrame(self.ClimaticDataPage)
        self.ClimaticDataSideButtons.setObjectName(u"ClimaticDataSideButtons")
        self.ClimaticDataSideButtons.setGeometry(QRect(10, 0, 191, 541))
        self.ClimaticDataSideButtons.setFrameShape(QFrame.Shape.NoFrame)
        self.ClimaticDataSideButtons.setFrameShadow(QFrame.Shadow.Raised)
        self.clearButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.clearButtonPage2.setObjectName(u"clearButtonPage2")
        self.clearButtonPage2.setGeometry(QRect(0, 440, 181, 81))
        self.clearButtonPage2.setFont(font)
        self.clearButtonPage2.setIcon(icon10)
        self.clearButtonPage2.setIconSize(QSize(60, 90))
        self.clearButtonPage2.setCheckable(False)
        self.clearButtonPage2.setFlat(False)
        self.climaticTreeButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.climaticTreeButtonPage2.setObjectName(u"climaticTreeButtonPage2")
        self.climaticTreeButtonPage2.setEnabled(False)
        self.climaticTreeButtonPage2.setGeometry(QRect(0, 300, 181, 81))
        icon12 = QIcon()
        icon12.addFile(u":/inactive/climatic.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.climaticTreeButtonPage2.setIcon(icon12)
        self.climaticTreeButtonPage2.setIconSize(QSize(70, 60))
        self.climaticTreeButtonPage2.setFlat(False)
        self.fileBrowserButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.fileBrowserButtonPage2.setObjectName(u"fileBrowserButtonPage2")
        self.fileBrowserButtonPage2.setGeometry(QRect(0, 20, 181, 81))
        self.fileBrowserButtonPage2.setFont(font)
        self.fileBrowserButtonPage2.setIcon(icon9)
        self.fileBrowserButtonPage2.setIconSize(QSize(60, 90))
        self.fileBrowserButtonPage2.setCheckable(False)
        self.fileBrowserButtonPage2.setFlat(False)
        self.statisticsButtonPage2 = QPushButton(self.ClimaticDataSideButtons)
        self.statisticsButtonPage2.setObjectName(u"statisticsButtonPage2")
        self.statisticsButtonPage2.setEnabled(False)
        self.statisticsButtonPage2.setGeometry(QRect(0, 160, 181, 81))
        self.statisticsButtonPage2.setIcon(icon8)
        self.statisticsButtonPage2.setIconSize(QSize(50, 50))
        self.statisticsButtonPage2.setFlat(False)
        self.ClimaticDataTabs = QFrame(self.ClimaticDataPage)
        self.ClimaticDataTabs.setObjectName(u"ClimaticDataTabs")
        self.ClimaticDataTabs.setGeometry(QRect(200, 0, 941, 541))
        self.ClimaticDataTabs.setFrameShape(QFrame.Shape.NoFrame)
        self.ClimaticDataTabs.setFrameShadow(QFrame.Shadow.Raised)
        self.tabWidget2 = QTabWidget(self.ClimaticDataTabs)
        self.tabWidget2.setObjectName(u"tabWidget2")
        self.tabWidget2.setGeometry(QRect(0, 0, 931, 541))
        self.ClimTab1 = QWidget()
        self.ClimTab1.setObjectName(u"ClimTab1")
        self.textEditClimStart = QTextBrowser(self.ClimTab1)
        self.textEditClimStart.setObjectName(u"textEditClimStart")
        self.textEditClimStart.setGeometry(QRect(0, 0, 931, 511))
        self.tabWidget2.addTab(self.ClimTab1, "")
        self.ClimTab2 = QWidget()
        self.ClimTab2.setObjectName(u"ClimTab2")
        self.textEditClimData = QTextBrowser(self.ClimTab2)
        self.textEditClimData.setObjectName(u"textEditClimData")
        self.textEditClimData.setGeometry(QRect(0, 0, 931, 231))
        self.graphicsViewClimData = QGraphicsView(self.ClimTab2)
        self.graphicsViewClimData.setObjectName(u"graphicsViewClimData")
        self.graphicsViewClimData.setGeometry(QRect(0, 260, 931, 251))
        self.tabWidget2.addTab(self.ClimTab2, "")
        self.ClimTab4 = QWidget()
        self.ClimTab4.setObjectName(u"ClimTab4")
        self.frameClimStats = QFrame(self.ClimTab4)
        self.frameClimStats.setObjectName(u"frameClimStats")
        self.frameClimStats.setGeometry(QRect(0, 0, 921, 551))
        self.frameClimStats.setFrameShape(QFrame.Shape.NoFrame)
        self.frameClimStats.setFrameShadow(QFrame.Shadow.Raised)
        self.StatisticsTitle = QLabel(self.frameClimStats)
        self.StatisticsTitle.setObjectName(u"StatisticsTitle")
        self.StatisticsTitle.setGeometry(QRect(0, 0, 541, 81))
        self.StatisticsTitle.setFont(font3)
        self.StatisticsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.ClimaticChartSettings = QLabel(self.frameClimStats)
        self.ClimaticChartSettings.setObjectName(u"ClimaticChartSettings")
        self.ClimaticChartSettings.setGeometry(QRect(30, 90, 340, 391))
        self.ClimaticChartSettings.setMinimumSize(QSize(340, 391))
        self.ClimaticChartSettings.setMaximumSize(QSize(340, 391))
        self.ClimaticChart_2 = QLabel(self.frameClimStats)
        self.ClimaticChart_2.setObjectName(u"ClimaticChart_2")
        self.ClimaticChart_2.setGeometry(QRect(390, 10, 520, 500))
        self.ClimaticChart_2.setMinimumSize(QSize(520, 500))
        self.ClimaticChart_2.setMaximumSize(QSize(520, 500))
        self.ClimaticChartSettingsTitle = QLabel(self.frameClimStats)
        self.ClimaticChartSettingsTitle.setObjectName(u"ClimaticChartSettingsTitle")
        self.ClimaticChartSettingsTitle.setGeometry(QRect(50, 90, 301, 51))
        self.ClimaticChartSettingsTitle.setFont(font3)
        self.ClimaticChartSettingsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.ClimaticChartSettingsTextAxisX = QLabel(self.frameClimStats)
        self.ClimaticChartSettingsTextAxisX.setObjectName(u"ClimaticChartSettingsTextAxisX")
        self.ClimaticChartSettingsTextAxisX.setGeometry(QRect(120, 150, 161, 31))
        self.ClimaticChartSettingsTextAxisX.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.ClimaticChartSettingsAxisX = QComboBox(self.frameClimStats)
        self.ClimaticChartSettingsAxisX.setObjectName(u"ClimaticChartSettingsAxisX")
        self.ClimaticChartSettingsAxisX.setGeometry(QRect(40, 180, 321, 31))
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.ClimaticChartSettingsAxisX.sizePolicy().hasHeightForWidth())
        self.ClimaticChartSettingsAxisX.setSizePolicy(sizePolicy2)
        self.ClimaticChartSettingsAxisX.setLayoutDirection(Qt.LayoutDirection.LeftToRight)
        self.ClimaticChartSettingsAxisX.setAutoFillBackground(False)
        self.ClimaticChartSettingsAxisX.setEditable(False)
        self.ClimaticChartSettingsAxisX.setInsertPolicy(QComboBox.InsertPolicy.InsertAtBottom)
        self.ClimaticChartSettingsAxisX.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.ClimaticChartSettingsAxisX.setFrame(True)
        self.ClimaticChartSettingsTextAxisY = QLabel(self.frameClimStats)
        self.ClimaticChartSettingsTextAxisY.setObjectName(u"ClimaticChartSettingsTextAxisY")
        self.ClimaticChartSettingsTextAxisY.setGeometry(QRect(120, 230, 161, 31))
        self.ClimaticChartSettingsTextAxisY.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.ClimaticChartSettingsAxisY = QComboBox(self.frameClimStats)
        self.ClimaticChartSettingsAxisY.setObjectName(u"ClimaticChartSettingsAxisY")
        self.ClimaticChartSettingsAxisY.setGeometry(QRect(40, 260, 321, 30))
        sizePolicy2.setHeightForWidth(self.ClimaticChartSettingsAxisY.sizePolicy().hasHeightForWidth())
        self.ClimaticChartSettingsAxisY.setSizePolicy(sizePolicy2)
        self.ClimaticChartSettingsAxisY.setEditable(False)
        self.ClimaticChartSettingsAxisY.setFrame(True)
        self.PlotTypesCombobox = QComboBox(self.frameClimStats)
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.addItem("")
        self.PlotTypesCombobox.setObjectName(u"PlotTypesCombobox")
        self.PlotTypesCombobox.setGeometry(QRect(40, 340, 321, 30))
        sizePolicy2.setHeightForWidth(self.PlotTypesCombobox.sizePolicy().hasHeightForWidth())
        self.PlotTypesCombobox.setSizePolicy(sizePolicy2)
        self.PlotTypesCombobox.setEditable(False)
        self.PlotTypesCombobox.setFrame(True)
        self.PlotsTypesComboBox = QLabel(self.frameClimStats)
        self.PlotsTypesComboBox.setObjectName(u"PlotsTypesComboBox")
        self.PlotsTypesComboBox.setGeometry(QRect(120, 310, 161, 31))
        self.PlotsTypesComboBox.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.climatePlotDownloadButton = QPushButton(self.frameClimStats)
        self.climatePlotDownloadButton.setObjectName(u"climatePlotDownloadButton")
        self.climatePlotDownloadButton.setGeometry(QRect(150, 410, 93, 28))
        self.tabWidget2.addTab(self.ClimTab4, "")
        self.ClimTab3 = QWidget()
        self.ClimTab3.setObjectName(u"ClimTab3")
        self.frameClimTab3 = QFrame(self.ClimTab3)
        self.frameClimTab3.setObjectName(u"frameClimTab3")
        self.frameClimTab3.setGeometry(QRect(0, 0, 931, 511))
        self.frameClimTab3.setFrameShape(QFrame.Shape.NoFrame)
        self.frameClimTab3.setFrameShadow(QFrame.Shadow.Raised)
        self.climaticTreesLabel = QLabel(self.frameClimTab3)
        self.climaticTreesLabel.setObjectName(u"climaticTreesLabel")
        self.climaticTreesLabel.setGeometry(QRect(4, 65, 911, 441))
        self.climaticTreesLabel.setMinimumSize(QSize(911, 441))
        self.climaticTreesLabel.setMaximumSize(QSize(911, 441))
        self.climaticTreescomboBox = QComboBox(self.frameClimTab3)
        self.climaticTreescomboBox.setObjectName(u"climaticTreescomboBox")
        self.climaticTreescomboBox.setGeometry(QRect(790, 20, 131, 22))
        self.downloadGraphButton2 = QPushButton(self.frameClimTab3)
        self.downloadGraphButton2.setObjectName(u"downloadGraphButton2")
        self.downloadGraphButton2.setGeometry(QRect(650, 10, 121, 41))
        self.preferencesButton = QPushButton(self.frameClimTab3)
        self.preferencesButton.setObjectName(u"preferencesButton")
        self.preferencesButton.setGeometry(QRect(10, 10, 121, 41))
        self.tabWidget2.addTab(self.ClimTab3, "")
        self.stackedWidget.addWidget(self.ClimaticDataPage)
        self.ResultsPage = QWidget()
        self.ResultsPage.setObjectName(u"ResultsPage")
        self.ResultsSideButtons = QFrame(self.ResultsPage)
        self.ResultsSideButtons.setObjectName(u"ResultsSideButtons")
        self.ResultsSideButtons.setGeometry(QRect(10, 0, 191, 541))
        self.ResultsSideButtons.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsSideButtons.setFrameShadow(QFrame.Shadow.Raised)
        self.settingsButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.settingsButtonPage3.setObjectName(u"settingsButtonPage3")
        self.settingsButtonPage3.setGeometry(QRect(0, 20, 181, 81))
        icon13 = QIcon()
        icon13.addFile(u":/inactive/settings.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.settingsButtonPage3.setIcon(icon13)
        self.settingsButtonPage3.setIconSize(QSize(50, 50))
        self.settingsButtonPage3.setFlat(False)
        self.submitButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.submitButtonPage3.setObjectName(u"submitButtonPage3")
        self.submitButtonPage3.setGeometry(QRect(0, 160, 181, 81))
        icon14 = QIcon()
        icon14.addFile(u":/inactive/submit.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.submitButtonPage3.setIcon(icon14)
        self.submitButtonPage3.setIconSize(QSize(50, 50))
        self.submitButtonPage3.setFlat(False)
        self.statisticsButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.statisticsButtonPage3.setObjectName(u"statisticsButtonPage3")
        self.statisticsButtonPage3.setGeometry(QRect(0, 300, 181, 81))
        self.statisticsButtonPage3.setIcon(icon8)
        self.statisticsButtonPage3.setIconSize(QSize(50, 50))
        self.statisticsButtonPage3.setFlat(False)
        self.clearButtonPage3 = QPushButton(self.ResultsSideButtons)
        self.clearButtonPage3.setObjectName(u"clearButtonPage3")
        self.clearButtonPage3.setGeometry(QRect(0, 440, 181, 81))
        self.clearButtonPage3.setFont(font)
        self.clearButtonPage3.setIcon(icon10)
        self.clearButtonPage3.setIconSize(QSize(60, 90))
        self.clearButtonPage3.setCheckable(False)
        self.clearButtonPage3.setFlat(False)
        self.ResultsTab = QFrame(self.ResultsPage)
        self.ResultsTab.setObjectName(u"ResultsTab")
        self.ResultsTab.setGeometry(QRect(200, 0, 931, 541))
        self.ResultsTab.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsTab.setFrameShadow(QFrame.Shadow.Raised)
        self.ResultsTitle = QLabel(self.ResultsTab)
        self.ResultsTitle.setObjectName(u"ResultsTitle")
        self.ResultsTitle.setGeometry(QRect(0, 0, 931, 61))
        self.ResultsTitle.setFont(font3)
        self.ResultsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.tabWidgetResult = QTabWidget(self.ResultsTab)
        self.tabWidgetResult.setObjectName(u"tabWidgetResult")
        self.tabWidgetResult.setGeometry(QRect(0, 30, 931, 461))
        self.resultTab = QWidget()
        self.resultTab.setObjectName(u"resultTab")
        self.textEditResults = QTextBrowser(self.resultTab)
        self.textEditResults.setObjectName(u"textEditResults")
        self.textEditResults.setGeometry(QRect(0, 0, 931, 461))
        self.tabWidgetResult.addTab(self.resultTab, "")
        self.resultStatTab = QWidget()
        self.resultStatTab.setObjectName(u"resultStatTab")
        self.frameResultsStats = QFrame(self.resultStatTab)
        self.frameResultsStats.setObjectName(u"frameResultsStats")
        self.frameResultsStats.setGeometry(QRect(0, 0, 921, 541))
        self.frameResultsStats.setFrameShape(QFrame.Shape.NoFrame)
        self.frameResultsStats.setFrameShadow(QFrame.Shadow.Raised)
        self.ResultsStatsFilters = QFrame(self.frameResultsStats)
        self.ResultsStatsFilters.setObjectName(u"ResultsStatsFilters")
        self.ResultsStatsFilters.setGeometry(QRect(550, 0, 371, 81))
        self.ResultsStatsFilters.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsStatsFilters.setFrameShadow(QFrame.Shadow.Raised)
        self.criteriaComboBox = QComboBox(self.ResultsStatsFilters)
        self.criteriaComboBox.setObjectName(u"criteriaComboBox")
        self.criteriaComboBox.setGeometry(QRect(160, 10, 201, 30))
        sizePolicy2.setHeightForWidth(self.criteriaComboBox.sizePolicy().hasHeightForWidth())
        self.criteriaComboBox.setSizePolicy(sizePolicy2)
        self.criteriaComboBox.setFrame(True)
        self.ResultsStatsListTitle = QLabel(self.ResultsStatsFilters)
        self.ResultsStatsListTitle.setObjectName(u"ResultsStatsListTitle")
        self.ResultsStatsListTitle.setGeometry(QRect(0, 10, 161, 31))
        self.ResultsStatsListTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.phyloTreescomboBox = QComboBox(self.ResultsStatsFilters)
        self.phyloTreescomboBox.setObjectName(u"phyloTreescomboBox")
        self.phyloTreescomboBox.setGeometry(QRect(160, 50, 201, 30))
        sizePolicy2.setHeightForWidth(self.phyloTreescomboBox.sizePolicy().hasHeightForWidth())
        self.phyloTreescomboBox.setSizePolicy(sizePolicy2)
        self.phyloTreescomboBox.setFrame(True)
        self.ResultsStatsListTitle_2 = QLabel(self.ResultsStatsFilters)
        self.ResultsStatsListTitle_2.setObjectName(u"ResultsStatsListTitle_2")
        self.ResultsStatsListTitle_2.setGeometry(QRect(0, 50, 161, 31))
        self.ResultsStatsListTitle_2.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.ResultsStatsTitle = QLabel(self.frameResultsStats)
        self.ResultsStatsTitle.setObjectName(u"ResultsStatsTitle")
        self.ResultsStatsTitle.setGeometry(QRect(0, 0, 541, 81))
        self.ResultsStatsTitle.setFont(font3)
        self.ResultsStatsTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.PhyloTreeLabel = QLabel(self.frameResultsStats)
        self.PhyloTreeLabel.setObjectName(u"PhyloTreeLabel")
        self.PhyloTreeLabel.setGeometry(QRect(0, 0, 921, 450))
        self.PhyloTreeLabel.setMinimumSize(QSize(921, 450))
        self.PhyloTreeLabel.setMaximumSize(QSize(921, 450))
        self.downloadResultsPlotButton = QPushButton(self.frameResultsStats)
        self.downloadResultsPlotButton.setObjectName(u"downloadResultsPlotButton")
        self.downloadResultsPlotButton.setGeometry(QRect(420, 30, 93, 28))
        self.PhyloTreeLabel.raise_()
        self.ResultsStatsFilters.raise_()
        self.ResultsStatsTitle.raise_()
        self.downloadResultsPlotButton.raise_()
        self.tabWidgetResult.addTab(self.resultStatTab, "")
        self.stackedWidget.addWidget(self.ResultsPage)
        self.ResultsStats = QWidget()
        self.ResultsStats.setObjectName(u"ResultsStats")
        self.ResultsStatsButtons = QFrame(self.ResultsStats)
        self.ResultsStatsButtons.setObjectName(u"ResultsStatsButtons")
        self.ResultsStatsButtons.setGeometry(QRect(0, 0, 201, 551))
        self.ResultsStatsButtons.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsStatsButtons.setFrameShadow(QFrame.Shadow.Raised)
        self.settingsButtonPage4 = QPushButton(self.ResultsStatsButtons)
        self.settingsButtonPage4.setObjectName(u"settingsButtonPage4")
        self.settingsButtonPage4.setGeometry(QRect(10, 20, 181, 81))
        self.settingsButtonPage4.setIcon(icon13)
        self.settingsButtonPage4.setIconSize(QSize(50, 50))
        self.settingsButtonPage4.setFlat(False)
        self.submitButtonPage4 = QPushButton(self.ResultsStatsButtons)
        self.submitButtonPage4.setObjectName(u"submitButtonPage4")
        self.submitButtonPage4.setGeometry(QRect(10, 160, 181, 81))
        self.submitButtonPage4.setIcon(icon14)
        self.submitButtonPage4.setIconSize(QSize(50, 50))
        self.submitButtonPage4.setFlat(False)
        self.statisticsButtonPage4 = QPushButton(self.ResultsStatsButtons)
        self.statisticsButtonPage4.setObjectName(u"statisticsButtonPage4")
        self.statisticsButtonPage4.setGeometry(QRect(10, 300, 181, 81))
        self.statisticsButtonPage4.setIcon(icon8)
        self.statisticsButtonPage4.setIconSize(QSize(50, 50))
        self.statisticsButtonPage4.setFlat(False)
        self.clearButtonPage4 = QPushButton(self.ResultsStatsButtons)
        self.clearButtonPage4.setObjectName(u"clearButtonPage4")
        self.clearButtonPage4.setGeometry(QRect(10, 440, 181, 81))
        self.clearButtonPage4.setFont(font)
        self.clearButtonPage4.setIcon(icon10)
        self.clearButtonPage4.setIconSize(QSize(60, 90))
        self.clearButtonPage4.setCheckable(False)
        self.clearButtonPage4.setFlat(False)
        self.ResultsStatsTab = QFrame(self.ResultsStats)
        self.ResultsStatsTab.setObjectName(u"ResultsStatsTab")
        self.ResultsStatsTab.setGeometry(QRect(200, -20, 921, 571))
        self.ResultsStatsTab.setFrameShape(QFrame.Shape.NoFrame)
        self.ResultsStatsTab.setFrameShadow(QFrame.Shadow.Raised)
        self.stackedWidget.addWidget(self.ResultsStats)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)

        self.helpButton.setDefault(False)
        self.stackedWidget.setCurrentIndex(1)
        self.tabWidget.setCurrentIndex(2)
        self.tabWidget2.setCurrentIndex(3)
        self.tabWidgetResult.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"iPhyloGeo", None))
#if QT_CONFIG(tooltip)
        self.helpButton.setToolTip(QCoreApplication.translate("MainWindow", u"How to use the application", None))
#endif // QT_CONFIG(tooltip)
        self.helpButton.setText("")
#if QT_CONFIG(tooltip)
        self.climaticDataButton.setToolTip(QCoreApplication.translate("MainWindow", u"Climatic Data", None))
#endif // QT_CONFIG(tooltip)
        self.climaticDataButton.setText("")
        self.appLogo.setText("")
#if QT_CONFIG(tooltip)
        self.geneticDataButton.setToolTip(QCoreApplication.translate("MainWindow", u"Genetic Data", None))
#endif // QT_CONFIG(tooltip)
        self.geneticDataButton.setText("")
#if QT_CONFIG(shortcut)
        self.geneticDataButton.setShortcut(QCoreApplication.translate("MainWindow", u"Down", None))
#endif // QT_CONFIG(shortcut)
        self.darkModeButton.setText("")
#if QT_CONFIG(tooltip)
        self.homeButton.setToolTip(QCoreApplication.translate("MainWindow", u"Home", None))
#endif // QT_CONFIG(tooltip)
        self.homeButton.setText("")
#if QT_CONFIG(shortcut)
        self.homeButton.setShortcut(QCoreApplication.translate("MainWindow", u"Down", None))
#endif // QT_CONFIG(shortcut)
        self.resultsButton.setText("")
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
        self.sequenceAlignmentButtonPage1.setText(QCoreApplication.translate("MainWindow", u" Alignment", None))
        self.statisticsButtonPage1.setText(QCoreApplication.translate("MainWindow", u" Statistics", None))
        self.fileBrowserButtonPage1.setText(QCoreApplication.translate("MainWindow", u" File Browser", None))
        self.clearButtonPage1.setText(QCoreApplication.translate("MainWindow", u" Clear", None))
        self.geneticTreeButtonPage1.setText(QCoreApplication.translate("MainWindow", u"Genetic Tree", None))
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
        self.seqAlignLabel.setText("")
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"Window size", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"Starting position", None))
        self.StartSequenceAlignmentButton.setText(QCoreApplication.translate("MainWindow", u"Start", None))
        self.geneticSettingsButton.setText(QCoreApplication.translate("MainWindow", u"Settings", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab3), QCoreApplication.translate("MainWindow", u"Sequence Alignment", None))
        self.downloadSimilarityButton.setText(QCoreApplication.translate("MainWindow", u"Download", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u"Starting Position", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"Reference", None))
        self.GenStatsTitle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\"\n"
"                        font-size:24pt; text-decoration: underline;\">Alignment\n"
"                        Chart</span></p></body></html>", None))
        self.textEditGenStats_2.setText("")
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"Window size", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab4), QCoreApplication.translate("MainWindow", u"Species Stats", None))
        self.GeneticTreeLabel.setText("")
        self.downloadGraphButton.setText(QCoreApplication.translate("MainWindow", u"Download Graph", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GenTab5), QCoreApplication.translate("MainWindow", u"Genetic Tree", None))
        self.clearButtonPage2.setText(QCoreApplication.translate("MainWindow", u" Clear", None))
        self.climaticTreeButtonPage2.setText(QCoreApplication.translate("MainWindow", u" Climatic Tree", None))
        self.fileBrowserButtonPage2.setText(QCoreApplication.translate("MainWindow", u" File Browser", None))
        self.statisticsButtonPage2.setText(QCoreApplication.translate("MainWindow", u" Statistics", None))
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
        self.ClimaticChartSettings.setText("")
        self.ClimaticChart_2.setText("")
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
        self.PlotTypesCombobox.setItemText(0, QCoreApplication.translate("MainWindow", u"Scatter Plot", None))
        self.PlotTypesCombobox.setItemText(1, QCoreApplication.translate("MainWindow", u"Line Plot", None))
        self.PlotTypesCombobox.setItemText(2, QCoreApplication.translate("MainWindow", u"Bar Graph", None))
        self.PlotTypesCombobox.setItemText(3, QCoreApplication.translate("MainWindow", u"Violin Plot", None))
        self.PlotTypesCombobox.setItemText(4, QCoreApplication.translate("MainWindow", u"Pie Plot", None))

#if QT_CONFIG(whatsthis)
        self.PlotsTypesComboBox.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.PlotsTypesComboBox.setText(QCoreApplication.translate("MainWindow", u"choose Plot Type", None))
        self.climatePlotDownloadButton.setText(QCoreApplication.translate("MainWindow", u"Download", None))
        self.tabWidget2.setTabText(self.tabWidget2.indexOf(self.ClimTab4), QCoreApplication.translate("MainWindow", u"Statistics", None))
        self.climaticTreesLabel.setText("")
        self.downloadGraphButton2.setText(QCoreApplication.translate("MainWindow", u"Download Graph", None))
        self.preferencesButton.setText(QCoreApplication.translate("MainWindow", u"Preferences", None))
        self.tabWidget2.setTabText(self.tabWidget2.indexOf(self.ClimTab3), QCoreApplication.translate("MainWindow", u"Climatic Tree", None))
        self.settingsButtonPage3.setText(QCoreApplication.translate("MainWindow", u" Settings", None))
        self.submitButtonPage3.setText(QCoreApplication.translate("MainWindow", u"Save", None))
        self.statisticsButtonPage3.setText(QCoreApplication.translate("MainWindow", u" Statistics", None))
        self.clearButtonPage3.setText(QCoreApplication.translate("MainWindow", u" Clear", None))
        self.ResultsTitle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\"\n"
"                  font-size:20pt;\n"
"                  font-weight:600;\">Results</span></p></body></html>", None))
        self.tabWidgetResult.setTabText(self.tabWidgetResult.indexOf(self.resultTab), QCoreApplication.translate("MainWindow", u"Results", None))
#if QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle.setText(QCoreApplication.translate("MainWindow", u"Condition", None))
#if QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle_2.setWhatsThis("")
#endif // QT_CONFIG(whatsthis)
        self.ResultsStatsListTitle_2.setText(QCoreApplication.translate("MainWindow", u"trees", None))
        self.ResultsStatsTitle.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\"\n"
"                    font-size:24pt; text-decoration:\n"
"                    underline;\">Statistics</span></p></body></html>", None))
        self.PhyloTreeLabel.setText("")
        self.downloadResultsPlotButton.setText(QCoreApplication.translate("MainWindow", u"Download", None))
        self.tabWidgetResult.setTabText(self.tabWidgetResult.indexOf(self.resultStatTab), QCoreApplication.translate("MainWindow", u"Statistics", None))
        self.settingsButtonPage4.setText(QCoreApplication.translate("MainWindow", u" Settings", None))
        self.submitButtonPage4.setText(QCoreApplication.translate("MainWindow", u"Save", None))
        self.statisticsButtonPage4.setText(QCoreApplication.translate("MainWindow", u" Statistics", None))
        self.clearButtonPage4.setText(QCoreApplication.translate("MainWindow", u" Clear", None))
    # retranslateUi

