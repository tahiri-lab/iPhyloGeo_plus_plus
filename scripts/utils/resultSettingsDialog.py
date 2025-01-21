from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QApplication, QDialog, QGridLayout, QLabel
from utils.settings import HoverLabel
from utils.settings import Params2
from utils.MyDumper import update_yaml_param


from aphylogeo.params import Params


class ResultSettingsDialog(QDialog):
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Settings')
        self.setGeometry(100, 100, 800, 600)
        self.initUI()
        self.apply_styles()
        
    def initUI(self):  # noqa: N803
        self.setObjectName("Dialog")
        self.resize(1000, 600)  # Decreased overall size
        
        layout = QGridLayout()

        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)


        button_layout = QGridLayout()
        button_layout.setHorizontalSpacing(30)
        button_layout.setVerticalSpacing(20)
        button_layout.setContentsMargins(50, 30, 50, 0)
        
        self.reset_button = QtWidgets.QPushButton()
        self.reset_button.setObjectName("reset_button")
        button_layout.addWidget(self.reset_button, 0, 0)
        
        self.ok_button = QtWidgets.QPushButton()
        self.ok_button.setObjectName("ok_button")
        button_layout.addWidget(self.ok_button, 0, 1)
        
        self.cancel_button = QtWidgets.QPushButton()
        self.cancel_button.setObjectName("cancel_button")
        button_layout.addWidget(self.cancel_button, 0, 2)
               
        layout.addItem(button_layout, 1, 0)
        
        
        content_layout = QGridLayout()
        
        self.userParams = QtWidgets.QGroupBox(self)
        content_layout.addWidget(self.userParams, 0, 0)
        content_layout.setHorizontalSpacing(50)
        
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)  
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.userParams.sizePolicy().hasHeightForWidth())
        self.userParams.setSizePolicy(sizePolicy)
        
        font = QtGui.QFont()
        font.setPointSize(11)
        self.userParams.setFont(font)
        
        self.userParams.setAlignment(QtCore.Qt.AlignCenter)
        self.userParams.setCheckable(False)
        self.userParams.setObjectName("userParams")

        self.gridLayout_2 = QtWidgets.QGridLayout(self.userParams)
        self.gridLayout_2.setContentsMargins(8, 20, 8, 8)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setObjectName("gridLayout_2")

        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setObjectName("verticalLayout")

        self.paramsDetails = QtWidgets.QGroupBox(self)
        content_layout.addWidget(self.paramsDetails, 0, 1)
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.paramsDetails.sizePolicy().hasHeightForWidth())
        self.paramsDetails.setSizePolicy(sizePolicy)
        
        font.setPointSize(11)
        self.paramsDetails.setFont(font)
        
        self.paramsDetails.setAlignment(QtCore.Qt.AlignCenter)
        self.paramsDetails.setCheckable(False)
        self.paramsDetails.setObjectName("paramsDetails")
        
        group_box_layout = QtWidgets.QVBoxLayout()
        self.paramsDetails.setLayout(group_box_layout)
        
        HoverLabel.image_label = QLabel()
        HoverLabel.image_label.setAlignment(Qt.AlignCenter)
        HoverLabel.image_label.setPixmap(QPixmap("./img/other/final.png"))
        group_box_layout.addWidget(HoverLabel.image_label)
        
        self.textEdit = QtWidgets.QTextEdit(self.paramsDetails)
        group_box_layout.addWidget(self.textEdit)
        
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.textEdit.sizePolicy().hasHeightForWidth())
        self.textEdit.setSizePolicy(sizePolicy)
        
        self.textEdit.setMinimumSize(QtCore.QSize(380, 130))
        
        font.setPointSize(8)
        self.textEdit.setFont(font)
        
        self.textEdit.setReadOnly(True)
        self.textEdit.setObjectName("textEdit")

        self.methodBox = QtWidgets.QGroupBox(self.userParams)
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.methodBox.sizePolicy().hasHeightForWidth())
        self.methodBox.setSizePolicy(sizePolicy)
        
        font.setPointSize(9)
        self.methodBox.setFont(font)
        
        self.methodBox.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.methodBox.setAlignment(QtCore.Qt.AlignCenter)
        self.methodBox.setFlat(True)
        self.methodBox.setCheckable(False)
        self.methodBox.setObjectName("methodBox")
        self.verticalLayout.addWidget(self.methodBox)

        self.gridLayout_3 = QtWidgets.QGridLayout(self.methodBox)
        self.gridLayout_3.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.gridLayout_3.setContentsMargins(5, 5, 5, 5)
        self.gridLayout_3.setSpacing(5)
        self.gridLayout_3.setObjectName("gridLayout_3")
        font.setPointSize(8)

        self.metrics = HoverLabel(
            "Calculus method",
            "Select the method for calculating phylogenetic distances. Options include:\n- Robinson & Foulds: Measures the difference in tree topologies.\n- Least Square: Minimizes the sum of squared differences between observed and predicted distances.\n- Euclidean Distance: Measures the straight-line distance between points in a multi-dimensional space.",
            self.textEdit,
            HoverLabel.image_label,
            "./img/other/calculus.png",
        )
        self.metrics.setFont(font)
        self.metrics.setIndent(10)
        self.metrics.setObjectName("metrics")
        self.gridLayout_3.addWidget(self.metrics, 0, 0, 1, 1)

        self.comboBox_metrics = QtWidgets.QComboBox(self.methodBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.comboBox_metrics.sizePolicy().hasHeightForWidth())
        self.comboBox_metrics.setSizePolicy(sizePolicy)
        font.setPointSize(8)
        self.comboBox_metrics.setFont(font)
        self.comboBox_metrics.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.comboBox_metrics.setIconSize(QtCore.QSize(20, 20))
        self.comboBox_metrics.setFrame(True)
        self.comboBox_metrics.setObjectName("comboBox_metrics")
        self.comboBox_metrics.addItem("")
        self.comboBox_metrics.addItem("")
        self.comboBox_metrics.addItem("")
        self.gridLayout_3.addWidget(self.comboBox_metrics, 0, 1, 1, 1)

        self.thresholdBox = QtWidgets.QGroupBox(self.userParams)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.thresholdBox.sizePolicy().hasHeightForWidth())
        self.thresholdBox.setSizePolicy(sizePolicy)
        font.setPointSize(9)
        self.thresholdBox.setFont(font)
        self.thresholdBox.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.thresholdBox.setAlignment(QtCore.Qt.AlignCenter)
        self.thresholdBox.setFlat(True)
        self.thresholdBox.setCheckable(False)
        self.thresholdBox.setObjectName("thresholdBox")

        self.gridLayout_4 = QtWidgets.QGridLayout(self.thresholdBox)
        self.gridLayout_4.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.gridLayout_4.setContentsMargins(5, 5, 5, 5)
        self.gridLayout_4.setSpacing(5)
        self.gridLayout_4.setObjectName("gridLayout_4")

        self.spinBox_bootstrap = QtWidgets.QSpinBox(self.thresholdBox)
        font.setPointSize(8)
        self.spinBox_bootstrap.setFont(font)
        self.spinBox_bootstrap.setMinimum(0)
        self.spinBox_bootstrap.setMaximum(1000)
        self.spinBox_bootstrap.setSingleStep(1)
        self.spinBox_bootstrap.setObjectName("spinBox_bootstrap")
        self.gridLayout_4.addWidget(self.spinBox_bootstrap, 0, 1, 1, 1)

        self.bootstrapValue = HoverLabel(
            "Bootstrap threshold",
            "Set the bootstrap threshold for phylogenetic tree support. Higher values provide more reliable results but require longer computation time. Bootstrap values are used to assess the reliability of inferred phylogenetic trees.",
            self.textEdit,
            HoverLabel.image_label,
            "./img/other/bootstrap.png",
        )
        font.setPointSize(8)
        self.bootstrapValue.setFont(font)
        self.bootstrapValue.setIndent(10)
        self.bootstrapValue.setObjectName("bootstrapValue")
        self.gridLayout_4.addWidget(self.bootstrapValue, 0, 0, 1, 1)

        self.spinBox_metricThreshold = QtWidgets.QSpinBox(self.thresholdBox)
        font.setPointSize(8)
        self.spinBox_metricThreshold.setFont(font)
        self.spinBox_metricThreshold.setProperty("value", 60)
        self.spinBox_metricThreshold.setObjectName("spinBox_metricThreshold")
        self.gridLayout_4.addWidget(self.spinBox_metricThreshold, 1, 1, 1, 1)

        self.metricThreshold = HoverLabel(
            "Metric threshold",
            "Set the threshold for distance metrics. Higher thresholds provide more accurate results but increase computational load. This parameter controls the cutoff for considering distances in the analysis.",
            self.textEdit,
            HoverLabel.image_label,
            "./img/other/metric.png",
        )
        font.setPointSize(8)
        self.metricThreshold.setFont(font)
        self.metricThreshold.setIndent(10)
        self.metricThreshold.setObjectName("metricThreshold")
        self.gridLayout_4.addWidget(self.metricThreshold, 1, 0, 1, 1)
        self.verticalLayout.addWidget(self.thresholdBox)
        self.gridLayout_2.addLayout(self.verticalLayout, 5, 0, 1, 1)

        self.comboBox_metrics.setProperty("value", Params.distance_method)
        self.spinBox_metricThreshold.setProperty("value", Params.dist_threshold)
        self.spinBox_bootstrap.setProperty("value", Params.bootstrap_threshold)

        layout.addItem(content_layout, 0, 0)
        
        self.ok_button.clicked.connect(self.saveData)
        self.ok_button.clicked.connect(self.close)
        self.cancel_button.clicked.connect(self.close)
        self.reset_button.clicked.connect(self.resetValues)

        self.setLayout(layout)
        self.retranslateUi()
        self.show()
        #  QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self):  # noqa: N803
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("Dialog", "Parameters"))
        self.cancel_button.setText(_translate("Dialog", "Cancel"))
        self.reset_button.setText(_translate("Dialog", "Reset"))
        self.ok_button.setText(_translate("Dialog", "OK"))
        self.userParams.setTitle(_translate("Dialog", "User Parameters"))
        self.methodBox.setTitle(_translate("Dialog", "Method Selection"))
        self.metrics.setText(_translate("Dialog", "Calculus method"))
        self.comboBox_metrics.setItemText(0, _translate("Dialog", "Least square"))
        self.comboBox_metrics.setItemText(1, _translate("Dialog", "Robinson and Foulds"))
        self.comboBox_metrics.setItemText(2, _translate("Dialog", "Euclidean distance"))
        self.thresholdBox.setTitle(_translate("Dialog", "Threshold Selection"))
        self.bootstrapValue.setText(_translate("Dialog", "Bootstrap threshold"))
        self.metricThreshold.setText(_translate("Dialog", "Metric threshold"))
        self.paramsDetails.setTitle(_translate("Dialog", "Details"))
        self.textEdit.setPlaceholderText(_translate("Dialog", "Hover over the parameter you want details about..."))

    def saveData(self):
        metrics = self.comboBox_metrics.currentIndex()
        bootstrap_value = self.spinBox_bootstrap.value()
        metric_threshold = self.spinBox_metricThreshold.value()

        update_yaml_param(Params, "./scripts/utils/params.yaml", "distance_method", str(metrics))
        update_yaml_param(Params, "./scripts/utils/params.yaml", "bootstrap_threshold", bootstrap_value)
        update_yaml_param(Params, "./scripts/utils/params.yaml", "dist_threshold", metric_threshold)

    def resetValues(self):
        self.comboBox_metrics.setProperty("value", Params2.distance_method)
        self.spinBox_metricThreshold.setProperty("value", Params2.dist_threshold)
        self.spinBox_bootstrap.setProperty("value", Params2.bootstrap_threshold)

    def apply_styles(self):  # noqa: N803
        self.setStyleSheet(
            """
            QDialog {
                background-color: #f8f9fa;
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            }
            QGroupBox {
                border: 1px solid #d1d1d1;
                border-radius: 5px;
                margin-top: 5px;
                background-color: #fff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 3px 5px;
                font-size: 13px;
                font-weight: bold;
                color: #343a40;
            }
            QLabel {
                font-size: 12px;
                color: #495057;
            }
            QPushButton {
                background-color: #007bff;
                color: white;
                border: none;
                border-radius: 5px;
                padding: 5px 10px;
                font-size: 12px;
            }
            QPushButton#cancel_button {
                background-color: #dc3545;
            }
            QPushButton:hover {
                background-color: #0056b3;
            }
            QPushButton#cancel_button:hover {
                background-color: #c82333;
            }
            QTextEdit {
                border: 1px solid #ced4da;
                border-radius: 5px;
                padding: 5px;
                background-color: #fff;
                color: #495057;
                font-size: 12px;
            }
            QComboBox {
                border: 1px solid #ced4da;
                border-radius: 5px;
                padding: 5px;
                background-color: #fff;
                color: #495057;
                font-size: 12px;
            }
            QSpinBox {
                border: 1px solid #ced4da;
                border-radius: 5px;
                padding: 5px;
                background-color: #fff;
                color: #495057;
                font-size: 12px;
            }
        """
        )


if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    window = ResultSettingsDialog()
    window.show()
    sys.exit(app.exec_())
