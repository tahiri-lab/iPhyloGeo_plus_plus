from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QLabel
from PyQt6.QtWidgets import QApplication, QDialog, QGridLayout

from utils.hover_label import HoverLabel
from utils.my_dumper import update_yaml_param

from aphylogeo.params import Params

class ResultSettingsDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Settings")
        self.setGeometry(100, 100, 800, 600)
        self.initUI()
        self.apply_styles()

    def initUI(self):  # noqa: N803
        self.setObjectName("Dialog")
        self.resize(1000, 600)  # Decreased overall size

        layout = QGridLayout()

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred, QtWidgets.QSizePolicy.Policy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)

        layout.addItem(self.setup_buttons(), 1, 0)

        content_layout = QGridLayout()
        content_layout.setHorizontalSpacing(50)

        self.userParams, self.userParamsLayout = create_box_and_layout(self, 11, "userParams")

        content_layout.addWidget(self.userParams, 0, 0)

        self.paramsDetails, paramsDetailsLayout = create_box_and_layout(self, 11, "paramsDetails")

        content_layout.addWidget(self.paramsDetails, 0, 1)
        
        HoverLabel.image_label = QLabel()
        HoverLabel.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        HoverLabel.image_label.setPixmap(QPixmap("./img/other/final.png"))
        paramsDetailsLayout.addWidget(HoverLabel.image_label)

        self.textEdit = QtWidgets.QTextEdit(self.paramsDetails)
        
        sizePolicy.setHeightForWidth(self.textEdit.sizePolicy().hasHeightForWidth())
        self.textEdit.setSizePolicy(sizePolicy)

        self.textEdit.setReadOnly(True)
        self.textEdit.setObjectName("textEdit")
        
        paramsDetailsLayout.addWidget(self.textEdit)

        self.methodBox, self.methodBoxLayout  = create_box_and_layout(self.userParams, 9, "methodBox")
  
        self.metrics = HoverLabel(
            "Calculus method",
            "Select the method for calculating phylogenetic distances. Options include:\n- Robinson & Foulds: Measures the difference in tree topologies.\n- Least Square: Minimizes the sum of squared differences between observed and predicted distances.\n- Euclidean Distance: Measures the straight-line distance between points in a multi-dimensional space.",
            self.textEdit,
            HoverLabel.image_label,
            "./img/other/calculus.png",
        )
        self.methodBoxLayout.addWidget(self.metrics, 0, 0, 1, 1)

        self.comboBox_metrics = QtWidgets.QComboBox(self.methodBox)
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred, QtWidgets.QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.comboBox_metrics.sizePolicy().hasHeightForWidth())  
        self.comboBox_metrics.setSizePolicy(sizePolicy)
        
        self.comboBox_metrics.setObjectName("comboBox_metrics")
        self.comboBox_metrics.addItem("")
        self.comboBox_metrics.addItem("")
        self.comboBox_metrics.addItem("")
        
        self.methodBoxLayout.addWidget(self.comboBox_metrics, 0, 1, 1, 1)

        self.thresholdBox, self.thresholdBoxLayout = create_box_and_layout(self.userParams, 9, "thresholdBox")
       
        self.spinBox_bootstrap = QtWidgets.QSpinBox(self.thresholdBox)
        
        self.spinBox_bootstrap.setMinimum(0)
        self.spinBox_bootstrap.setMaximum(1000)
        self.spinBox_bootstrap.setObjectName("spinBox_bootstrap")
        self.spinBox_bootstrap.setFixedHeight(25)
        
        self.thresholdBoxLayout.addWidget(self.spinBox_bootstrap, 0, 1, 1, 1)

        self.bootstrapValue = HoverLabel(
            "Bootstrap threshold",
            "Set the bootstrap threshold for phylogenetic tree support. Higher values provide more reliable results but require longer computation time. Bootstrap values are used to assess the reliability of inferred phylogenetic trees.",
            self.textEdit,
            HoverLabel.image_label,
            "./img/other/bootstrap.png",
        )
        self.thresholdBoxLayout.addWidget(self.bootstrapValue, 0, 0, 1, 1)

        self.spinBox_metricThreshold = QtWidgets.QSpinBox(self.thresholdBox)
        self.spinBox_metricThreshold.setObjectName("spinBox_metricThreshold")
        self.spinBox_metricThreshold.setFixedHeight(25)
        
        self.thresholdBoxLayout.addWidget(self.spinBox_metricThreshold, 1, 1, 1, 1)

        self.metricThreshold = HoverLabel(
            "Metric threshold",
            "Set the threshold for distance metrics. Higher thresholds provide more accurate results but increase computational load. This parameter controls the cutoff for considering distances in the analysis.",
            self.textEdit,
            HoverLabel.image_label,
            "./img/other/metric.png",
        )
        self.thresholdBoxLayout.addWidget(self.metricThreshold, 1, 0, 1, 1)
        
        
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetNoConstraint)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setObjectName("verticalLayout")
        
        self.verticalLayout.addWidget(self.methodBox)
        self.verticalLayout.addWidget(self.thresholdBox)
        self.userParamsLayout.addLayout(self.verticalLayout, 5, 0, 1, 1)

        self.resetValues()

        layout.addItem(content_layout, 0, 0)

        self.setLayout(layout)
        self.retranslateUi()
        self.show()
        
    def setup_buttons(self):
        
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
        
        
        self.ok_button.clicked.connect(self.saveData)
        self.ok_button.clicked.connect(self.close)
        self.cancel_button.clicked.connect(self.close)
        self.reset_button.clicked.connect(self.resetValues)
        
        return button_layout
    
        

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
        self.textEdit.setText(_translate("Dialog", "Hover over the parameter you want details about..."))

    def saveData(self):
        metrics = self.comboBox_metrics.currentIndex()
        bootstrap_value = self.spinBox_bootstrap.value()
        metric_threshold = self.spinBox_metricThreshold.value()

        update_yaml_param(Params, "params.yaml", "distance_method", str(metrics))
        update_yaml_param(Params, "params.yaml", "bootstrap_threshold", bootstrap_value)
        update_yaml_param(Params, "params.yaml", "dist_threshold", metric_threshold)

    def resetValues(self):
        self.comboBox_metrics.setCurrentIndex(int(Params.distance_method))
        self.spinBox_metricThreshold.setProperty("value", Params.dist_threshold)
        self.spinBox_bootstrap.setProperty("value", Params.bootstrap_threshold)

    def apply_styles(self):  # noqa: N803
        self.setStyleSheet(
            """
            QDialog {
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            }
            QGroupBox {
                border: 1px solid #d1d1d1;
                border-radius: 5px;
                margin-top: 5px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 3px 5px;
                font-size: 13px;
                font-weight: bold;
            }
            QLabel {
                font-size: 12px;
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
                font-size: 12px;
            }
            QComboBox {
                border: 1px solid #ced4da;
                border-radius: 5px;
                padding: 5px;
                font-size: 12px;
            }
            QSpinBox {
                border: 1px solid #ced4da;
                border-radius: 5px;
                padding: 5px;
                font-size: 12px;
            }
        """
        )



def create_box_and_layout(parent, fontSize, name):
    
    groupBox = QtWidgets.QGroupBox(parent)

    sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred, QtWidgets.QSizePolicy.Policy.Preferred)
    sizePolicy.setHorizontalStretch(0)
    sizePolicy.setVerticalStretch(0)
    sizePolicy.setHeightForWidth(groupBox.sizePolicy().hasHeightForWidth())
    groupBox.setSizePolicy(sizePolicy)    

    font = QtGui.QFont()
    font.setPointSize(fontSize)
    groupBox.setFont(font)

    groupBox.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
    groupBox.setCheckable(False)
    groupBox.setObjectName(name)
    
    layout = QtWidgets.QGridLayout(groupBox)
    if fontSize == 11:
        layout.setContentsMargins(8, 20, 8, 8)
    else:
        layout.setContentsMargins(5, 5, 5, 5)
        
    layout.setSpacing(5)
    layout.setObjectName(name + "layout")
    
    return groupBox, layout
    
if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    window = ResultSettingsDialog()
    window.show()
    sys.exit(app.exec())
