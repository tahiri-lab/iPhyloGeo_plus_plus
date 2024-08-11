import os
import yaml
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QLabel
from aphylogeo.params import Params
from yaml.loader import SafeLoader

class Params2:
    PARAMETER_KEYS = {
        "bootstrap_threshold": 0,
        "dist_threshold": 60,
        "window_size": 20,
        "step_size": 100,
        "data_names": ["ALLSKY_SFC_SW_DWN_newick", "T2M_newick", "QV2M_newick", "PRECTOTCORR_newick", "WS10M_newick"],
        "reference_gene_dir": "./datasets/example",
        "reference_gene_file": "sequences.fasta",
        "file_name": "./datasets/example/geo.csv",
        "specimen": "id",
        "names": ["id", "ALLSKY_SFC_SW_DWN", "T2M", "PRECTOTCORR", "QV2M", "WS10M"],
        "makeDebugFiles": False,
        "bootstrap_amount": 100,
        "alignment_method": "2",
        "distance_method": "0",
        "fit_method": "1",
        "tree_type": "2",
        "rate_similarity": 90,
        "method_similarity": "1",
    }

    @classmethod
    def load_from_file(cls, params_file=os.path.join(os.path.dirname(__file__), "./utils/params.yaml")):
        with open(params_file) as f:
            params = yaml.load(f, Loader=SafeLoader)
            cls.validate_and_set_params(params)

    @classmethod
    def update_from_dict(cls, params_content):
        cls.validate_and_set_params(params_content)

    @classmethod
    def validate_and_set_params(cls, params_dict):
        for key, value in params_dict.items():
            if key in cls.PARAMETER_KEYS:
                setattr(cls, key, value)
            else:
                raise ValueError(f"Invalid parameter: {key}")

        if hasattr(cls, "reference_gene_dir") and hasattr(cls, "reference_gene_file"):
            cls.reference_gene_filepath = os.path.join(cls.reference_gene_dir, cls.reference_gene_file)
        else:
            cls.reference_gene_filepath = None

Params.load_from_file("./utils/params.yaml")
Params2.load_from_file("./utils/params_default.yaml")

class MyDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(MyDumper, self).increase_indent(flow, False)

    def represent_list(self, data):
        return self.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

yaml.add_representer(list, MyDumper.represent_list, Dumper=MyDumper)

def update_yaml_param(params, file_path, property_name, new_value):
    if isinstance(new_value, list):
        new_value = [element.strip() for element in new_value]
    params.update_from_dict({property_name: new_value})

    with open(file_path, "r") as yaml_file:
        data = yaml.safe_load(yaml_file)

    if property_name in data:
        data[property_name] = new_value
    else:
        print(f"Warning: Property '{property_name}' not found in '{file_path}'.")

    with open(file_path, "w") as yaml_file:
        yaml.dump(data, yaml_file, default_flow_style=None, Dumper=MyDumper, sort_keys=False)

class HoverLabel(QLabel):
    def __init__(self, text, hover_text, text_edit, image_label, hover_image_path, *args, **kwargs):
        super().__init__(text, *args, **kwargs)
        self.default_text = text
        self.hover_text = hover_text
        self.text_edit = text_edit
        self.image_label = image_label
        self.hover_image_path = hover_image_path

    def enterEvent(self, event):
        self.text_edit.clear()
        self.image_label.clear()
        self.text_edit.setText(self.hover_text)
        self.image_label.setPixmap(QPixmap(self.hover_image_path))
        super().enterEvent(event)

class Settings(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(700, 400)  # Decreased overall size
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        Dialog.setMinimumSize(QtCore.QSize(700, 400))  # Adjusted minimum size
        Dialog.setMaximumSize(QtCore.QSize(700, 400))  # Adjusted maximum size

        self.apply_styles(Dialog)

        self.reset_button = QtWidgets.QPushButton(Dialog)
        self.reset_button.setGeometry(QtCore.QRect(180, 350, 100, 30))
        self.reset_button.setObjectName("reset_button")
        self.ok_button = QtWidgets.QPushButton(Dialog)
        self.ok_button.setGeometry(QtCore.QRect(290, 350, 100, 30))
        self.ok_button.setObjectName("ok_button")
        self.cancel_button = QtWidgets.QPushButton(Dialog)
        self.cancel_button.setGeometry(QtCore.QRect(400, 350, 100, 30))
        self.cancel_button.setObjectName("cancel_button")

        self.userParams = QtWidgets.QGroupBox(Dialog)
        self.userParams.setGeometry(QtCore.QRect(10, 10, 260, 330))
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
        self.gridLayout_2.setContentsMargins(8, 8, 8, 8)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setObjectName("gridLayout_2")

        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setObjectName("verticalLayout")

        self.paramsDetails = QtWidgets.QGroupBox(Dialog)
        self.paramsDetails.setGeometry(QtCore.QRect(280, 10, 410, 330))
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
        HoverLabel.image_label.setFixedSize(380, 150)
        HoverLabel.image_label.setAlignment(Qt.AlignCenter)
        HoverLabel.image_label.setPixmap(QPixmap("../img/other/final.png"))
        group_box_layout.addWidget(HoverLabel.image_label)

        self.textEdit = QtWidgets.QTextEdit(self.paramsDetails)
        group_box_layout.addWidget(self.textEdit)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.textEdit.sizePolicy().hasHeightForWidth())
        self.textEdit.setSizePolicy(sizePolicy)
        self.textEdit.setMinimumSize(QtCore.QSize(380, 130))
        self.textEdit.setMaximumSize(QtCore.QSize(380, 130))
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

        self.metrics = HoverLabel("Calculus method", "Select the method for calculating phylogenetic distances. Options include:\n- Robinson & Foulds: Measures the difference in tree topologies.\n- Least Square: Minimizes the sum of squared differences between observed and predicted distances.\n- Euclidean Distance: Measures the straight-line distance between points in a multi-dimensional space.", self.textEdit, HoverLabel.image_label, "../img/other/calculus.png")
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

        self.bootstrapValue = HoverLabel("Bootstrap threshold", "Set the bootstrap threshold for phylogenetic tree support. Higher values provide more reliable results but require longer computation time. Bootstrap values are used to assess the reliability of inferred phylogenetic trees.", self.textEdit, HoverLabel.image_label, "../img/other/bootstrap.png")
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

        self.metricThreshold = HoverLabel("Metric threshold", "Set the threshold for distance metrics. Higher thresholds provide more accurate results but increase computational load. This parameter controls the cutoff for considering distances in the analysis.", self.textEdit, HoverLabel.image_label, "../img/other/metric.png")
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

        self.ok_button.clicked.connect(self.saveData)
        self.ok_button.clicked.connect(Dialog.close)
        self.cancel_button.clicked.connect(Dialog.close)
        self.reset_button.clicked.connect(self.resetValues)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Parameters"))
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

        update_yaml_param(Params, "./utils/params.yaml", "distance_method", str(metrics))
        update_yaml_param(Params, "./utils/params.yaml", "bootstrap_threshold", bootstrap_value)
        update_yaml_param(Params, "./utils/params.yaml", "dist_threshold", metric_threshold)

    def resetValues(self):
        self.comboBox_metrics.setProperty("value", Params2.distance_method)
        self.spinBox_metricThreshold.setProperty("value", Params2.dist_threshold)
        self.spinBox_bootstrap.setProperty("value", Params2.bootstrap_threshold)

    def apply_styles(self, Dialog):
        Dialog.setStyleSheet("""
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
        """)

def comboBoxSelected(self, index):
    if index != 0:
        HoverLabel.Settings.selected_index = index

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Settings()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
