import yaml
from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets
from aphylogeo.params import Params

Params.load_from_file("params.yaml")

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

class UiDialog(object):
    selected_index = 0

    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(600, 800)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        Dialog.setMinimumSize(QtCore.QSize(600, 800))
        Dialog.setMaximumSize(QtCore.QSize(600, 800))

        self.formLayout = QtWidgets.QFormLayout(Dialog)
        self.formLayout.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setLabelAlignment(QtCore.Qt.AlignLeft)
        self.formLayout.setFormAlignment(QtCore.Qt.AlignTop)
        self.formLayout.setSpacing(20)

        self.createLabelAndLineEdit(Dialog, "File Name", "file_name", Params.file_name, read_only=True)
        self.createLabelAndLineEdit(Dialog, "Specimen", "specimen", Params.specimen)
        self.createLabelAndLineEdit(Dialog, "Names", "names", ", ".join(Params.names), read_only=True)
        self.createLabelAndSpinBox(Dialog, "Bootstrap Threshold", "bootstrap_threshold", Params.bootstrap_threshold)
        self.createLabelAndSpinBox(Dialog, "Dist Threshold", "dist_threshold", Params.dist_threshold)
        self.createLabelAndSpinBox(Dialog, "Window Size", "window_size", Params.window_size)
        self.createLabelAndSpinBox(Dialog, "Step Size", "step_size", Params.step_size)
        self.createLabelAndSpinBox(Dialog, "Bootstrap Amount", "bootstrap_amount", Params.bootstrap_amount)
        self.createLabelAndLineEdit(Dialog, "Data Names", "data_names", ", ".join(Params.data_names), read_only=True)
        self.createLabelAndLineEdit(Dialog, "Reference Gene Dir", "reference_gene_dir", Params.reference_gene_dir, read_only=True)
        self.createLabelAndLineEdit(Dialog, "Reference Gene File", "reference_gene_file", Params.reference_gene_file, read_only=True)
        self.createLabelAndCheckBox(Dialog, "Make Debug Files", "makeDebugFiles", Params.makeDebugFiles)
        self.createLabelAndLineEdit(Dialog, "Alignment Method", "alignment_method", Params.alignment_method)
        self.createLabelAndLineEdit(Dialog, "Distance Method", "distance_method", str(Params.distance_method))
        self.createLabelAndLineEdit(Dialog, "Fit Method", "fit_method", Params.fit_method)
        self.createLabelAndLineEdit(Dialog, "Tree Type", "tree_type", Params.tree_type)
        self.createLabelAndSpinBox(Dialog, "Rate Similarity", "rate_similarity", Params.rate_similarity)
        self.createLabelAndLineEdit(Dialog, "Method Similarity", "method_similarity", Params.method_similarity)

        self.cancel_button = QtWidgets.QPushButton(Dialog)
        self.cancel_button.setText("Cancel")
        self.cancel_button.clicked.connect(Dialog.close)
        self.formLayout.addRow(self.cancel_button)

        self.ok_button = QtWidgets.QPushButton(Dialog)
        self.ok_button.setText("Ok")
        self.ok_button.clicked.connect(self.saveData)
        self.ok_button.clicked.connect(Dialog.close)
        self.formLayout.addRow(self.ok_button)

        self.translateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def createLabelAndLineEdit(self, Dialog, label_text, obj_name, value, read_only=False):
        label = QtWidgets.QLabel(Dialog)
        label.setText(label_text)
        line_edit = QtWidgets.QLineEdit(Dialog)
        line_edit.setObjectName(obj_name)
        line_edit.setText(value)
        if read_only:
            line_edit.setReadOnly(True)
        line_edit.setMinimumWidth(300)
        self.formLayout.addRow(label, line_edit)

    def createLabelAndSpinBox(self, Dialog, label_text, obj_name, value):
        label = QtWidgets.QLabel(Dialog)
        label.setText(label_text)
        spin_box = QtWidgets.QSpinBox(Dialog)
        spin_box.setObjectName(obj_name)
        spin_box.setValue(value)
        spin_box.setMinimumWidth(300)
        self.formLayout.addRow(label, spin_box)

    def createLabelAndCheckBox(self, Dialog, label_text, obj_name, checked):
        label = QtWidgets.QLabel(Dialog)
        label.setText(label_text)
        check_box = QtWidgets.QCheckBox(Dialog)
        check_box.setObjectName(obj_name)
        check_box.setChecked(checked)
        self.formLayout.addRow(label, check_box)

    def translateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Parameters"))

    def saveData(self):
        Dialog = self.formLayout.parentWidget()

        fields = {
            "file_name": Dialog.findChild(QtWidgets.QLineEdit, "file_name").text(),
            "specimen": Dialog.findChild(QtWidgets.QLineEdit, "specimen").text(),
            "names": Dialog.findChild(QtWidgets.QLineEdit, "names").text().split(", "),
            "bootstrap_threshold": Dialog.findChild(QtWidgets.QSpinBox, "bootstrap_threshold").value(),
            "dist_threshold": Dialog.findChild(QtWidgets.QSpinBox, "dist_threshold").value(),
            "window_size": Dialog.findChild(QtWidgets.QSpinBox, "window_size").value(),
            "step_size": Dialog.findChild(QtWidgets.QSpinBox, "step_size").value(),
            "bootstrap_amount": Dialog.findChild(QtWidgets.QSpinBox, "bootstrap_amount").value(),
            "data_names": Dialog.findChild(QtWidgets.QLineEdit, "data_names").text().split(", "),
            "reference_gene_dir": Dialog.findChild(QtWidgets.QLineEdit, "reference_gene_dir").text(),
            "reference_gene_file": Dialog.findChild(QtWidgets.QLineEdit, "reference_gene_file").text(),
            "makeDebugFiles": Dialog.findChild(QtWidgets.QCheckBox, "makeDebugFiles").isChecked(),
            "alignment_method": Dialog.findChild(QtWidgets.QLineEdit, "alignment_method").text(),
            "distance_method": int(Dialog.findChild(QtWidgets.QLineEdit, "distance_method").text()),
            "fit_method": Dialog.findChild(QtWidgets.QLineEdit, "fit_method").text(),
            "tree_type": Dialog.findChild(QtWidgets.QLineEdit, "tree_type").text(),
            "rate_similarity": Dialog.findChild(QtWidgets.QSpinBox, "rate_similarity").value(),
            "method_similarity": Dialog.findChild(QtWidgets.QLineEdit, "method_similarity").text(),
        }

        for property_name, new_value in fields.items():
            update_yaml_param(Params, "params.yaml", property_name, new_value)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = UiDialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())