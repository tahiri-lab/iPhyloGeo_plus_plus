import os

import yaml
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication, QLabel
from yaml.loader import SafeLoader
from aphylogeo.params import Params


class Params2:
    PARAMETER_KEYS = {
        "file_name": "./datasets/simplot.csv",
        "specimen": "id",
        "names": ["id", "ALLSKY_SFC_SW_DWN", "T2M", "QV2M", "PRECTOTCORR", "WS10M", "LAT", "LONG"],
        "bootstrap_threshold": 0,
        "dist_threshold": 60,
        "window_size": 200,
        "step_size": 100,
        "bootstrap_amount": 100,
        "data_names": ["ALLSKY_SFC_SW_DWN", "T2M", "QV2M", "PRECTOTCORR", "WS10M", "LAT", "LONG"],
        "reference_gene_dir": "./datasets",
        "reference_gene_file": "simplot.fasta",
        "makeDebugFiles": True,
        "alignment_method": "4",
        "distance_method": "3",
        "fit_method": "1",
        "tree_type": "1",
        "rate_similarity": 70,
        "method_similarity": "1",
    }

    @classmethod
    def load_from_file(cls, params_file=os.path.join(os.path.dirname(__file__), "./scripts/utils/params.yaml")):  # noqa: B008
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


try:
    Params.load_from_file("./scripts/utils/params.yaml")
    Params2.load_from_file("./scripts/utils/params_default.yaml")
except FileNotFoundError:
    Params.validate_and_set_params(Params2.PARAMETER_KEYS)
    Params2.validate_and_set_params(Params2.PARAMETER_KEYS)

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

if __name__ == "__main__":
    import sys
    
    app = QApplication(sys.argv)
    sys.exit(app.exec_())
