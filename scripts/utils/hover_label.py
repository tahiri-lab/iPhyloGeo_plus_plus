
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QLabel


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
