from PyQt5 import QtGui, QtCore, QtWidgets


def style_buttons(buttons, dark_mode):
    for button in buttons:
        button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        button.setStyleSheet(get_button_style(dark_mode))
        button.setGraphicsEffect(create_shadow_effect(10, 140))


def get_button_style(dark_mode):
    if dark_mode:
        background_color = "#646464"
        hover_color = "#B7B7B6"
    else:
        background_color = "#DEDDDA"
        hover_color = "#B7B7B6"

    return f"""
        QPushButton {{
            padding: 10px 20px;
            font-weight: bold;
            background-color: {background_color};
            border-radius: 20px;
            transition: background-color 0.3s ease;
        }}
        QPushButton:hover {{
            background-color: {hover_color};
        }}
    """


def create_shadow_effect(blur_radius, alpha):
    shadow_effect = QtWidgets.QGraphicsDropShadowEffect()
    shadow_effect.setBlurRadius(blur_radius)
    shadow_effect.setColor(QtGui.QColor(0, 0, 0, alpha))
    shadow_effect.setOffset(3, 3)
    return shadow_effect
