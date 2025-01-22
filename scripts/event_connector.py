from enum import Enum
from typing import List, Union

from scripts.main import UiMainWindow


class QtEvents(str, Enum):
    # Buttons
    clicked = "clicked"
    pressed = "pressed"
    released = "released"
    toggled = "toggled"
    enterEvent = "enterEvent"
    leaveEvent = "leaveEvent"

    # LineEdit
    textChanged = "textChanged"
    textEdited = "textEdited"
    returnPressed = "returnPressed"
    editingFinished = "editingFinished"

    # ComboBox
    currentIndexChanged = "currentIndexChanged"
    currentTextChanged = "currentTextChanged"

    # Slider
    valueChanged = "valueChanged"


def connect_event(widgets: Union[str, List[str]], event: QtEvents):
    """Connects an event to a widget or a list of widgets"""

    if isinstance(widgets, str):
        widgets = [widgets]

    def decorator(func):
        func._widget_info = {"widget_name": widgets, "event": event}
        return func

    return decorator


def connect_decorated_methods(widget: UiMainWindow):
    """Connects all decorated methods to their respective widgets"""
    for attr_name in dir(widget):
        connector_method = getattr(widget, attr_name)
        if hasattr(connector_method, "_widget_info"):
            widget_info = connector_method._widget_info
            for name in widget_info["widget_name"]:
                found_widget = getattr(widget, name)
                signal = getattr(found_widget, widget_info["event"])

                signal.connect(connector_method)
