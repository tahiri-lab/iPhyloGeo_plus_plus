import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel
from PyQt5.QtCore import QTimer

class LoadingAnimation(QWidget):
    def __init__(self):
        super().__init__()

        # Set up the layout
        self.layout = QVBoxLayout()

        # Create a QLabel for the loading animation
        self.loading_label = QLabel("Loading")
        self.layout.addWidget(self.loading_label)

        # Set the layout to the main window
        self.setLayout(self.layout)

        # Initialize the timer
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_loading_text)

        # Start the timer with an interval of 500 milliseconds
        self.timer.start(500)

        # Initialize the loading states
        self.loading_states = ["Loading.", "Loading..", "Loading..."]
        self.current_state_index = 0

    def update_loading_text(self):
        # Update the text of the QLabel
        self.loading_label.setText(self.loading_states[self.current_state_index])

        # Update the index to the next state
        self.current_state_index = (self.current_state_index + 1) % len(self.loading_states)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = LoadingAnimation()
    window.show()
    sys.exit(app.exec_())
