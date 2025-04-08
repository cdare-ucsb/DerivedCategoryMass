
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton
)

from PyQt5.QtWebEngineWidgets import QWebEngineView



class PlotPage(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        layout = QVBoxLayout()

        self.view = QWebEngineView()
        self.back_button = QPushButton("Back to Home")
        self.back_button.clicked.connect(lambda: self.parent.stack.setCurrentWidget(self.parent.home_page))

        layout.addWidget(self.view)
        layout.addWidget(self.back_button)
        self.setLayout(layout)

    def update_html(self, html):
        self.view.setHtml(html)

