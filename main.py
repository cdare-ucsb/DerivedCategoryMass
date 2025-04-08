
import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QStackedWidget
)

from app.windows.home_page import HomePage
from app.windows.form_page import FormPage
from app.windows.plot_page import PlotPage
from app.windows.intersection_form_page import IntersectionFormPage

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Derived Category Mass App")

        self.stack = QStackedWidget()

        self.home_page = HomePage(self)
        self.form_page = FormPage(self)
        self.plot_page = PlotPage(self)
        self.intersection_form_page = IntersectionFormPage(self)


        self.stack.addWidget(self.home_page)
        self.stack.addWidget(self.form_page)
        self.stack.addWidget(self.plot_page)
        self.stack.addWidget(self.intersection_form_page)


        self.setCentralWidget(self.stack)
        self.resize(1000, 700)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())