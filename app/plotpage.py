import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QLabel, QVBoxLayout,
    QComboBox, QPushButton, QStackedWidget, QLineEdit, QMessageBox, QHBoxLayout, QGraphicsOpacityEffect
)
from PyQt5.QtCore import Qt, QPropertyAnimation, QEasingCurve, QTimer
from PyQt5.QtWebEngineWidgets import QWebEngineView



import plotly.graph_objects as go
import plotly.io as pio



GEOMETRIC_CONTEXTS = ["Local P1", "Local P2", "P1", "P2", "K3 Surface"]

class HomePage(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent

        input_width = 300

        # ========= Outer Layout =========
        outer_layout = QVBoxLayout()
        outer_layout.setAlignment(Qt.AlignCenter)

        # ========= Title Label (Fades In First) =========
        self.welcome_label = QLabel("<h1>Welcome to derived category mass</h1>")
        self.welcome_label.setAlignment(Qt.AlignCenter)

        title_opacity = QGraphicsOpacityEffect()
        self.welcome_label.setGraphicsEffect(title_opacity)
        self.title_fade = QPropertyAnimation(title_opacity, b"opacity")
        self.title_fade.setDuration(800)
        self.title_fade.setStartValue(0.0)
        self.title_fade.setEndValue(1.0)
        self.title_fade.setEasingCurve(QEasingCurve.InOutQuad)

        # ========= Content Widget (Fades In Second) =========
        self.content_widget = QWidget()
        content_layout = QVBoxLayout()
        content_layout.setSpacing(20)
        content_layout.setAlignment(Qt.AlignCenter)

        # Geometry context label
        self.context_label = QLabel("Choose your geometric context:")
        self.context_label.setAlignment(Qt.AlignCenter)
        content_layout.addWidget(self.context_label)

        # Dropdown
        self.dropdown = QComboBox()
        self.dropdown.addItems(["Local P1", "Local P2", "P1", "P2", "K3 Surface"])
        self.dropdown.setMaximumWidth(input_width)

        dropdown_row = QHBoxLayout()
        dropdown_row.setAlignment(Qt.AlignCenter)
        dropdown_row.addWidget(self.dropdown)
        content_layout.addLayout(dropdown_row)

        # K3-specific fields
        self.k3_fields_container = QWidget()
        k3_layout = QVBoxLayout()
        k3_layout.setAlignment(Qt.AlignCenter)

        self.picard_input = QLineEdit()
        self.picard_input.setPlaceholderText("Picard Rank (e.g. 1)")
        self.picard_input.setMaximumWidth(input_width)

        self.polarization_input = QLineEdit()
        self.polarization_input.setPlaceholderText("Polarization (e.g. H)")
        self.polarization_input.setMaximumWidth(input_width)

        picard_row = QHBoxLayout()
        picard_row.setAlignment(Qt.AlignCenter)
        picard_row.addWidget(self.picard_input)

        polar_row = QHBoxLayout()
        polar_row.setAlignment(Qt.AlignCenter)
        polar_row.addWidget(self.polarization_input)

        k3_layout.addLayout(picard_row)
        k3_layout.addLayout(polar_row)

        self.k3_fields_container.setLayout(k3_layout)
        content_layout.addWidget(self.k3_fields_container)
        self.k3_fields_container.hide()

        # Continue button
        self.continue_button = QPushButton("Continue")
        self.continue_button.setMaximumWidth(input_width)

        button_row = QHBoxLayout()
        button_row.setAlignment(Qt.AlignCenter)
        button_row.addWidget(self.continue_button)
        content_layout.addLayout(button_row)

        # Set signals
        self.dropdown.currentTextChanged.connect(self._update_k3_fields_visibility)
        self.continue_button.clicked.connect(self.go_to_form)

        self.content_widget.setLayout(content_layout)

        # Content fade-in
        content_opacity = QGraphicsOpacityEffect()
        content_opacity.setOpacity(0.0)  # ❗ Make invisible to start
        self.content_widget.setGraphicsEffect(content_opacity)
        
        self.content_fade = QPropertyAnimation(content_opacity, b"opacity")
        self.content_fade.setDuration(800)
        self.content_fade.setStartValue(0.0)
        self.content_fade.setEndValue(1.0)
        self.content_fade.setEasingCurve(QEasingCurve.InOutQuad)

        # Add to outer layout
        outer_layout.addWidget(self.welcome_label)
        outer_layout.addWidget(self.content_widget)
        self.setLayout(outer_layout)

        # Start animations
        self.title_fade.start()

        # Delay content fade by 1000ms
        QTimer.singleShot(1000, self.content_fade.start)

    def _update_k3_fields_visibility(self, text):
        if text == "K3 Surface":
            self.k3_fields_container.show()
        else:
            self.k3_fields_container.hide()

    def go_to_form(self):
        selected_context = self.dropdown.currentText()

        # If K3 surface, validate fields
        if selected_context == "K3 Surface":
            picard = self.picard_input.text().strip()
            polarization = self.polarization_input.text().strip()

            print(picard)
            print(polarization)

            if not picard or not polarization:
                QMessageBox.warning(
                    self,
                    "Missing K3 Data",
                    "Please fill in both Picard rank and polarization before continuing."
                )
                return  # Do not proceed

        # Everything is valid — move forward
        self.parent.form_page.set_geometry_context(
            selected_context,
            picard_rank=self.picard_input.text(),
            polarization=self.polarization_input.text()
        )
        self.parent.stack.setCurrentWidget(self.parent.form_page)

class FormPage(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent

        # Data placeholders
        self.geometry_context = None
        self.picard_rank = None
        self.polarization = None

        # Layout
        self.layout = QVBoxLayout()

        self.info_label = QLabel("Geometry context not selected yet.")
        self.layout.addWidget(self.info_label)

        self.generate_button = QPushButton("Generate Sample Plot")
        self.generate_button.clicked.connect(self.generate_plot)
        self.layout.addWidget(self.generate_button)

        self.setLayout(self.layout)

    def set_geometry_context(self, context_name, picard_rank=None, polarization=None):
        self.geometry_context = context_name
        self.picard_rank = picard_rank
        self.polarization = polarization

        info = f"<b>Selected Geometry Context:</b> {context_name}"

        if context_name == "K3 Surface":
            info += f"<br><b>Picard Rank:</b> {picard_rank or '(unspecified)'}"
            info += f"<br><b>Polarization:</b> {polarization or '(unspecified)'}"

        self.info_label.setText(info)

    def generate_plot(self):
        # Generate a dummy plot (replace with real derived category logic)
        fig = go.Figure()
        fig.add_trace(go.Bar(x=["ch₀", "ch₁", "ch₂"], y=[1, 2, 3]))
        plot_html = pio.to_html(fig, full_html=False, include_plotlyjs='cdn')

        # Optional K3-specific info
        k3_info = ""
        if self.geometry_context == "K3 Surface":
            k3_info = f"""
            <p><b>Picard Rank:</b> {self.picard_rank or '(unspecified)'}<br>
               <b>Polarization:</b> {self.polarization or '(unspecified)'}</p>
            """

        # Wrap full HTML page
        full_html = f"""
        <html>
        <head>
          <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
          <script id="MathJax-script" async
            src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
          </script>
        </head>
        <body>
          <h2>{self.geometry_context} Plot</h2>
          {k3_info}
          <p>This is a placeholder for a derived object’s Chern character:</p>
          <p>$$\\mathrm{{ch}}(E) = r + c_1 + \\frac{{1}}{{2}}c_1^2 - c_2$$</p>
          {plot_html}
        </body>
        </html>
        """

        self.parent.plot_page.update_html(full_html)
        self.parent.stack.setCurrentWidget(self.parent.plot_page)


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


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Derived Category Mass App")

        self.stack = QStackedWidget()

        self.home_page = HomePage(self)
        self.form_page = FormPage(self)
        self.plot_page = PlotPage(self)

        self.stack.addWidget(self.home_page)
        self.stack.addWidget(self.form_page)
        self.stack.addWidget(self.plot_page)

        self.setCentralWidget(self.stack)
        self.resize(1000, 700)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
