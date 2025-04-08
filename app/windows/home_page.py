
from PyQt5.QtWidgets import (
    QWidget, QLabel, QVBoxLayout, QComboBox, QPushButton,  QLineEdit, QMessageBox, QHBoxLayout, QGraphicsOpacityEffect
)
from PyQt5.QtCore import Qt, QPropertyAnimation, QEasingCurve, QTimer


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
        self.welcome_label = QLabel("<h1>Welcome to DerivedCategoryMass</h1>")
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
        self.dropdown.addItems(GEOMETRIC_CONTEXTS)
        self.dropdown.setMaximumWidth(input_width)

        dropdown_row = QHBoxLayout()
        dropdown_row.setAlignment(Qt.AlignCenter)
        dropdown_row.addWidget(self.dropdown)
        content_layout.addLayout(dropdown_row)

        # K3-specific fields
        self.k3_fields_container = QWidget()
        k3_layout = QVBoxLayout()
        k3_layout.setAlignment(Qt.AlignCenter)

        self.picard_basis_input = QLineEdit()
        self.picard_basis_input.setPlaceholderText("Enter a basis for the Picard Group (e.g. H, C, D)")
        self.picard_basis_input.setMaximumWidth(input_width)

        self.polarization_input = QLineEdit()
        self.polarization_input.setPlaceholderText("What is the choice of polarization? (e.g. H)")
        self.polarization_input.setMaximumWidth(input_width)

        picard_row = QHBoxLayout()
        picard_row.setAlignment(Qt.AlignCenter)
        picard_row.addWidget(self.picard_basis_input)

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

        if selected_context == "K3 Surface":
            picard = self.picard_basis_input.text().strip()
            polarization = self.polarization_input.text().strip()

            print(polarization)

            if not picard or not polarization:
                QMessageBox.warning(
                    self,
                    "Missing K3 Data",
                    "Please fill in both Picard rank and polarization before continuing."
                )
                return

            # Parse the polarization input as a list of basis divisors
            # E.g., "H, C, D" -> ["H", "C", "D"]
            basis = [s.strip() for s in picard.split(",") if s.strip()]
            print(basis)


            if len(basis) == 0:
                QMessageBox.warning(
                    self,
                    "Invalid Basis",
                    "Polarization input must contain at least one divisor label."
                )
                return

            self.parent.intersection_form_page.set_basis(basis)
            self.parent.form_page.set_geometry_context(
                selected_context,
                picard_rank=self.picard_basis_input.text(),
                polarization=self.polarization_input.text()
            )

            self.parent.stack.setCurrentWidget(self.parent.intersection_form_page)
        else:
            self.parent.form_page.set_geometry_context(
                selected_context,
                picard_rank=self.picard_basis_input.text(),
                polarization=self.polarization_input.text()
            )
            self.parent.stack.setCurrentWidget(self.parent.form_page)
