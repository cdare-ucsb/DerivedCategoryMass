
from PyQt5.QtWidgets import (
    QWidget, QLabel, QVBoxLayout,
    QPushButton, QHBoxLayout
)

from PyQt5.QtCore import Qt


from app.widgets.intersection_matrix_widget import IntersectionMatrixWidget


class IntersectionFormPage(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.layout.addStretch(1)  # top space

        # === Vertical stack of label + matrix
        self.center_stack = QVBoxLayout()
        self.center_stack.setAlignment(Qt.AlignCenter)

        self.info_label = QLabel("Please enter the intersection numbers for your divisor basis:")
        self.info_label.setAlignment(Qt.AlignCenter)
        self.center_stack.addWidget(self.info_label)

        self.matrix_wrapper = QWidget()  # will get inserted later
        self.center_stack.addWidget(self.matrix_wrapper)

        # Wrap center stack inside a horizontal layout to center it fully
        center_block = QWidget()
        center_block.setLayout(self.center_stack)

        center_hbox = QHBoxLayout()
        center_hbox.setAlignment(Qt.AlignCenter)
        center_hbox.addWidget(center_block)

        self.layout.addLayout(center_hbox)
        self.layout.addStretch(1)  # bottom space

        # === Continue button
        self.continue_button = QPushButton("Continue")
        self.continue_button.setFixedWidth(300)

        button_row = QHBoxLayout()
        button_row.setAlignment(Qt.AlignCenter)
        button_row.addWidget(self.continue_button)

        self.layout.addLayout(button_row)

        self.continue_button.clicked.connect(self.submit)

    def set_basis(self, basis_list):
        if self.matrix_wrapper:
            self.center_stack.removeWidget(self.matrix_wrapper)
            self.matrix_wrapper.setParent(None)

        self.matrix_widget = IntersectionMatrixWidget(basis_list)

        self.matrix_wrapper = QWidget()
        wrapper_layout = QHBoxLayout()
        wrapper_layout.setAlignment(Qt.AlignCenter)
        wrapper_layout.addWidget(self.matrix_widget)
        self.matrix_wrapper.setLayout(wrapper_layout)

        self.center_stack.addWidget(self.matrix_wrapper)


    def submit(self):
        intersection_data = self.matrix_widget.get_intersection_form()
        # Pass the data back to form_page
        self.parent.form_page.set_k3_intersection_data(intersection_data)
        self.parent.stack.setCurrentWidget(self.parent.form_page)

