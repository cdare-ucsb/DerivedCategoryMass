
from PyQt5.QtWidgets import (
    QWidget, QLabel, 
    QLineEdit, 
    QGridLayout, QSizePolicy
)
from PyQt5.QtCore import Qt 



class IntersectionMatrixWidget(QWidget):
    def __init__(self, basis, parent=None):
        super().__init__(parent)
        self.basis = basis
        self.entries = {}  # Maps (i,j) -> QLineEdit

        self._build_ui()

    def _build_ui(self):
        layout = QGridLayout()
        layout.setSpacing(8)

        n = len(self.basis)

        # Add top header
        for j, name in enumerate(self.basis):
            header = QLabel(f"<b>{name}</b>")
            header.setAlignment(Qt.AlignCenter)
            layout.addWidget(header, 0, j + 1)

        # Add side header and input fields
        for i, row_label in enumerate(self.basis):
            header = QLabel(f"<b>{row_label}</b>")
            header.setAlignment(Qt.AlignCenter)
            layout.addWidget(header, i + 1, 0)

            for j, col_label in enumerate(self.basis):
                # Create field only for upper triangle and diagonal
                if j >= i:
                    field = QLineEdit()
                    field.setAlignment(Qt.AlignCenter)
                    field.setFixedWidth(60)
                    layout.addWidget(field, i + 1, j + 1)
                    self.entries[(i, j)] = field

                    # If not diagonal, mirror value to symmetric cell
                    if j != i:
                        mirror_field = QLineEdit()
                        mirror_field.setReadOnly(True)
                        mirror_field.setAlignment(Qt.AlignCenter)
                        mirror_field.setFixedWidth(60)
                        layout.addWidget(mirror_field, j + 1, i + 1)
                        self.entries[(j, i)] = mirror_field

                        # Connect to keep symmetry
                        def make_sync_func(i1, j1):
                            return lambda text: self.entries[(j1, i1)].setText(text)

                        field.textChanged.connect(make_sync_func(i, j))

        self.setLayout(layout)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

    def get_intersection_form(self):
        """
        Returns the intersection form as a dictionary mapping
        (Div1, Div2) -> int, based on current inputs.
        """
        result = {}
        n = len(self.basis)
        for i in range(n):
            for j in range(n):
                field = self.entries.get((i, j))
                if field:
                    text = field.text().strip()
                    try:
                        val = int(text)
                    except ValueError:
                        val = 0  # Default or raise if you prefer
                    result[(self.basis[i], self.basis[j])] = val
        return result
