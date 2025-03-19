#!/bin/bash

VENV_DIR="venv"

# Create virtual environment if not present
if [ ! -d "$VENV_DIR" ]; then
    echo "🔧 Creating virtual environment..."
    python3 -m venv "$VENV_DIR"
fi

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Upgrade pip/setuptools/wheel
pip install --upgrade pip setuptools wheel

# Function to check if package is installed, install if not
function check_and_install {
    PACKAGE=$1
    if ! pip show "$PACKAGE" >/dev/null 2>&1; then
        echo -e "\033[33m📦 $PACKAGE not found, installing...\033[0m"
        pip install "$PACKAGE"
    else
        echo -e "\033[32m✅ $PACKAGE already installed.\033[0m"
    fi
}

# List your dependencies here:
check_and_install "Flask"
check_and_install "Flask-SocketIO"
check_and_install "torch"
check_and_install "numpy"
check_and_install "pandas"
check_and_install "plotly"

echo -e "\033[32m\n\n\n\n\t\tStarting Flask server...\033[0m\n\n\n\n"
python driver.py

deactivate