#!/bin/bash

VENV_DIR="venv"

# Desired Python Version
MAX_SUPPORTED_PYTHON_MAJOR=3
MAX_SUPPORTED_PYTHON_MINOR=12

# Get current python version
PYTHON_VERSION=$(python3 --version | awk '{print $2}')
PYTHON_MAJOR=$(echo "$PYTHON_VERSION" | cut -d'.' -f1)
PYTHON_MINOR=$(echo "$PYTHON_VERSION" | cut -d'.' -f2)

# Check Python version
if  { [ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -gt 12 ]; }; then
    echo -e "\033[31mYour Python 3 version ($PYTHON_MAJOR.$PYTHON_MINOR) is too new.\033[0m"
    echo -e "\033[33mPlease downgrade to Python 3.12 using:\033[0m"
    echo -e "  brew uninstall python@3.13"
    echo -e "  brew install python@3.12"
    echo -e "Then recreate your virtual environment:"
    echo -e "  rm -rf venv"
    echo -e "  python3.12 -m venv venv"
    exit 1
fi

# Create virtual environment if not present
if [ ! -d "$VENV_DIR" ]; then
    echo "🔧 Creating virtual environment..."
    python3 -m venv "$VENV_DIR"
fi

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Upgrade pip/setuptools/wheel
pip install --upgrade pip setuptools wheel

# Function to install exact version of numpy
install_numpy() {
    INSTALLED_NUMPY=$(pip show numpy 2>/dev/null | grep Version | awk '{print $2}')
    if [[ "$INSTALLED_NUMPY" != "1.25.1" ]]; then
        echo -e "\033[33m📦 Installing numpy==1.25.1...\033[0m"
        pip uninstall numpy -y >/dev/null 2>&1
        pip install numpy==1.25.1
    else
        echo -e "\033[32m✅ numpy==1.25.1 already installed.\033[0m"
    fi
}

# Function to install exact version of plotly
install_plotly() {
    INSTALLED_plotly=$(pip show plotly 2>/dev/null | grep Version | awk '{print $2}')
    if [[ "$INSTALLED_NUMPY" != "5.13.1" ]]; then
        echo -e "\033[33m📦 Installing plotly==5.13.1...\033[0m"
        pip uninstall plotly -y >/dev/null 2>&1
        pip install plotly==5.13.1
    else
        echo -e "\033[32m✅ plotly==5.13.1 already installed.\033[0m"
    fi
}

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
check_and_install "pandas"
check_and_install "python-dotenv"
check_and_install "matplotlib"
check_and_install "tqdm"

install_numpy
install_plotly

echo -e "\033[32m\n\n\n\n\t\tStarting Flask server...\033[0m\n\n\n\n"
python driver.py

deactivate