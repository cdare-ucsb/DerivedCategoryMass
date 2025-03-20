#!/bin/bash

VENV_DIR="venv"

# Required Python version
REQUIRED_PYTHON="3.11"

# Check if python3.11 is installed
if ! command -v python3.11 &> /dev/null; then
    echo -e "\033[31m❌ Python 3.11 is NOT installed.\033[0m"
    echo -e "\033[33m➡️  Please install Python 3.11 first:\033[0m"
    echo -e "   🔹 On macOS (Homebrew): \033[32mbrew install python@3.11\033[0m"
    echo -e "   🔹 On Ubuntu/Debian (APT): \033[32msudo add-apt-repository ppa:deadsnakes/ppa && sudo apt update && sudo apt install python3.11\033[0m"
    exit 1
fi

# Check the installed version to ensure it's correct
PYTHON_VERSION=$(python3.11 --version 2>&1 | awk '{print $2}')
PYTHON_MAJOR=$(echo "$PYTHON_VERSION" | cut -d'.' -f1)
PYTHON_MINOR=$(echo "$PYTHON_VERSION" | cut -d'.' -f2)

if [ "$PYTHON_MAJOR" -ne 3 ] || [ "$PYTHON_MINOR" -ne 11 ]; then
    echo -e "\033[31m❌ Detected Python version: $PYTHON_VERSION\033[0m"
    echo -e "\033[33m➡️  Please ensure Python 3.11 is installed and set correctly.\033[0m"
    exit 1
fi

echo -e "\033[32m✅ Python 3.11 is correctly installed.\033[0m"

# Create virtual environment if not present
if [ ! -d "$VENV_DIR" ]; then
    echo "🔧 Creating virtual environment..."
    python3.11 -m venv "$VENV_DIR"
fi

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Ensure essential system dependencies are installed
echo -e "\033[33m🔧 Installing system dependencies for NumPy...\033[0m"
pip install --upgrade pip setuptools wheel cython

# Function to install exact version of numpy
install_numpy() {
    INSTALLED_NUMPY=$(pip show numpy 2>/dev/null | grep Version | awk '{print $2}')
    if [[ "$INSTALLED_NUMPY" != "1.25.1" ]]; then
        echo -e "\033[33m📦 Installing numpy==1.25.1 with prebuilt wheel...\033[0m"
        pip uninstall numpy -y >/dev/null 2>&1
        pip install --no-cache-dir --no-build-isolation numpy==1.25.1
    else
        echo -e "\033[32m✅ numpy==1.25.1 already installed.\033[0m"
    fi
}

# Function to install exact version of plotly (Fixed typo)
install_plotly() {
    INSTALLED_PLOTLY=$(pip show plotly 2>/dev/null | grep Version | awk '{print $2}')
    if [[ "$INSTALLED_PLOTLY" != "5.13.1" ]]; then
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
