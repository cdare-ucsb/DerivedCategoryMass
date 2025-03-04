import webbrowser
from threading import Timer
import os
import importlib
import sys
from app import create_app


# Define required dependencies and their versions
REQUIRED_DEPENDENCIES = {
    "flask": "3.1.0",
    "torch": "2.2.2",
    "numpy": "1.25.1",
    "pandas": "2.2.3",
    "plotly": "5.13.1",
}


app = create_app()

# Automatically open the app in the default web browser
def open_browser():
    webbrowser.open_new('http://127.0.0.1:5000/')


def check_dependencies():
    """
    Helper function to check if the required dependencies are installed and up-to-date.
    If not all of the required dependencies are installed, the function will print a warning
    message and exit the program.
    """
    missing_dependencies = []
    for package, required_version in REQUIRED_DEPENDENCIES.items():
        try:
            module = importlib.import_module(package)
            installed_version = getattr(module, '__version__', None)
            if installed_version is None:
                print(f" \033[91m Warning: Could not determine version for {package}.\033[0m")
            elif installed_version != required_version:
                print(f"\033[91m Version mismatch: {package} (Installed: {installed_version}, Required: {required_version}) \033[0m")
                missing_dependencies.append(f"{package}=={required_version}")
        except ImportError:
            print(f"\033[91m Missing package: {package}\033[0m")
            missing_dependencies.append(f"{package}=={required_version}")

    if missing_dependencies:
        print("\nTo install the correct dependencies, run:\n")
        print(f"\t\033[1;92m pip install -r requirements.txt \033[0m")
        sys.exit(1)



if __name__ == "__main__":

    check_dependencies()
    
    # Only open the browser if not reloading (prevents opening twice in debug mode)
    if not app.debug or os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        # Use a timer to give the server time to start before opening the browser
        Timer(0.5, open_browser).start()

    app.run(debug=False)
