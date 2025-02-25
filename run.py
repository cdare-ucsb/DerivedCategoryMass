import webbrowser
from threading import Timer
import os
from app import create_app

app = create_app()

# Automatically open the app in the default web browser
def open_browser():
    webbrowser.open_new('http://127.0.0.1:5000/')

if __name__ == "__main__":
    
    # Only open the browser if not reloading (prevents opening twice in debug mode)
    if not app.debug or os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        # Use a timer to give the server time to start before opening the browser
        Timer(0.5, open_browser).start()

    app.run(debug=False)
