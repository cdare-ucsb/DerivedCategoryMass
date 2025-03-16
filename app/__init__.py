from flask import Flask
from flask_socketio import SocketIO
import os

# Initialize SocketIO globally (without app yet)
socketio = SocketIO(cors_allowed_origins="*")  # CORS is required for cross-origin requests

def create_app(config_class=None):
    """
    Factory function to create and configure the Flask app.
    """
    app = Flask(__name__)

    # Load configuration settings
    if config_class:
        app.config.from_object(config_class)
    else:
        app.config.from_mapping(
            SECRET_KEY=os.environ.get('SECRET_KEY') or 'your-secret-key',
            MODEL_PATH='app/models/trained_model.pth'
        )

    # Import and register routes
    with app.app_context():
        from app import routes  # Ensure routes are registered
        app.register_blueprint(routes.bp)

    # Attach SocketIO to the app
    socketio.init_app(app, async_mode="threading")  

    return app
