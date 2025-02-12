from flask import Flask
import torch  # Import PyTorch for later use
import os

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

    return app
