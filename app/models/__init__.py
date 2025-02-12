import torch
from app.models.model import SimpleModel
import os

# Global variable to store the loaded model
model = None

def load_model(model_path=None):
    """
    Loads the PyTorch model from the specified path.
    If no path is provided, it uses the default from the config.
    """
    global model
    if model is None:
        model_path = model_path or os.environ.get('MODEL_PATH', 'app/models/trained_model.pth')
        print(f"Loading model from: {model_path}")
        
        # Initialize and load the model
        model = SimpleModel()  # Replace SimpleModel with your actual model class
        model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
        model.eval()  # Set model to evaluation mode
    return model


def preprocess_input(input_data):
    """
    Converts and preprocesses the input data to a PyTorch tensor.
    Adjust this function based on your model's expected input format.
    """
    return torch.tensor(input_data, dtype=torch.float32).unsqueeze(0)  # Example for single input


def predict(model, input_tensor):
    """
    Performs prediction using the given model and input tensor.
    """
    with torch.no_grad():
        output = model(input_tensor)
    return output
