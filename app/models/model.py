import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F  # Correct import for relu, softmax, etc.

class SimpleModel(nn.Module):
    def __init__(self):
        super(SimpleModel, self).__init__()
        self.fc = nn.Linear(10, 1)

    def forward(self, x):
        return self.fc(x)



class AdvancedModel(nn.Module):
    def __init__(self, input_size=10, hidden_sizes=[64, 128, 64], output_size=1, dropout_rate=0.2):
        """
        Initializes a fully connected neural network with multiple hidden layers.
        
        Args:
            input_size (int): Number of input features.
            hidden_sizes (list of int): Sizes of hidden layers.
            output_size (int): Number of output features.
            dropout_rate (float): Dropout rate for regularization.
        """
        super(AdvancedModel, self).__init__()

        self.hidden_layers = nn.ModuleList()  # List of hidden layers
        self.dropout = nn.Dropout(p=dropout_rate)  # Dropout for regularization
        
        # Create the hidden layers
        prev_size = input_size
        for hidden_size in hidden_sizes:
            self.hidden_layers.append(nn.Linear(prev_size, hidden_size))
            prev_size = hidden_size

        # Output layer
        self.output_layer = nn.Linear(prev_size, output_size)

    def forward(self, x):
        """
        Defines the forward pass of the model.
        
        Args:
            x (torch.Tensor): Input tensor.
        
        Returns:
            torch.Tensor: Output tensor.
        """
        for layer in self.hidden_layers:
            x = F.relu(layer(x))  # Apply ReLU activation after each hidden layer
            x = self.dropout(x)   # Apply dropout for regularization

        # Pass through the output layer without activation (e.g., for regression)
        x = self.output_layer(x)
        return x
