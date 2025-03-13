import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F  # Correct import for relu, softmax, etc.
import numpy as np
import plotly.graph_objects as go


from .SphericalTwist import SphericalTwist
from .CoherentSheaf import LineBundle


def generate_training_data_single_twist(line_bundle_1, line_bundle_2, degree=1):
    """
    Generate training data for the model.
    """

    _TRAIN_VECTOR_SIZE = 20000
    
    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                          LineBundle(line_bundle_2, catagory='K3'),
                          degree=degree)

    
    x_train = np.random.uniform(-10, 10, (_TRAIN_VECTOR_SIZE, 1))
    y_train = np.random.uniform(-10, 10, (_TRAIN_VECTOR_SIZE, 1))
    # Compute z_train properly as a NumPy array
    z_train = np.array([sph.mass(x_train[i, 0],
                                y_train[i, 0],
                                degree) for i in range(_TRAIN_VECTOR_SIZE)])

    # Reshape to ensure it has shape (1000, 1) instead of (1000,)
    z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

    # Convert to PyTorch tensors
    x_train_tensor = torch.tensor(x_train, dtype=torch.float32)
    y_train_tensor = torch.tensor(y_train, dtype=torch.float32)
    z_train_tensor = torch.tensor(z_train, dtype=torch.float32)

    return x_train_tensor, y_train_tensor, z_train_tensor


class FNN_model_1(nn.Module):
    def __init__(self):
        super(FNN_model_1, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(2, 64),  # Input: (x, y)
            nn.ReLU(),
            nn.Linear(64, 64),
            nn.ReLU(),
            nn.Linear(64, 1)  # Output: z
        )

    def forward(self, xy):
        return self.model(xy)
    
class FNN_model_2(nn.Module):
    def __init__(self):
        super(FNN_model_2, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(2, 128),  # Input: (x, y)
            nn.SiLU(),  # Swish activation (better than ReLU)
            nn.Linear(128, 128),
            nn.BatchNorm1d(128),  # Normalization for stable training
            nn.SiLU(),
            nn.Linear(128, 128),
            nn.Dropout(0.2),  # Prevent overfitting
            nn.SiLU(),
            nn.Linear(128, 1)  # Output: predicted z (mass)
        )

    def forward(self, xy):
        return self.model(xy)
    
class FNN_model_3(nn.Module):
    def __init__(self):
        super(FNN_model_3, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(2, 64),  # Input: (x, y)
            nn.SiLU(),
            nn.Linear(64, 64),
            nn.SiLU(),
            nn.Linear(64, 1)  # Output: z
        )

    def forward(self, xy):
        return self.model(xy)


if __name__ == '__main__':
    # Generate training data
    x_train_tensor, y_train_tensor, z_train_tensor = generate_training_data_single_twist(1, 2, degree=1)

    # Prepare input data as (x, y) pairs
    xy_train_tensor = torch.cat([x_train_tensor, y_train_tensor], dim=1)

    # Initialize model
    model = FNN_model_3()

    # Define loss function and optimizer
    criterion = nn.MSELoss()  # Mean Squared Error loss
    optimizer = optim.Adam(model.parameters(), lr=0.01)

    # Use learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=20, verbose=True)


    # Training loop
    epochs = 5000
    for epoch in range(epochs):
        optimizer.zero_grad()
        
        # Forward pass
        z_pred = model(xy_train_tensor)
        
        # Compute loss
        loss = criterion(z_pred, z_train_tensor)
        
        # Backward pass
        loss.backward()
        optimizer.step()

        # Step the scheduler
        scheduler.step(loss)
        
        # Print loss occasionally
        if epoch % 500 == 0:
            print(f'Epoch {epoch}, Loss: {loss.item():.6f}')


    # Save the trained model
    torch.save(model.state_dict(), 'trained_model.pth')



    # PLOT PREDICTED VALUES
    

    # Define x values (spread a  a region)
    x_vals = np.linspace(-5, 5, 150)  # X values from -5 to 5
    _NUM_Y_VALS_ = 100

    # Generate y values satisfying y > x^2
    y_vals = []
    for x in x_vals:
        y_range = np.linspace(0.1, 5, _NUM_Y_VALS_)  # 50 points per x value
        y_vals.append(y_range)

    # Convert to numpy array
    y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # Repeat x values to match the shape of y
    x_vals = np.repeat(x_vals, _NUM_Y_VALS_)  # Each x value repeats 10 times

    # Convert x and y values to tensors
    x_vals_tensor = torch.tensor(x_vals, dtype=torch.float32).reshape(-1, 1)
    y_vals_tensor = torch.tensor(y_vals, dtype=torch.float32).reshape(-1, 1)

    xy_vals_tensor = torch.cat([x_vals_tensor, y_vals_tensor], dim=1)

    # Set model to evaluation mode
    model.eval()

    # Disable gradient tracking for faster computation
    with torch.no_grad():
        z_vals_tensor = model(xy_vals_tensor)  # Get predictions

    z_vals = z_vals_tensor.numpy().flatten()  # Convert to NumPy and flatten

    # Plot the surface
    fig = go.Figure(data=[go.Scatter3d(z=z_vals, x=x_vals, y=y_vals,
                                    mode='markers', marker=dict(size=3, color=z_vals, colorscale='viridis'))])

    fig.update_layout(
        title="",
        autosize=True,
        margin=dict(l=0, r=0, b=0, t=30),
        scene=dict(
            
            bgcolor="white",  # Changes the 3D plot background ,

            xaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground =True,
                zerolinecolor="white",),
            yaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground =True,
                zerolinecolor="white"),
            zaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground =True,
                zerolinecolor="white"
            )
        )
    )

    fig.show()
