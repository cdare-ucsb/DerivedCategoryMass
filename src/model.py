import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import numpy as np
import plotly.graph_objects as go
import pandas as pd
import itertools
from tqdm import tqdm


from .SphericalTwist import SphericalTwist, DoubleSphericalTwist
from .CoherentSheaf import LineBundle
from .LocalP2 import LePotier


from dotenv import load_dotenv
import os


# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']
NN_MODEL_MODES = os.getenv("NN_MODEL_MODES").split(",") # ['mass', 'disc']





class SphericalTwistNeuralNetwork():
    r"""!
    This class implements the functionality required to train a neural network model to predict the mass or the discrete Laplacian of a single spherical twist,
    so that PyTorch does not directly need to be imported into the main application file for the Flask app. The SingleTwistModel class
    acts as a wrapper for a general PyTorch neural network model, and saves the training data as a member variabe. The class also provides
    methods to train the model, save the model, load the model, and plot the predictions of the model using Plotly.
    """

    def __init__(self, line_bundle_1,
                line_bundle_2,
                catagory, 
                line_bundle_3 = None,
                degree=1,
                x_min=-5, x_max=5, y_min=0, y_max=5,
                data_size=20000,
                mode='mass',
                x_granularity=0.1,
                y_granularity=0.05):
        r"""!
        The constructor for the SingleTwistModel class. The constructor initializes the model, the training data, and the model mode.
        By default, the model mode is set to 'mass', which means that the model is trained to predict the mass of the spherical twist.
        The model mode can also be set to 'disc', which means that the model is trained to predict the discrete Laplacian of the spherical twist. 
        The constructor also checks the validity of the input arguments and raises a ValueError if any of the input arguments are invalid

        \param line_bundle_1 The first line bundle in the spherical twist
        \param line_bundle_2 The second line bundle in the spherical twist
        \param catagory The catagory of the model; this is either 'P1', 'P2', or 'K3'
        \param degree The degree of the K3 surface; this is only used when the catagory is 'K3'
        \param x_min The minimum x value for the input data
        \param x_max The maximum x value for the input data
        \param y_min The minimum y value for the input data
        \param y_max The maximum y value for the input data
        \param data_size The number of data points to generate for the training data
        \param mode The mode of the model; this is either 'mass' for when the model should be predicting the mass, or 'disc' for when the model should be predicting the discrete Laplacian
        \param x_granularity The granularity of the x values for the discrete Laplacian; this is only used when the mode is 'disc'
        \param y_granularity The granularity of the y values for the discrete Laplacian; this is only used when the mode is 'disc'
        """
        
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f'Invalid catagory. Choose from {IMPLEMENTED_CATAGORIES}')
        if x_min >= x_max:
            raise ValueError('x_min must be less than x_max')
        if y_min >= y_max:
            raise ValueError('y_min must be less than y_max')
        if data_size <= 0:
            raise ValueError('data_size must be a positive integer')
        if not isinstance(degree, int):
            raise ValueError('degree must be an integer')
        if not isinstance(line_bundle_1, int):
            raise ValueError('line_bundle_1 must be an integer')
        if not isinstance(line_bundle_2, int):
            raise ValueError('line_bundle_2 must be an integer')
        if line_bundle_3 is not None and not isinstance(line_bundle_3, int):
            raise ValueError('line_bundle_3 must be an integer')
        if not isinstance(x_min, (int, float)):
            raise ValueError('x_min must be a number')
        if not isinstance(x_max, (int, float)):
            raise ValueError('x_max must be a number')
        if not isinstance(y_min, (int, float)):
            raise ValueError('y_min must be a number')
        if not isinstance(y_max, (int, float)):
            raise ValueError('y_max must be a number')
        if not isinstance(data_size, int):
            raise ValueError('data_size must be an integer')
        
        if catagory != 'P2' and y_min < 0:
            y_min = 0

        if mode not in NN_MODEL_MODES:
            raise ValueError(f'Invalid mode. Choose from {NN_MODEL_MODES}')
        
        self.mode = mode ## The mode of the model dictates what the model is trained to predict; this is either 'mass' for when the model should be predicting the mass, or 'disc' for when the model should be predicting the discrete Lagrangian

        self.catagory = catagory ## The catagory of the model; this is either 'P1', 'P2', or 'K3'
        
        self.degree = degree ## An optional argument for the degree of the K3 surface; this is only used when the catagory is 'K3'
    
        self.catagory = catagory ## The catagory of the model; this is either 'P1', 'P2', or 'K3'

        self.x_min = x_min ## The minimum x value for the input data

        self.x_max = x_max ## The maximum x value for the input data

        self.y_min = y_min ## The minimum y value for the input data
         
        self.y_max = y_max ## The maximum y value for the input data
    
        # Construct the xy-input data
        x_train, y_train = _generate_xy_grid(x_min, x_max, y_min, y_max, catagory=catagory, data_size=data_size, random_dist=True)
        
        x_tensor = torch.tensor(x_train, dtype=torch.float32).unsqueeze(1) # shape (N, 1)
        y_tensor = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1)

        
        sph = None
        if line_bundle_3 is not None:
            sph = DoubleSphericalTwist(LineBundle(line_bundle_1, catagory=catagory),
                            LineBundle(line_bundle_2, catagory=catagory),
                            LineBundle(line_bundle_3, catagory=catagory),
                            degree=degree)
        else:    
            sph = SphericalTwist(LineBundle(line_bundle_1, catagory=catagory),
                            LineBundle(line_bundle_2, catagory=catagory),
                            degree=degree)
        
            
        z_train = []
        
        if mode == 'mass':
            # Standard mode, create output tensor for mass prediction
            if catagory == 'K3':
                # Format input data for K3 case

                z_train = np.array([sph.mass(x_train[i],
                                            y_train[i],
                                            degree) for i in range(data_size)])
                z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

                self.output_tensor = torch.tensor(z_train, dtype=torch.float32) 

                self.input_tensor = torch.cat([x_tensor, y_tensor], dim=1) 
            elif catagory == 'P2':
                # Format input data for P2 case

                z_train = np.array([sph.mass(x_train[i], y_train[i]) for i in range(data_size)])
                z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

                self.output_tensor = torch.tensor(z_train, dtype=torch.float32)
                self.input_tensor = torch.cat([x_tensor, y_tensor], dim=1)
            elif catagory == 'P1':
                # Format input data for P1 case

                z_train = np.array([sph.mass(complex(x_train[i], y_train[i])) for i in range(data_size)])
                z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

                self.output_tensor = torch.tensor(z_train, dtype=torch.float32)
                self.input_tensor = torch.cat([x_tensor, y_tensor], dim=1)
            else:
                raise ValueError('Somehow an invalid catagory was passed to generate_training_data_single_twist; this should have been caught earlier')
        

        elif mode == 'disc':


            def _disc_Lapl(x, y):
                """
                Helper method to compute the discrete Laplacian of the spherical twist at a given point (x, y)
                """
                if self.catagory == 'P2':
                    return sph.mass(x + x_granularity, y) + sph.mass(x - x_granularity, y) + sph.mass(x, y + y_granularity) + sph.mass(x, y - y_granularity) - 4*sph.mass(x, y)
                elif self.catagory == 'P1':
                    return sph.mass(complex(x + x_granularity, y)) + sph.mass(complex(x - x_granularity, y)) + sph.mass(complex(x, y + y_granularity)) + sph.mass(complex(x, y - y_granularity)) - 4*sph.mass(complex(x, y))
                elif self.catagory == 'K3':
                    return sph.mass(x + x_granularity, y, self.degree) + sph.mass(x - x_granularity, y, self.degree) + sph.mass(x, y + y_granularity, self.degree) + sph.mass(x, y - y_granularity, self.degree) - 4*sph.mass(x, y, self.degree)
                else:
                    raise ValueError("Invalid catagory")

            # Discrete Laplacian mode, create output tensor for discrete Laplacian prediction
            z_train = np.array([_disc_Lapl(x_train[i], y_train[i]) for i in range(data_size)])
            z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

            self.output_tensor = torch.tensor(z_train, dtype=torch.float32)

            self.input_tensor = torch.cat([x_tensor, y_tensor], dim=1)
        else:
            raise ValueError('Somehow an invalid mode was passed to generate_training_data_single_twist; this should have been caught earlier')


        class FNN_model_3(nn.Module):
            r"""!
            The PyTorch neural network model for the SingleTwistModel class. The model is a simple feedforward neural network with 3 layers.
            """

            def __init__(self):
                r"""!
                Initialize the neural network model with 3 layers. The input layer has 2 nodes, the hidden layer has 64 nodes, and the output layer has 1 node.
                The activation function used is the SiLU activation function, which is the Sigmoid-weighted Linear Unit activation function.
                """
                super(FNN_model_3, self).__init__()
                self.model = nn.Sequential(
                    nn.Linear(2, 64),  # Input: (x, y)
                    nn.SiLU(),
                    nn.Linear(64, 64),
                    nn.SiLU(),
                    nn.Linear(64, 1)  # Output: z
                )

            def forward(self, xy):
                """
                The forward method of the neural network model. This method takes the input tensor xy and passes it through the neural network model.
                """
                return self.model(xy)

        self.model = FNN_model_3() ## The PyTorch neural network model for the SingleTwistModel class. The model is a simple feedforward neural network with 3 layers.

        



    def train(self, num_epochs):
        r"""!
        Method to train the neural network model for the SingleTwistModel class. The method uses the Adam optimizer and the Mean Squared Error loss function.
        The method also uses a learning rate scheduler to adjust the learning rate during training. The method prints the loss every 10% of the total number of epochs.

        \param num_epochs The number of epochs to train the model for
        """
        

        # Define loss function and optimizer
        criterion = nn.MSELoss()  # Mean Squared Error loss
        optimizer = optim.Adam(self.model.parameters(), lr=0.01)

        # Use learning rate scheduler
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                        factor=0.5, patience=20,
                                                        verbose=True)

        print_loss_num = int(num_epochs / 10)
        outer_progress = tqdm(range(num_epochs), desc="Training Progress", unit="epoch")

        for epoch in outer_progress:
            # Training loop
            
            optimizer.zero_grad()
            
            # Forward pass
            pred = self.model(self.input_tensor)
            
            # Compute loss
            loss = criterion(pred, self.output_tensor.view(-1, 1))
            
            # Backward pass
            loss.backward()
            optimizer.step()

            # Step the scheduler
            scheduler.step(loss)
            
            # Print loss occasionally
            if epoch % print_loss_num == 0:
                print(f'Epoch {epoch}, Loss: {loss.item():.6f}')

    def save_model(self, filename):
        r"""!
        Method to save the PyTorch neural network's state to a file specified by the filename argument. This simply calls the torch.save method
        with the model's state dictionary and the filename argument.

        \param filename The name of the file to save the model's state to
        """
        torch.save(self.model.state_dict(), filename)

    def load_model(self, filename):
        r"""!
        Method to load the PyTorch neural network's state from a file specified by the filename argument. This simply calls the torch.load method
        with the filename argument and sets the model's state dictionary to the loaded state dictionary.

        \param filename The name of the file to load the model's state from
        """
        self.model.load_state_dict(torch.load(filename))

    

    def predictions_to_plotly(self, color_scale='viridis'):
        r"""!
        Method to plot the predictions of the neural network model using Plotly. The method generates a grid of x and y values and passes them
        through the neural network model to get the predicted z values. The method then creates a Plotly 3D scatter plot of the predicted z values
        with the x and y values as the x and y coordinates. The method returns the Plotly figure object.

        \param color_scale The color scale to use for the Plotly 3D scatter plot; the default is 'viridis'

        \return The Plotly figure object of the 3D scatter plot of the predicted z values
        """

        
        x_vals, y_vals = _generate_xy_grid(self.x_min, self.x_max, self.y_min, self.y_max, catagory=self.catagory, data_size=20000, random_dist=False)

        # Convert x and y values to tensors
        x_vals_tensor = torch.tensor(x_vals, dtype=torch.float32).reshape(-1, 1)
        y_vals_tensor = torch.tensor(y_vals, dtype=torch.float32).reshape(-1, 1)

        xy_vals_tensor = torch.cat([x_vals_tensor, y_vals_tensor], dim=1)

        # Set model to evaluation mode
        self.model.eval()

        # Disable gradient tracking for faster computation
        with torch.no_grad():
            z_vals_tensor = self.model(xy_vals_tensor)  # Get predictions

        z_vals = z_vals_tensor.numpy().flatten()  # Convert to NumPy and flatten

        # Plot the surface
        fig = go.Figure(data=[go.Scatter3d(z=z_vals, x=x_vals, y=y_vals,
                                        mode='markers', marker=dict(size=3, color=z_vals, colorscale=color_scale))])

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

        return fig



            

class SingleTwistCollectionModel():
    r"""!
    This class implements the functionality required to train a neural network model to predict the mass or the discrete Laplacian of a collection of single spherical twists,
    so that PyTorch does not directly need to be imported into the main application file for the Flask app. The SingleTwistCollectionModel class
    acts as a wrapper for a general PyTorch neural network model, and saves the training data as a member variabe. The class also provides
    methods to train the model, save the model, load the model, and plot the predictions of the model using Plotly. Unlike the SingleTwistModel, this
    class does not actually create any of the data, but instead takes in a filename to a CSV file that contains the training data. The format of the CSV
    file should be as follows:

    x, y, line bundle 1, line bundle 2, (optional) degree, mass

    The  create_training_data_single_twist_collection method in the SphericalTwist module can be used to generate the training data and save it to a CSV file.
    """

    def __init__(self, filename, catagory):
        r"""!
        Creates the SingleTwistCollectionModel object with the specified filename and catagory. The constructor reads the CSV file and creates a PyTorch Dataset
        and DataLoader object from the data. The constructor also initializes the neural network model with the specified catagory.

        \param filename The name of the CSV file that contains the training data
        \param catagory The catagory of the model; this is either 'P1', 'P2', or 'K3'
        """

        class DataClass(Dataset):
            r"""!
            The PyTorch neural network model for the SingleTwistModel class. The model is a simple feedforward neural network with 3 layers.
            """

            def __init__(self, csv_file):
                self.csv_file = csv_file  # Just store the filename, NOT the data
                self.total_samples = sum(1 for _ in open(csv_file)) - 1  # Count lines (excluding header)

            def __len__(self):
                return self.total_samples  # Total number of rows

            def __getitem__(self, index):
                # Read only the required row using `skiprows`
                row = pd.read_csv(self.csv_file, skiprows=index+1, nrows=1).iloc[0]

                # Explicitly use `.iloc[]` to avoid FutureWarning
                features = torch.tensor(row.iloc[:-1].values, dtype=torch.float32)
                label = torch.tensor(row.iloc[-1], dtype=torch.float32)

                return features, label
            
        dataset = DataClass(filename)
        # Create DataLoader
        self.dataloader = DataLoader(dataset, batch_size=64, shuffle=True, num_workers=4, pin_memory=True) 

        class FNN_model_3(nn.Module):
            r"""!
            The PyTorch neural network model for the SingleTwistModel class. The model is a simple feedforward neural network with 3 layers.
            The input size of the model is determined by the catagory of the model.
            """
            def __init__(self, catagory):

                if catagory == 'K3':
                    self._input_size = 5
                elif catagory == 'P1':
                    self._input_size = 4
                elif catagory == 'P2':
                    self._input_size = 4
                else:
                    raise ValueError('Invalid catagory')
                
                super(FNN_model_3, self).__init__()
                self.model = nn.Sequential(
                    nn.Linear(self._input_size, 64),  # Input: (x, y)
                    nn.SiLU(),
                    nn.Linear(64, 64),
                    nn.SiLU(),
                    nn.Linear(64, 1)  # Output: z
                )

            def forward(self, xy):
                return self.model(xy)
            
        self.catagory = catagory ## The catagory of the model; this is either 'P1', 'P2', or 'K3'
        self.model = FNN_model_3(catagory) ## The PyTorch neural network model for the SingleTwistModel class. The model is a simple feedforward neural network with 3 layers.

    def train(self, num_epochs=50):
        r"""!
        Method to train the neural network model for the SingleTwistCollectionModel class. The method uses the Adam optimizer and the Mean Squared Error loss function.
        The method also uses a learning rate scheduler to adjust the learning rate during training. The method prints the loss every 10% of the total number of epochs.

        \param num_epochs The number of epochs to train the model for (default is 50)

        """
        
        #Define loss function and optimizer
        criterion = nn.MSELoss()  # Mean Squared Error loss
        optimizer = optim.Adam(self.model.parameters(), lr=0.01)

        # Use learning rate scheduler
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=20, verbose=True)
        print_rounds = int(num_epochs / 10)
        outer_progress = tqdm(range(num_epochs), desc="Training Progress", unit="epoch")

        for epoch in outer_progress:
            inner_progress = tqdm(self.dataloader, desc=f"Epoch {epoch+1}/{num_epochs}", unit="batch", leave=False)

            for batch in inner_progress:
                input_tensor, output_tensor = batch
                optimizer.zero_grad()
                
                # Forward pass
                pred = self.model(input_tensor)
                
                # Compute loss
                loss = criterion(pred, output_tensor.view(-1, 1))
                
                # Backward pass
                loss.backward()
                optimizer.step()

                # Step the scheduler
                scheduler.step(loss)
            
            # Print loss occasionally
            if epoch % print_rounds == 0:
                print(f'Epoch {epoch}, Loss: {loss.item():.6f}')




 
def create_training_data_single_twist_collection(filename, catagory, degree=1,
                                        x_min=-10, x_max=10, y_min=0, y_max=10,
                                        lb_1_min=-2, lb_1_max=2, lb_2_min=-3, lb_2_max=3,
                                        data_size=5000):
    r"""!
    Method to construct a CSV file with training data for a collection of single spherical twists. The method generates data for all possible combinations
    of line bundles lb_1 and lb_2 in the specified ranges. The method generates data for the specified catagory and degree, and saves the data to the specified
    filename. The method uses the _append_training_data_to_CSV_single_twist method to append the training data to the CSV file.

    WARNING: Every iteration of append_training_data_to_CSV_single_twist with the default 20,000 samples
    makes a file that is roughly 1.5MB and takes about 1.83s (on 8 Intel I9 cores). Thus, ranging over
    all spherical twists Tw_a O(b)  with -5≤ a ≤ 5 and -5 ≤ b ≤ 5 will produce a file that is roughly 150MB and
    will take roughly 3 and a half minutes (on 8 Intel I9 cores).


    \param filename The name of the CSV file to save the training data to
    \param catagory The catagory of the model; this is either 'P1', 'P2', or 'K3'
    \param degree The degree of the K3 surface; this is only used when the catagory is 'K3'
    \param x_min The minimum x value for the input data
    \param x_max The maximum x value for the input data 
    \param y_min The minimum y value for the input data
    \param y_max The maximum y value for the input data
    \param lb_1_min The minimum value for the first line bundle in the spherical twist
    \param lb_1_max The maximum value for the first line bundle in the spherical twist
    \param lb_2_min The minimum value for the second line bundle in the spherical twist
    \param lb_2_max The maximum value for the second line bundle in the spherical twist
    \param data_size The number of data points to generate for the training data
    
    """

    if catagory not in IMPLEMENTED_CATAGORIES:
        raise ValueError(f'Invalid catagory. Choose from {IMPLEMENTED_CATAGORIES}')
    if x_min >= x_max:
        raise ValueError('x_min must be less than x_max')
    if y_min >= y_max:
        raise ValueError('y_min must be less than y_max')
    if data_size <= 0:
        raise ValueError('data_size must be a positive integer')
    if not isinstance(degree, int):
        raise ValueError('degree must be an integer')
    if not isinstance(lb_1_min, int):
        raise ValueError('lb_1_min must be an integer')
    if not isinstance(lb_1_max, int):
        raise ValueError('lb_1_max must be an integer')
    if not isinstance(lb_2_min, int):
        raise ValueError('lb_2_min must be an integer')
    if not isinstance(lb_2_max, int):
        raise ValueError('lb_2_max must be an integer')
    if not isinstance(x_min, (int, float)):
        raise ValueError('x_min must be a number')
    if not isinstance(x_max, (int, float)):
        raise ValueError('x_max must be a number')
    if not isinstance(y_min, (int, float)):
        raise ValueError('y_min must be a number')
    if not isinstance(y_max, (int, float)):
        raise ValueError('y_max must be a number')
    if not isinstance(data_size, int):
        raise ValueError('data_size must be an integer')
    
    if catagory != 'P2' and y_min < 0:
        y_min = 0

    lb_1_range = range(lb_1_min, lb_1_max+1)
    lb_2_range = range(lb_2_min, lb_2_max+1)

    total_iterations = len(lb_1_range) * len(lb_2_range)  # Total combinations

    for lb_1, lb_2 in tqdm(itertools.product(lb_1_range, lb_2_range), total=total_iterations, desc="Generating Data", unit="pair"):
        _append_training_data_to_CSV_single_twist(filename, line_bundle_1=lb_1, line_bundle_2=lb_2,
                                                catagory=catagory, degree=degree, 
                                                x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
                                                data_size=data_size)








# class FNN_model_1(nn.Module):
#     def __init__(self, catagory):

#         if catagory == 'K3':
#             self._input_size = 5
#         elif catagory == 'P1':
#             self._input_size = 4
#         elif catagory == 'P2':
#             self._input_size = 4
#         else:
#             raise ValueError('Invalid catagory')

#         super(FNN_model_1, self).__init__()
#         self.model = nn.Sequential(
#             nn.Linear(self._input_size, 64),  # Input: (x, y)
#             nn.ReLU(),
#             nn.Linear(64, 64),
#             nn.ReLU(),
#             nn.Linear(64, 1)  # Output: z
#         )

#     def forward(self, xy):
#         return self.model(xy)
    
# class FNN_model_2(nn.Module):
#     def __init__(self, catagory):

#         if catagory == 'K3':
#             self._input_size = 5
#         elif catagory == 'P1':
#             self._input_size = 4
#         elif catagory == 'P2':
#             self._input_size = 4
#         else:
#             raise ValueError('Invalid catagory')
        
#         super(FNN_model_2, self).__init__()
#         self.model = nn.Sequential(
#             nn.Linear(self._input_size, 128),  # Input: (x, y)
#             nn.SiLU(),  # Swish activation (better than ReLU)
#             nn.Linear(128, 128),
#             nn.BatchNorm1d(128),  # Normalization for stable training
#             nn.SiLU(),
#             nn.Linear(128, 128),
#             nn.Dropout(0.2),  # Prevent overfitting
#             nn.SiLU(),
#             nn.Linear(128, 1)  # Output: predicted z (mass)
#         )

#     def forward(self, xy):
#         return self.model(xy)
    

    
# class FNN_model_3(nn.Module):
#     def __init__(self, catagory):

#         if catagory == 'K3':
#             self._input_size = 5
#         elif catagory == 'P1':
#             self._input_size = 4
#         elif catagory == 'P2':
#             self._input_size = 4
#         else:
#             raise ValueError('Invalid catagory')
        
#         super(FNN_model_3, self).__init__()
#         self.model = nn.Sequential(
#             nn.Linear(self._input_size, 64),  # Input: (x, y)
#             nn.SiLU(),
#             nn.Linear(64, 64),
#             nn.SiLU(),
#             nn.Linear(64, 1)  # Output: z
#         )

#     def forward(self, xy):
#         return self.model(xy)
    



###############################################
#         STATIC HELPER FUNCTIONS             #
###############################################




def _generate_xy_grid(x_min, x_max, y_min, y_max, catagory, random_dist = True, data_size=20000):
    r"""!
    Helper function to create a grid of x and y values for the input data. The function generates the x and y values
    either randomly or from a linear grid, depending on the value of the random_dist argument. The function also checks
    the validity of the input arguments and raises a ValueError if any of the input arguments are invalid

    \param x_min The minimum x value for the input data
    \param x_max The maximum x value for the input data
    \param y_min The minimum y value for the input data
    \param y_max The maximum y value for the input data
    \param catagory The catagory of the model; this is either 'P1', 'P2', or 'K3'
    \param random_dist A boolean flag to determine whether the x and y values should be generated randomly or from a linear grid
    \param data_size The number of data points to generate for the training data

    \return x_vals The x values for the input data
    \return y_vals The y values for the input data
    """
    if catagory not in IMPLEMENTED_CATAGORIES:
        raise ValueError(f'Invalid catagory. Choose from {IMPLEMENTED_CATAGORIES}')
    if x_min >= x_max:
        raise ValueError('x_min must be less than x_max')
    if y_min >= y_max:
        raise ValueError('y_min must be less than y_max')
    if data_size <= 0:
        raise ValueError('data_size must be a positive integer')
    

    if catagory != 'P2':
        # Can create data from random points on the upper half plane.

        if random_dist:
            # Randomly create x and y values from uniform distribution
            x_vals = np.random.uniform(x_min, x_max, size=(data_size, ))
            y_vals = np.random.uniform(y_min, y_max, size=(data_size, ))

            return x_vals, y_vals
        else:
            # The x and y values are not random, but are generated from a grid that should
            # be linearly spaced
            N = int(np.sqrt(data_size))
            y_vals = []
            x_vals = np.linspace(x_min, x_max, N)

            for x in x_vals:
                y_range = np.linspace(y_min, y_max, N)  
                y_vals.append(y_range)

            # Convert to numpy array
            y_vals = np.array(y_vals).flatten()  # Flatten the y array

            # Repeat x values to match the shape of y
            x_vals = np.repeat(x_vals, N)  # Each x value repeats N times

            return x_vals, y_vals

    else:


        DLP = LePotier()

        if random_dist:

            largest_abs_x = max(abs(x_min), abs(x_max))
            y_max = max(y_max, DLP.curve_estimate(largest_abs_x))
            

            x_vals, y_vals = [], []
            while len(x_vals) < data_size:
                x = np.random.uniform(x_min, x_max)
                y = np.random.uniform(-1, y_max)

                if DLP.is_above_curve(x, y):
                    x_vals.append(x)
                    y_vals.append(y)

            return x_vals, y_vals
        else:
            # The x and y values are not random, but are generated from a grid that should
            # be linearly spaced

            largest_abs_x = max(abs(x_min), abs(x_max))
            y_max = max(y_max, DLP.curve_estimate(largest_abs_x)) # Ensure that the values are large enough to cover the curve

            N = int(np.sqrt(data_size))
            y_vals = []
            x_vals = np.linspace(x_min, x_max, N)

            for x in x_vals:
                y_start = max(DLP.curve_estimate(x), y_min)

                y_range = np.linspace(y_start, y_max, N)  
                y_vals.append(y_range)

            # Convert to numpy array
            y_vals = np.array(y_vals).flatten()

            # Repeat x values to match the shape of y
            x_vals = np.repeat(x_vals, N)
            
            return x_vals, y_vals



def _append_training_data_to_CSV_single_twist(filename, line_bundle_1, line_bundle_2, catagory, degree=1,
                                            x_min=-5, x_max=5, y_min=0, y_max=10,
                                            data_size=20000):
    r"""!
    Helper function to append training data for a single spherical twist to a CSV file. The function generates the input data
    and the output data for the spherical twist and appends the data to the CSV file. The function checks the validity of the input
    arguments and raises a ValueError if any of the input arguments are invalid. The function uses the SphericalTwist class to compute
    the mass of the spherical twist at each point in the grid specified by the input data.

    \param filename The name of the CSV file to save the training data to
    \param line_bundle_1 The first line bundle in the spherical twist
    \param line_bundle_2 The second line bundle in the spherical twist  
    \param catagory The catagory of the model; this is either 'P1', 'P2', or 'K3'
    \param degree The degree of the K3 surface; this is only used when the catagory is 'K3'
    \param x_min The minimum x value for the input data
    \param x_max The maximum x value for the input data
    \param y_min The minimum y value for the input data
    \param y_max The maximum y value for the input data
    \param data_size The number of data points to generate for the training data    
    """
    
    if catagory not in IMPLEMENTED_CATAGORIES:
        raise ValueError(f'Invalid catagory. Choose from {IMPLEMENTED_CATAGORIES}')
    
    # Construct the xy-input data
    x_train, y_train = _generate_xy_grid(x_min, x_max, y_min, y_max, catagory=catagory, data_size=data_size, random_dist=True)

    # Constant columns for the line bundle data
    line_bundle_1_col = np.full((len(x_train), 1), line_bundle_1, dtype=np.int32)
    line_bundle_2_col = np.full((len(x_train), 1), line_bundle_2, dtype=np.int32)

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory=catagory),
                          LineBundle(line_bundle_2, catagory=catagory),
                          degree=degree)
    
    if catagory == 'K3':
        # Format input data for K3 case

        z_train = np.array([sph.mass(x_train[i, 0],
                                    y_train[i, 0],
                                    degree) for i in range(data_size)])
        z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)
        
        degree_col = np.full((len(x_train), 1), degree, dtype=np.int32)

        # Sample DataFrame
        df = pd.DataFrame({
            'x': x_train.flatten(),
            'y': y_train.flatten(),
            'LB1': line_bundle_1_col.flatten(),
            'LB2': line_bundle_2_col.flatten(),
            'Degree': degree_col.flatten(),
            'Mass': z_train.flatten()
        })

        file_exists = os.path.isfile(filename)
        # Append DataFrame to CSV (or create a new one)
        df.to_csv(filename, mode='a', header=not file_exists, index=False)  # Append data to CSV

        
    elif catagory == 'P1':
        # Format input data for P1 case

        z_train = np.array([sph.mass(complex(x_train[i, 0], y_train[i, 0])) for i in range(data_size)]),
        z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

        # Sample DataFrame
        df = pd.DataFrame({
            'x': x_train.flatten(),
            'y': y_train.flatten(),
            'LB1': line_bundle_1_col.flatten(),
            'LB2': line_bundle_2_col.flatten(),
            'Mass': z_train.flatten()
        })

        file_exists = os.path.isfile(filename)
        # Append DataFrame to CSV (or create a new one)
        df.to_csv(filename, mode='a', header=not file_exists, index=False)  # Append data to CSV

    elif catagory == 'P2':
        # Format input data for P2 case

        z_train = np.array([sph.mass(x_train[i, 0], y_train[i, 0]) for i in range(data_size)]),
        z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)
        
        # Sample DataFrame
        df = pd.DataFrame({
            'x': x_train.flatten(),
            'y': y_train.flatten(),
            'LB1': line_bundle_1_col.flatten(),
            'LB2': line_bundle_2_col.flatten(),
            'Mass': z_train.flatten()
        })

        file_exists = os.path.isfile(filename)
        # Append DataFrame to CSV (or create a new one)
        df.to_csv(filename, mode='a', header=not file_exists, index=False)  # Append data to CSV

        
    else:
        raise ValueError('Somehow an invalid catagory was passed to generate_training_data_single_twist; this should have been caught earlier')




###############################################
#                MAIN SCRIPT                  #
###############################################

if __name__ == '__main__':


    lb_1_range = range(-1, 15+1)
    lb_2_range = range(-15, 15+1)

    for lb1, lb2 in itertools.product(lb_1_range, lb_2_range):

        print("\n\n\n\n---------------------------------------------")
        print(f'Traning model for discontinuities of Tw_{lb1} O({lb2})')
        print("---------------------------------------------\n\n")

        st_model = SphericalTwistNeuralNetwork(line_bundle_1=lb1,
                                    line_bundle_2=lb2,
                                    catagory='K3',
                                    data_size=25000,
                                    mode='disc',
                                    x_min=-10,
                                    x_max=10,
                                    y_min=0,
                                    y_max=10)
        
        st_model.train(num_epochs=5000)
        
        filename =  os.path.abspath(os.path.join(os.path.dirname(__file__), f'../app/models/disc_K3_deg{1}_{lb1}_{lb2}.pth'))
        print(f"Saving model at {filename} .............")
        st_model.save_model(filename)
        print("Model saved successfully!")
        print("\n\n\n\n---------------------------------------------")

    



    # Generate training data
    # filename =  os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/training_data_K3.csv'))
    # create_K3_training_data_single_twist(filename)
    

    
    

   