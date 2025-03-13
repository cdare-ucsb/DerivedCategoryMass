import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F  # Correct import for relu, softmax, etc.
from torch.utils.data import DataLoader, Dataset
import numpy as np
import plotly.graph_objects as go
import pandas as pd
import itertools
from tqdm import tqdm


from SphericalTwist import SphericalTwist
from CoherentSheaf import LineBundle
from LocalP2 import LePotier


from dotenv import load_dotenv
import os


# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']





class SingleTwistModel():

    def __init__(self, line_bundle_1,
                line_bundle_2,
                catagory, 
                degree=1,
                x_min=-5, x_max=5, y_min=0, y_max=5,
                  data_size=20000):
        
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
        
        
    
        # Construct the xy-input data
        x_train, y_train = [], []
        if catagory == 'P2':
            x_train, y_train = _generate_xy_tensors_P2(x_min, x_max, y_min, y_max, data_size)
        else:
            x_train = np.random.uniform(x_min, x_max, (data_size, 1))
            y_train = np.random.uniform(y_min, y_max, (data_size, 1))
        
        x_tensor = torch.tensor(x_train, dtype=torch.float32)
        y_tensor = torch.tensor(y_train, dtype=torch.float32)

        

        sph = SphericalTwist(LineBundle(line_bundle_1, catagory=catagory),
                          LineBundle(line_bundle_2, catagory=catagory),
                          degree=degree)
        
        if catagory == 'K3':
        # Format input data for K3 case

            z_train = np.array([sph.mass(x_train[i, 0],
                                        y_train[i, 0],
                                        degree) for i in range(data_size)])
            z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

            self.output_tensor = torch.tensor(z_train, dtype=torch.float32)
            self.input_tensor = torch.cat([x_tensor, y_tensor], dim=1)
        elif catagory == 'P2':
            # Format input data for P2 case

            z_train = np.array([sph.mass(x_train[i, 0], y_train[i, 0]) for i in range(data_size)])
            z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

            self.output_tensor = torch.tensor(z_train, dtype=torch.float32)
            self.input_tensor = torch.cat([x_tensor, y_tensor], dim=1)
        elif catagory == 'P1':
            # Format input data for P1 case

            z_train = np.array([sph.mass(complex(x_train[i, 0], y_train[i, 0])) for i in range(data_size)])
            z_train = z_train.reshape(-1, 1)  # Converts (1000,) → (1000,1)

            self.output_tensor = torch.tensor(z_train, dtype=torch.float32)
            self.input_tensor = torch.cat([x_tensor, y_tensor], dim=1)
        else:
            raise ValueError('Somehow an invalid catagory was passed to generate_training_data_single_twist; this should have been caught earlier')
        
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




        self.model = FNN_model_3()
        self.catagory = catagory
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max



    def train(self, num_epochs):
        

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
        torch.save(self.model.state_dict(), filename)

    def load_model(self, filename):
        self.model.load_state_dict(torch.load(filename))

    def predictions_to_plotly(self, color_scale='viridis'):

         # Define x values (spread a  a region)
        x_vals = np.linspace(self.x_min, self.x_max, 150)  # X values from -5 to 5
        _NUM_Y_VALS_ = 100

        # Generate y values satisfying y > x^2
        y_vals = []
        for x in x_vals:
            y_range = np.linspace(self.y_min, self.y_max, _NUM_Y_VALS_)  # 50 points per x value
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

    def __init__(self, filename, catagory):

        class DataClass(Dataset):

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
            
        self.catagory = catagory
        self.model = FNN_model_3(catagory)

    def train(self, num_epochs=50):
        
        #Define loss function and optimizer
        criterion = nn.MSELoss()  # Mean Squared Error loss
        optimizer = optim.Adam(model.parameters(), lr=0.01)

        # Use learning rate scheduler
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=20, verbose=True)
        print_rounds = int(num_epochs / 10)
        outer_progress = tqdm(range(num_epochs), desc="Training Progress", unit="epoch")

        for epoch in outer_progress:
            inner_progress = tqdm(self.dataloader, desc=f"Epoch {epoch+1}/{epochs}", unit="batch", leave=False)

            for batch in inner_progress:
                input_tensor, output_tensor = batch
                optimizer.zero_grad()
                
                # Forward pass
                pred = model(input_tensor)
                
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
    """
    WARNING: Every iteration of append_training_data_to_CSV_single_twist with the default 20,000 samples
    makes a file that is roughly 1.5MB and takes about 1.83s (on 8 Intel I9 cores). Thus, ranging over
    all spherical twists Tw_a O(b)  with -5≤ a ≤ 5 and -5 ≤ b ≤ 5 will produce a file that is roughly 150MB and
    will take roughly 3 and a half minutes (on 8 Intel I9 cores).
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








        




class FNN_model_1(nn.Module):
    def __init__(self, catagory):

        if catagory == 'K3':
            self._input_size = 5
        elif catagory == 'P1':
            self._input_size = 4
        elif catagory == 'P2':
            self._input_size = 4
        else:
            raise ValueError('Invalid catagory')

        super(FNN_model_1, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(self._input_size, 64),  # Input: (x, y)
            nn.ReLU(),
            nn.Linear(64, 64),
            nn.ReLU(),
            nn.Linear(64, 1)  # Output: z
        )

    def forward(self, xy):
        return self.model(xy)
    
class FNN_model_2(nn.Module):
    def __init__(self, catagory):

        if catagory == 'K3':
            self._input_size = 5
        elif catagory == 'P1':
            self._input_size = 4
        elif catagory == 'P2':
            self._input_size = 4
        else:
            raise ValueError('Invalid catagory')
        
        super(FNN_model_2, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(self._input_size, 128),  # Input: (x, y)
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
    



###############################################
#         STATIC HELPER FUNCTIONS             #
###############################################

def _generate_xy_tensors_P2(x_min, x_max, y_max, data_size):
    """
    Generate x and y tensors for the P2 category.
    """
    DLP = LePotier()

    largest_abs_x = max(abs(x_min), abs(x_max))
    if y_max < DLP.curve_estimate(largest_abs_x):
        raise ValueError(f'y_max is too small for the range {(x_min, x_max)}')

    x_train, y_train = [], []
    while len(x_train) < data_size:
        x = np.random.uniform(x_min, x_max)
        y = np.random.uniform(-1, y_max)

        if DLP.is_above_curve(x, y):
            x_train.append(x)
            y_train.append(y)

    return x_train, y_train


def _append_training_data_to_CSV_single_twist(filename, line_bundle_1, line_bundle_2, catagory, degree=1,
                                            x_min=-5, x_max=5, y_min=0, y_max=10,
                                            data_size=20000):
    
    if catagory not in IMPLEMENTED_CATAGORIES:
        raise ValueError(f'Invalid catagory. Choose from {IMPLEMENTED_CATAGORIES}')
    
    # Construct the xy-input data
    x_train, y_train = [], []
    if catagory == 'P2':
        x_train, y_train = _generate_xy_tensors_P2(x_min, x_max, y_min, y_max, data_size)
    else:
        x_train = np.random.uniform(x_min, x_max, (data_size, 1))
        y_train = np.random.uniform(y_min, y_max, (data_size, 1))


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



if __name__ == '__main__':

    st_model = SingleTwistModel(line_bundle_1=1, line_bundle_2=2, catagory='K3', data_size=25000)
    
    st_model.train(num_epochs=5000)
    fig = st_model.predictions_to_plotly()
    fig.show()





    # Generate training data
    # filename =  os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/training_data_K3.csv'))
    # create_K3_training_data_single_twist(filename)
    

    
    

   