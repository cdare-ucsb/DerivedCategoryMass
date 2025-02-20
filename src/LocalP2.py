import math
import cmath
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np




class StabilityCondition():

    def __init__(self, s, q):
        self.s = s
        self.q = q


    def plot_central_charge(self, objects, labels=None):
        """
        Plots the central charge on the complex plane and draws a line from the origin (0,0) to each complex number.

        Args:
            central_charges (list of complex): A list of complex numbers to plot.
            labels (list of str, optional): Labels corresponding to each complex number.
        """

        central_charges = [obj.central_charge(self.s, self.q) for obj in objects]

        plt.figure(figsize=(6, 6))
        
        for i, z in enumerate(central_charges):

            stab_color = 'blue' if objects[i].is_semistable(self.s, self.q) else 'red'
            line_color = 'cyan' if objects[i].is_semistable(self.s, self.q) else 'red'

            # Plot point
            plt.scatter(z.real, z.imag, color=stab_color, marker='o')
            
            # Draw line from (0,0) to the complex number
            plt.plot([0, z.real], [0, z.imag], linestyle='-', color=line_color, alpha=0.7)

            # Add labels if provided
            if labels and i < len(labels):
                plt.text(z.real, z.imag, f" {labels[i]}", fontsize=12, verticalalignment='bottom')

        # Draw x-axis and y-axis
        plt.axhline(0, color='black', linewidth=1)
        plt.axvline(0, color='black', linewidth=1)
        plt.grid(True, linestyle="--", alpha=0.6)

        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        plt.title("Central Charge on Complex Plane")

        plt.show()




        

if __name__ == "__main__":
    


    sph = SphericalTwist(LineBundle(1), LineBundle(2))

    # Define x values (spread around a region)
    x_vals = np.linspace(-5, 5, 100)  # X values from -2 to 2

    # Generate y values satisfying y > x^2
    y_vals = []
    for x in x_vals:
        y_min = x**2 -0.4  # Slightly above x^2
        y_max = 25  # Arbitrary upper bound
        y_range = np.linspace(y_min, y_max, 50)  # 50 points per x value
        y_vals.append(y_range)

    # Convert to numpy array
    y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # Repeat x values to match the shape of y
    x_vals = np.repeat(x_vals, 50)  # Each x value repeats 10 times

    masses = np.array([sph.mass(x, y) for x, y in zip(x_vals, y_vals)])

    # Plot the surface
    fig = go.Figure(data=[go.Scatter3d(z=masses, x=x_vals, y=y_vals,
                                    mode='markers', marker=dict(size=4, color=masses, colorscale='viridis'))])

    fig.update_layout(
        title="Interactive 3D Mesh Plot",
        autosize=True,
        margin=dict(l=0, r=0, b=0, t=30)
    )

    fig.show()
    




