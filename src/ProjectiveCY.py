import numpy as np
import cmath
import math
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
import cmath
import plotly.graph_objects as go
from plotly.graph_objs import *
import plotly.utils
import json

from .SphericalTwist import SphericalTwist, DoubleSphericalTwist
from .CoherentSheaf import LineBundle






class K3GeometricChamber():
    """
    This class represents the geometric chamber of a general projective K3 surface of picard rank 1. Since the geometric chamber consists of a Cantor-set of walls associated to the spherical vector bundles on it, there is no accurate way to simultaneously represent all walls at once. Thus, this class provides a method to plot such walls in the (α,β) plane up to some accuracy (i.e. granularity).

    Attributes
    ----------
    degree : int
        The degree of the K3 surface
    granularity : int
        The granularity of the walls to be plotted in the the (α,β) plane
    """
    

    def __init__(self, degree=1, granularity=15):
        """
        Constructor for the K3GeometricChamber class. This class represents the geometric chamber of a general projective K3 surface of picard rank 1. Since the geometric chamber consists of a Cantor-set of walls associated to the spherical vector bundles on it, there is no accurate way to simultaneously represent all walls at once. Thus, this class provides a method to plot such walls in the (α,β) plane up to some accuracy (i.e. granularity).

        Parameters
        ----------
        degree : int
            The degree of the K3 surface
        granularity : int
            The granularity of the walls to be plotted in the the (α,β) plane

        Raises
        ------
        ValueError
            If the degree is not a positive integer
        
        """

        if not isinstance(degree, int):
            raise ValueError("Degree must be an integer")
        if degree < 1:
            raise ValueError("Degree must be a positive integer")

        self.degree = degree
        self.granularity = granularity


    def create_mukai_vectors(self):
        """
        Method to create the mukai vectors of a spherical vector bundle on a K3 surface of degree d.


        Parameters
        ----------
        granularity : int
            The granularity of the mukai vectors to be computed. This roughly indicates how close the mukai vectors should get to integer points
            on the real axis

        Returns
        -------
        mukai_vectors : list
            A list of tuples representing the mukai vectors of the spherical vector bundles
        """

        mukai_vectors = []

        for n in range(-self.granularity, self.granularity):

            numerator = self.degree*n**2 + 1
            divisors = self._npDivs(numerator)
            for div in divisors:
                r = div
                s = numerator // div

                mukai_vectors.append((r, n, s))

        return mukai_vectors

    def plot_alpha_beta_plane(self, return_json=False):
        """
        Method to plot the walls of the geometric chamber of a K3 surface in the (α,β) plane. As this method is quite computationally intensive at the standard granularity of 15, we only plot the walls between -1.6 and 1.6 in the α direction and between 0 and √d in the β direction.

        Parameters
        ----------
        return_json : bool
            A flag indicating whether the plot should be returned as a JSON string or displayed in the browser

        
        """

        x_vals = []
        y_vals = []

        mukai_vectors = self.create_mukai_vectors()

        for r, l, s in mukai_vectors:
            alpha, beta = self._convert_to_alpha_beta(r, l, s)
            x_vals.append(alpha)
            y_vals.append(beta)

            x_vals.append(alpha+1)
            y_vals.append(beta)

            x_vals.append(alpha-1)
            y_vals.append(beta)
        
        fig = go.Figure(data=[go.Scatter(x=x_vals, y=y_vals, mode='markers')])

        # Add vertical lines from x-axis to each point
        for x_val, y_val in zip(x_vals, y_vals):
            fig.add_shape(
                type="line",
                x0=x_val,
                y0=0,
                x1=x_val,
                y1=y_val,
                line=dict(
                    color="LightSeaGreen",
                    width=1
                )
            )

        fig.update_layout(
            title="Mukai Vectors",
            xaxis_title="α",
            yaxis_title="β",
            showlegend=False
        )
        fig.update_xaxes(range=[-1.6, 1.6])
        fig.update_yaxes(range=[-0.2, math.sqrt(self.degree)+0.2])

        if return_json:
            return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        else:
            fig.show()



    def _npDivs(self, N):
        """
        Helper method to compute the divisors of a number N.

        Parameters
        ----------
        N : int
            The number whose divisors are to be computed

        Returns
        -------
        divs : numpy.ndarray
            An array of divisors of N
        """

        divs = np.arange(1,int(N**0.5)+1)  # potential divisors up to √N
        divs = divs[N % divs==0]           # divisors
        comp = N//divs[::-1]               # complement quotients
        return np.concatenate((divs,comp[divs[-1]==comp[0]:])) # combined


    def _convert_to_alpha_beta(self, r, l, s):
        """
        Helper method to compute the coordinates of a hole in the (α,β) plane corresponding to a mukai vector.


        Suppose (r, lH, s) is the mukai vector of a spherical vector bundle on a degree d K3 surface (H^2 = 2d). Then
        # dl^2 - rs = -1. For a corresponding hole in the (α,β) plane, we must have the real and imaginary part of 
        # < exp(αΗ + ι βΗ) , (r, lH, s) > vanish. This gives

        0 = <(1, αH, (α^2 - β^2)d), (r, lH, s)> = 2dal - s - rd(α^2 - β^2)
        0 = <(0, βΗ, 2αβd), (r, lH, s)> = 2βdl - 2rdαβ = 2dβ(l - rα)
        
        The second equation gives β = 0 or l = rα. Since βΗ is necessarily ample, we must have β > 0.
        Therefore, l = rα. Substituting into the first equation simplifies to

        dl^2 - rs + dr^2 β^2 = 0

        Using dl^2 - rs = -1, we find that β = 1/r√d. 

        Parameters
        ----------
        r : int
            The first component of the mukai vector
        l : int
            The second component of the mukai vector
        s : int
            The third component of the mukai vector

        Returns
        -------
        alpha : float
            The x-coordinate of the hole in the (α,β) plane
        beta : float
            The y-coordinate of the hole in the (α,β) plane
        """


        alpha = float(l) / r
        beta = float(1 / r*math.sqrt(self.degree))

        return alpha, beta



def _calcZ1(x, y, k, n):
    return cmath.exp(1j*(2*cmath.pi*k/n)) * (cmath.cos(x+y*1j)**(2/n))

def _calcZ2(x, y, k, n):
    return cmath.exp(1j*(2*cmath.pi*k/n)) * (cmath.sin(x+y*1j)**(2/n))

def _calcZ1Real(x, y, k, n):
    return (_calcZ1(x, y, k, n)).real

def _calcZ2Real(x, y, k, n):
    return (_calcZ2(x, y, k, n)).real

def _calcZ(x, y, k1_, k2_, n, a_):
    z1 = _calcZ1(x, y, k1_, n)
    z2 = _calcZ2(x, y, k2_, n)
    return z1.imag * math.cos(a_) + z2.imag*math.sin(a_)

def _calcZ_alt(x,y,k1_, k2_, n, a_):
    z1 = _calcZ1(x, y, k1_, n)
    z2 = _calcZ2(x, y, k2_, n)
    return z1.imag * math.sin(a_) - z2.imag*math.cos(a_)



def complex_hypersurface_matplot_animation_ex1(degree, filename='hypersurf',
                                            to_gif=False, to_html=False,
                                            y_granularity=30, x_granularity=30,
                                            nframes=100, t_interval=175):
    
    if to_gif and to_html:
        raise ValueError(f"Cannot make {filename} both a .gif and .html file. One value must be set False.")
    

    # set param range
    x = np.linspace(0, math.pi/2, x_granularity)
    y = np.linspace(-math.pi/2, math.pi/2, y_granularity)
    x, y = np.meshgrid(x, y)

    fig = plt.figure(dpi=100, figsize=(10, 5))

    ax = fig.add_subplot(projection='3d')

    def update(t):
        ax.cla()

        for k1 in range(degree):
            for k2 in range(degree):

                    X = np.frompyfunc(_calcZ1Real, 4, 1)(x, y, k1, degree).astype('float32')
                    Y = np.frompyfunc(_calcZ2Real, 4, 1)(x, y, k2, degree).astype('float32')
                    Z1 = np.frompyfunc(_calcZ, 6, 1)(x, y, k1, k2, degree, t/10).astype('float32')
                    Z2 = np.frompyfunc(_calcZ_alt, 6, 1)(x, y, k1, k2, degree, t/10).astype('float32')

                    ax.plot_surface(X, Y, Z2, alpha=0.8, cmap=cm.afmhot)
        

        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_zlim(-2, 2)

    ani = animation.FuncAnimation(fig = fig, func = update, frames = nframes, interval = t_interval)
    
    if to_gif:
        writergif = animation.PillowWriter(fps=30) 
        ani.save(filename, writer=writergif)
    elif to_html:
        with open(filename, "w") as f:
            print(ani.to_html5_video(), file=f)
    else:
        plt.show()  







#####################################################
#            STATIC PLOTTING EXAMPLES               #
#####################################################



def complex_hypersurface_plotly_ex1(degree, filename="hypersurface.html", 
                                    to_html=False, return_json=False,
                                    y_granularity=30, x_granularity=30):
    
    """
    First static example for plotting a degree d complex hypersurface; specifically, this just restricts the graph of
    
    z0^d + z1^d + ... + zn^d = 1

    to the region 0 <= z0, z1 <= pi/2. The plot is done using Plotly.

    Args:
    ----------------
    degree (int): The degree of the hypersurface
    filename (str): The name of the file to save the plot to
    to_html (bool): A flag indicating whether the plot should be saved to an HTML file
    return_json (bool): A flag indicating whether the plot should be returned as a JSON string
    y_granularity (int): The granularity of the y-axis
    x_granularity (int): The granularity of the x-axis
    """
    

    x = np.linspace(0, math.pi/2, x_granularity)
    y = np.linspace(-math.pi/2, math.pi/2, y_granularity)
    x, y = np.meshgrid(x, y)

    fig = go.Figure()

    for k1 in range(degree):
        for k2 in range(degree):
            
            X = np.frompyfunc(_calcZ1Real, 4, 1)(x, y, k1, degree).astype('float32')
            Y = np.frompyfunc(_calcZ2Real, 4, 1)(x, y, k2, degree).astype('float32')
            Z1 = np.frompyfunc(_calcZ, 6, 1)(x, y, k1, k2, degree, 0).astype('float32')
            Z2 = np.frompyfunc(_calcZ_alt, 6, 1)(x, y, k1, k2, degree, 0).astype('float32')

            fig.add_trace(go.Surface(x=X, y=Y, z=Z1, showscale=False, colorscale='blues'))

    if to_html:
        fig.write_html(filename, full_html=False)
    if not return_json:
        fig.show()
    else:
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    



#####################################################
#                  FLASK METHODS                    #
#####################################################





def ints_to_mass_plot_K3_sing_twist(line_bundle_1, line_bundle_2, degree=1, return_json=False):
    """
    Helper method to plot the mass of a single spherical twist on a K3 surface, and potentially pass the output as a JSON string. This is primarily used for the Flask application, which feeds it to an HTML template.

    Args:
    ----------------
    line_bundle_1 (int): The line bundle of the object being twisted
    line_bundle_2 (int): The line bundle of the object being twisted by
    degree (int): The degree of the K3 surface
    return_json (bool): A flag indicating whether the plot should be returned as a JSON string

    Returns:
    ----------------
    str: A JSON string representation of the plot

    Raises:
    ----------------
    ValueError: If the input data is not an integer
    """

    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Input data must be integers")

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                          LineBundle(line_bundle_2, catagory='K3'),
                          degree=degree)


    

    # Define x values (spread around a region)
    x_vals = np.linspace(-5, 5, 200)  # X values from -2 to 2

    # Generate y values satisfying y > x^2
    y_vals = []
    for x in x_vals:
        y_range = np.linspace(0.1, 10, 160)  # 50 points per x value
        y_vals.append(y_range)

    # Convert to numpy array
    y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # Repeat x values to match the shape of y
    x_vals = np.repeat(x_vals, 160)  # Each x value repeats 10 times

    masses = np.array([sph.mass(x, y, degree) for x, y in zip(x_vals, y_vals)])

    # Plot the surface
    fig = go.Figure(data=[go.Scatter3d(z=masses, x=x_vals, y=y_vals,
                                    mode='markers', marker=dict(size=3, color=masses, colorscale='viridis'))])

    fig.update_layout(
        title="",
        autosize=True,
        margin=dict(l=0, r=0, b=0, t=30),
        scene=dict(
            
            bgcolor="white",  # Changes the 3D plot background,

            xaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white",),
            yaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white"),
            zaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white"
            )
        )
    )

    if return_json:
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    else:
        fig.show()



def ints_to_mass_plot_K3_double_twist(line_bundle_1, line_bundle_2, line_bundle_3, degree=1, return_json=False):
    """
    Helper method to plot the mass of a double spherical twist on a K3 surface, and potentially pass the output as a JSON string. This is primarily used for the Flask application, which feeds it to an HTML template.

    Args:
    ----------------
    line_bundle_1 (int): The line bundle of the object being twisted
    line_bundle_2 (int): The line bundle of the object being twisted by
    line_bundle_3 (int): The line bundle of the object being twisted by
    degree (int): The degree of the K3 surface
    return_json (bool): A flag indicating whether the plot should be returned as a JSON string

    Returns:
    ----------------
    str: A JSON string representation of the plot

    Raises:
    ----------------
    ValueError: If the input data is not an integer
    """

    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int) or not isinstance(line_bundle_3, int):
        raise ValueError("Input data must be integers")

    sph = DoubleSphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                          LineBundle(line_bundle_2, catagory='K3'),
                          LineBundle(line_bundle_3, catagory='K3'),
                          degree=degree)


    

    # Define x values (spread around a region)
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

    masses = np.array([sph.mass(x, y, degree) for x, y in zip(x_vals, y_vals)])

    # Plot the surface
    fig = go.Figure(data=[go.Scatter3d(z=masses, x=x_vals, y=y_vals,
                                    mode='markers', marker=dict(size=3, color=masses, colorscale='viridis'))])

    fig.update_layout(
        title="",
        autosize=True,
        margin=dict(l=0, r=0, b=0, t=30),
        scene=dict(
            
            bgcolor="white",  # Changes the 3D plot background,

            xaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white",),
            yaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white"),
            zaxis = dict(
                backgroundcolor="white",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="white"
            )
        )
    )

    if return_json:
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    else:
        fig.show()




def single_twist_triangle_to_json_K3(line_bundle_1, line_bundle_2, degree=1):
    """
    Helper function to convert the data of a spherical twist triangle to a JSON string. The data includes the
     sheaf vectors, shift vectors, and dimension vectors of the objects in the triangle. This is primarily used
     for the front-end visualization of the spherical twist triangle in a Flask application.

    Args:
    ----------------
    line_bundle_1 (int): The line bundle of the first object in the spherical twist triangle
    line_bundle_2 (int): The line bundle of the second object in the spherical twist triangle

    Returns:
    ----------------
    str: A JSON string representation of the spherical twist triangle data
    """
    
    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Input data must be integers")

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                          LineBundle(line_bundle_2, catagory='K3'),
                          degree=degree)
    first_sheaf_vector = []

    if len(sph.defining_triangle.object1.sheaf_vector) == 1:
        first_sheaf_vector = [line_bundle_1]
    else:
        first_sheaf_vector = [line_bundle_1, line_bundle_1]


    object1 = {
        "sheaf_vector" : first_sheaf_vector,
        "shift_vector" : sph.defining_triangle.object1.shift_vector,
        "dimension_vector" : sph.defining_triangle.object1.dimension_vector
    }

    object2 = {
        "sheaf_vector" : [line_bundle_2],
        "shift_vector" : sph.defining_triangle.object2.shift_vector,
        "dimension_vector" : sph.defining_triangle.object2.dimension_vector
    }

    chain_complex_data = {
        "object1" : object1,
        "object2" : object2
    }

    
    return json.dumps(chain_complex_data)        


def double_twist_triangle_to_json_K3(line_bundle_1, line_bundle_2, line_bundle_3, degree=1):
    """
    Helper function to convert the data of a spherical twist triangle to a JSON string. The data includes the
     sheaf vectors, shift vectors, and dimension vectors of the objects in the triangle. This is primarily used
     for the front-end visualization of the spherical twist triangle in a Flask application.

    Args:
    ----------------
    line_bundle_1 (int): The line bundle of the first object in the spherical twist triangle
    line_bundle_2 (int): The line bundle of the second object in the spherical twist triangle

    Returns:
    ----------------
    str: A JSON string representation of the spherical twist triangle data
    """
    
    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int) or not isinstance(line_bundle_3, int):
        raise ValueError("Input data must be integers")

    sph = DoubleSphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                          LineBundle(line_bundle_2, catagory='K3'),
                          LineBundle(line_bundle_3, catagory='K3'),
                          degree=degree)
    first_sheaf_vector = []

    if len(sph.defining_triangle.object1) == 1:
        first_sheaf_vector = [line_bundle_1]
    elif len(sph.defining_triangle.object1) == 2:
        first_sheaf_vector = [line_bundle_1, line_bundle_1]
    else:
        first_sheaf_vector = [line_bundle_1, line_bundle_1, line_bundle_1]


    object1 = {
        "sheaf_vector" : first_sheaf_vector,
        "shift_vector" : sph.defining_triangle.object1.shift_vector,
        "dimension_vector" : sph.defining_triangle.object1.dimension_vector
    }


    chain_complex_data = {
        "object1" : object1,
        "degrees" : [line_bundle_1, line_bundle_2, line_bundle_3]
    }

    
    return json.dumps(chain_complex_data)        






if __name__ == "__main__":

    K3 = K3GeometricChamber(degree=1)
    
    complex_hypersurface_matplot_animation_ex1(filename='CY3fold.gif', degree=5, to_gif=True)

    # complex_hypersurface_plotly_ex1(degree=4)
