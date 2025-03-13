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
from numbers import Number

from SphericalTwist import SphericalTwist
from CoherentSheaf import LineBundle






class K3GeometricChamber():
    """!
    This class represents the geometric chamber of a general projective K3 surface of picard rank 1. Since the geometric chamber consists of a Cantor-set of walls associated to the spherical vector bundles on it, there is no accurate way to simultaneously represent all walls at once. Thus, this class provides a method to plot such walls in the (α,β) plane up to some accuracy (i.e. granularity).

    """

    def __init__(self, degree=1, granularity=15):
        r"""!
        Constructor for the K3GeometricChamber class. This class represents the geometric chamber of a general projective K3 surface of picard rank 1. Since the geometric chamber consists of a Cantor-set of walls associated to the spherical vector bundles on it, there is no accurate way to simultaneously represent all walls at once. Thus, this class provides a method to plot such walls in the (α,β) plane up to some accuracy (i.e. granularity).

        \param int degree The degree of the K3 surface

        \param int granularity The granularity of the walls to be plotted in the the (α,β) plane

        \throws ValueError If the degree is not a positive integer
        """

        if not isinstance(degree, int):
            raise ValueError("Degree must be an integer")
        if degree < 1:
            raise ValueError("Degree must be a positive integer")

        self.degree = degree ## The degree of the K3 surface
        
        self.granularity = granularity ## The granularity of the walls to be plotted in the the (α,β) plane


    def create_mukai_vectors(self):
        r"""!
        Method to create the mukai vectors of a spherical vector bundle on a K3 surface of degree d.

        \param int granularity The granularity of the mukai vectors to be computed. This roughly indicates how close the mukai vectors should get to integer points on the real axis

        \return list A list of tuples representing the mukai vectors of the spherical vector bundles
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
        """!
        Method to plot the walls of the geometric chamber of a K3 surface in the (α,β) plane. As this method is quite computationally intensive at the standard granularity of 15, we only plot the walls between -1.6 and 1.6 in the α direction and between 0 and √d in the β direction.

        \param bool return_json A flag indicating whether the plot should be returned as a JSON string or displayed in the browser
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
        r"""!
        Helper method to compute the divisors of a number N.

        \param int N The number whose divisors are to be computed

        \return numpy.ndarray An array of divisors of N
        """

        divs = np.arange(1,int(N**0.5)+1)  # potential divisors up to √N
        divs = divs[N % divs==0]           # divisors
        comp = N//divs[::-1]               # complement quotients
        return np.concatenate((divs,comp[divs[-1]==comp[0]:])) # combined


    def _convert_to_alpha_beta(self, r, l, s):
        r"""!
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


        \param int r The first component of the mukai vector
        \param int l The second component of the mukai vector
        \param int s The third component of the mukai vector

        \return tuple The coordinates of the hole in the (α,β) plane
        """


        alpha = float(l) / r
        beta = float(1 / r*math.sqrt(self.degree))

        return alpha, beta



def _calcZ1(x, y, k, n):
    r"""!
    Computes a single component of the complex hypersurface in the (x,y) plane.

    \param float x The x-coordinate of the point
    \param float y The y-coordinate of the point
    \param int k The index of the component
    \param int n The degree of the hypersurface

    \return complex The k-th component of the hypersurface at the point (x,y)
    """
    return cmath.exp(1j*(2*cmath.pi*k/n)) * (cmath.cos(x+y*1j)**(2/n))

def _calcZ2(x, y, k, n):
    r"""!
    Computes a single component of the complex hypersurface in the (x,y) plane.

    \param float x The x-coordinate of the point
    \param float y The y-coordinate of the point
    \param int k The index of the component
    \param int n The degree of the hypersurface

    \return complex The k-th component of the hypersurface at the point (x,y)
    """

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
    
    r"""!
    Example of a method to plot a complex degree d hypersurface using matplotlib, and save the animation either
    as an HTML <video> tag or as a GIF file. The hypersurface is given by the equation

    z0^d + z1^d + ... + zn^d = 1

    and is restricted to the region 0 <= z0, z1 <= pi/2. The plot is done using Matplotlib.

    \param int degree The degree of the hypersurface
    \param str filename The name of the file to save the plot to
    \param bool to_gif A flag indicating whether the plot should be saved as a GIF file
    \param bool to_html A flag indicating whether the plot should be saved as an HTML file
    \param int y_granularity The granularity of the y-axis
    \param int x_granularity The granularity of the x-axis
    \param int nframes The number of frames in the animation
    \param int t_interval: The time interval between frames

    \throws ValueError If both to_gif and to_html are both simultaneously set to True
    """
    
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
    
    """!
    First static example for plotting a degree d complex hypersurface; specifically, this just restricts the graph of
    
    z0^d + z1^d + ... + zn^d = 1

    to the region 0 <= z0, z1 <= pi/2. The plot is done using Plotly.

    \param int degree The degree of the hypersurface
    \param str filename The name of the file to save the plot to
    \param bool to_html A flag indicating whether the plot should be saved to an HTML file
    \param bool return_json A flag indicating whether the plot should be returned as a JSON string
    \param int y_granularity The granularity of the y-axis
    \param int x_granularity The granularity of the x-axis
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


    


        

    

def find_discontinuities_disc_Lapl_in_single_twist_mass_K3(line_bundle_1, line_bundle_2, degree=1,
                                                plot=False, colorscale='plasma',
                                                x_min=-5, x_max=5, x_granularity=0.1,
                                                y_min=0, y_max=5, y_granularity=0.05):
    

    r"""!
    Method to find the discontinuities in the mass of a single spherical twist on a K3 surface using the discretized
    Laplacian. The method computes the discretized Laplacian values at each point in the (x,y) plane, and then
    computes the average and standard deviation of the discretized Laplacian values. The method then highlights
    the points in the (x,y) plane that have a discretized Laplacian value greater than the average plus the standard
    deviation.

    \param int line_bundle_1 The line bundle of the object being twisted
    \param int line_bundle_2 The line bundle of the object being twisted by
    \param int degree The degree of the K3 surface
    \param plot A flag indicating whether the discontinuities should be plotted
    \param str colorscale The colorscale to be used in the plot
    \param float x_min The minimum x value for the plot
    \param float x_max The maximum x value for the plot
    \param float x_granularity The granularity of the x-axis
    \param float y_min The minimum y value for the plot
    \param float y_max The maximum y value for the plot
    \param float y_granularity The granularity of the y-axis

    \throws ValueError If the input data is not an integer

    \return list A list of tuples representing the discontinuities in the mass of the spherical twist
    """
    
    ##############################
    #   Input data validation    #
    ##############################

    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Line bundles must be integers")

    if not isinstance(degree, int):
        raise ValueError("Degree must be an integer")
    
    if not isinstance(plot, bool):
        raise ValueError("Plot must be a boolean")
    
    if not isinstance(x_granularity, Number) or not isinstance(y_granularity, Number):
        raise ValueError("Granularity must be a number")
    
    if y_granularity <= 0 or x_granularity <= 0:
        raise ValueError("Granularity must be a positive number")

    
    if y_min - y_granularity <=0:
        # We will ultimately need to compute the neighbors for the discrete Laplacian
        y_min = y_granularity

    
    ##############################
    #        Build Grid          #
    ##############################     
    
    # Define x values (spread evenly over numpy aregion)
    x_vals = np.arange(x_min, x_max, x_granularity)  
    y_vals = []
    
    for x in x_vals:
        y_range = np.arange(y_min, y_max, y_granularity)  
        y_vals.append(y_range)

    y_vals = np.array(y_vals).flatten()  # Flatten the y array to be 1 dimensional
    x_vals = np.repeat(x_vals, len(np.arange(y_min, y_max, y_granularity)))  # Each x value repeats y_vals_len number of times to make a grid.



    ##############################
    #    Compute Discretized     #
    #         Laplacian          #
    ##############################

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                          LineBundle(line_bundle_2, catagory='K3'),
                          degree=degree)

    def _disc_Lapl(x, y, degree):
        return sph.mass(x + x_granularity, y, degree) + sph.mass(x - x_granularity, y, degree) + sph.mass(x, y + y_granularity, degree) + sph.mass(x, y - y_granularity, degree) - 4*sph.mass(x, y, degree)
    

    disc_vals = np.array([_disc_Lapl(x, y, degree) for x, y in zip(x_vals, y_vals)])
    disc_vals_mean = np.mean(disc_vals)
    disc_vals_std = np.std(disc_vals)


    ##############################
    #    Find Discontinuities    #
    ##############################

    discontinuities = []
    for x, y, disc_val in zip(x_vals, y_vals, disc_vals):
        # Our threshold will be 1 standard deviation above the mean
        if abs(disc_val - disc_vals_mean) > disc_vals_std:
            discontinuities.append((x, y))




    if plot:

        ##############################
        #         Plot Data          #
        ##############################

        colors = np.abs(disc_vals - disc_vals_mean)  # Highlight stronger discontinuities
        fig = go.Figure(data=[go.Scatter3d(z=disc_vals, x=x_vals, y=y_vals,
                                        mode='markers', marker=dict(size=3, color=colors, colorscale=colorscale))])

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
    
    return discontinuities


        
        

    

    
    


    






if __name__ == "__main__":

    # K3 = K3GeometricChamber(degree=1)
    
    # complex_hypersurface_matplot_animation_ex1(filename='CY3fold.gif', degree=5, to_gif=True)

    # complex_hypersurface_plotly_ex1(degree=4)

    y_granularity=0.05
    x_granularity=0.1

    # discontinuities, colors = find_discontinuities_std_dev_in_single_twist_mass_K3(1, 2, degree=1)


    sph = SphericalTwist(LineBundle(1, catagory='K3'),
                          LineBundle(2, catagory='K3'),
                          degree=1)
    
    # Define x values (spread a  a region)
    x_vals = np.arange(-5, 5, x_granularity)  # X values from -5 to 5

    y_vals = []
    

    for x in x_vals:
        y_range = np.arange(0.1, 5, y_granularity)  # 50 points per x value
        y_vals_len = len(y_range)

        y_vals.append(y_range)

    y_vals = np.array(y_vals).flatten()  # Flatten the y array
    x_vals = np.repeat(x_vals, y_vals_len)  # Each x value repeats 10 times

    
    masses = np.array([sph.mass(x, y, 1) for x, y in zip(x_vals, y_vals)])

    # Plot the surface
    fig = go.Figure(data=[go.Scatter3d(z=masses, x=x_vals, y=y_vals,
                                    mode='markers', marker=dict(size=3, color=masses, colorscale='viridis' ))])

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


    discont = find_discontinuities_disc_Lapl_in_single_twist_mass_K3(1, 2, degree=1,
                                                                    plot=True, colorscale='Viridis')

    print(discont)
