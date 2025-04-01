import numpy as np
from numbers import Number
import plotly.graph_objects as go
import json
import plotly


from ..SphericalTwist.SphericalTwist import SphericalTwist, DoubleSphericalTwist
from ..CoherentSheaf.CoherentSheaf import LineBundle
from ..LocalP2 import LePotier



from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']




class MassPlot():
    r"""!
    As the purpose of the library is to compute the mass of spherical twists in certain chambers of 
    stability manifolds, a useful functionality is to plot a certain 3D slice of the graph of 
    the mass function. Since the SphericalTwist class allows significant flexibiliy with its 
    catagory parameter (which effectively dictates the stability manifold we wish to graph), 
    the majority of the plotting functionality will look the same for each catagory (currently
    'P1', 'P2', and 'K3'). The only difference will be the way the x and y values are generated
    for the plot. Thus, this class acts as wrapper to remove that redundancy and allow for easy
    plotting of the mass function for a given spherical twist in a given chamber of a stability manifold.
    """
    

    def __init__(self, line_bundle_1, line_bundle_2, catagory, # Required arguments
                line_bundle_3 = None, degree=1,                # Optional for double twist on K3s
                x_min=-5, x_max=5, x_granularity=0.1,          # Optional for x-axis range and granularity
                y_min=0, y_max=5, y_granularity=0.05):         # Optional for y-axis range and granularity
        r"""!
        Initialize the MassPlot object with the given parameters. The required parameters are the line bundles
        for the spherical twist, the catagory of the stability manifold, and the x and y ranges and granularities
        for the plot. The optional parameters are the degree of the K3 surface (default is 1), and a third line
        bundle for the double spherical twist (default is None).

        \param int line_bundle_1 The degree of the outermost line bundle in the spherical twist
        \param int line_bundle_2 The degree of the line bundle being twisted (if line_bundle_3 is None) or the degree of the middle line bundle in the double twist
        \param str catagory The catagory of the stability manifold
        \param int line_bundle_3 The degree of the line bundle being twisted in the double twist (default is None; setting this will cause the object to compute a double twist)
        \param int degree The degree of the K3 surface (default is 1)
        \param Number x_min The minimum x value for the plot (default is -5)
        \param Number x_max The maximum x value for the plot (default is 5)
        \param Number x_granularity The distance between points horizontally for the plot (default is 0.1)
        \param Number y_min The minimum y value for the plot (default is 0)
        \param Number y_max The maximum y value for the plot (default is 5)
        \param Number y_granularity The distance between points vertically for the plot (default is 0.05)
        """
        

        if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
            raise ValueError("Line bundles must be integers")
        
        self.line_bundle_1 = line_bundle_1 ## The degree of the outermost line bundle in the spherical twist

        self.line_bundle_2 = line_bundle_2 ## Either the degree of the line bundle being twisted (if line_bundle_3 is None) or the degree of the middle line bundle in the double twist

        if line_bundle_3 is not None and not isinstance(line_bundle_3, int):
            raise ValueError("Line bundle 3 must be an integer")
        if line_bundle_3 is not None and catagory not in __CURRENT_DOUBLE_TWIST_IMPLEMENTED__:
            raise ValueError("Double spherical twist not implemented for this catagory")
        elif line_bundle_3 is not None:
            self._num_twists = 2 ## Number of twists in the spherical twist
        else:
            self._num_twists = 1

        self.line_bundle_3 = line_bundle_3 ## Optional argument that allows the flexibility of double spherical twists

        if not isinstance(degree, int):
            raise ValueError("Degree must be an integer")
        if degree < 1:
            raise ValueError("Degree must be a positive integer")
        
        self.degree = degree ## Optional argument used for the degree of a K3 surface; default is 1. Not used for P1 and P2

        
        if not isinstance(x_granularity, Number) or not isinstance(y_granularity, Number):
            raise ValueError("Granularity must be a number")
        
        if y_granularity <= 0 or x_granularity <= 0:
            raise ValueError("Granularity must be a positive number")
        

        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError("Catagory not implemented")
        
        self.catagory = catagory  ## Catagory of the stability manifold

        if y_min - y_granularity <=0:
            # We will ultimately need to compute the neighbors for the discrete Laplacian
            y_min = y_granularity

        if y_min >= y_max:
            raise ValueError("y_min must be less than y_max")
        if x_min >= x_max:
            raise ValueError("x_min must be less than x_max")

        
        self.x_min = x_min  ## The minimum x value for the plot

        self.x_max = x_max ## The maximum x value for the plot

        self.x_granularity = x_granularity  ## The distance between points horizontally for the plot

        self.y_min = y_min ## The minimum y value for the plot

        self.y_max = y_max ## The maximum y value for the plot

        self.y_granularity = y_granularity ## The distance between points vertically for the plot

        self._xvals = np.arange(x_min, x_max, x_granularity) ## The list of x values for the plot

        self._yvals = [] ## The list of y values for the plot

        if self.catagory != 'P2':
            # Catagory is P1 or K3, use a rectangular grid
            for _ in self._xvals:
                self._yvals.append(  np.arange(y_min, y_max, y_granularity)  )

            self._yvals = np.array(self._yvals).flatten()  # Flatten the y array to be 1 dimensional
            self._xvals = np.repeat(self._xvals, len(np.arange(y_min, y_max, y_granularity)))  # Each x value repeats y_vals_len number of times to make a grid.
        else:
            # Catagory is P2, have to use DLP curve
            # We will no longer use y_min and y_max here since we know the shape of the DLP curve

            DLP = LePotier()
            
            max_abs_x = max(abs(x_min), abs(x_max))
            y_max = max(DLP.curve_estimate(max_abs_x), y_max)  # y_max must be greater than the Drezet-LePotier curve over the x-values chosen

            
            _lin_sp_const = int( (y_max + 0.5) / y_granularity )  # Linear space constant
            
            for x in self._xvals:
                y = float(DLP.curve_estimate(x)) # Get corresponding boundary point on Drezet-LePotier Curve
                if y_max - y < 0:
                    raise ValueError("y_max must be greater than the Drezet-LePotier curve over the x-values chosen")

            
                self._yvals.append(np.linspace(y, y_max, _lin_sp_const))  # Linear space from y to y_max

            # Convert to numpy array
            self._yvals = np.array(self._yvals).flatten()  # Flatten the y array

            # Repeat x values to match the shape of y
            self._xvals = np.repeat(self._xvals, _lin_sp_const)  # Each x value repeats _lin_sp_const times




    def get_mass_values(self):
        r"""!
            Compute the mass values for each point in the grid formed by the x and y values.
            This is done by creating an instance of the SphericalTwist class and then calling
            the mass method for each point in the grid.

            \return np.array The mass values for each point in the grid
        """

        sph = None
        
        if self._num_twists == 1:
            sph = SphericalTwist(LineBundle(self.line_bundle_1, catagory=self.catagory),
                          LineBundle(self.line_bundle_2, catagory=self.catagory),
                          degree=self.degree)

        elif self._num_twists == 2:
            sph = DoubleSphericalTwist(LineBundle(self.line_bundle_1, catagory=self.catagory),
                        LineBundle(self.line_bundle_2, catagory=self.catagory),
                        LineBundle(self.line_bundle_3, catagory=self.catagory),
                        degree=self.degree)
        else:
            raise ValueError("Invalid number of twists")
        
        if self.catagory == 'P2':
            return np.array([sph.mass(x, y) for x, y in zip(self._xvals, self._yvals)])
        elif self.catagory == 'P1':
            return np.array([sph.mass(complex(x, y)) for x, y in zip(self._xvals, self._yvals)])
        elif self.catagory == 'K3':
            return np.array([sph.mass(x, y, self.degree) for x, y in zip(self._xvals, self._yvals)])
        else:
            raise ValueError("Invalid catagory")
    

    def to_plotly_figure(self, color_scale='viridis', marker_size=3):
        r"""!
        The primariy functionality of this class; generate a plotly figure of the mass function
        for the given spherical twist in the given chamber of the stability manifold. The mass values
        are computed for each point in the grid formed by the x and y values, and then plotted as a 3D
        scatter plot.

        \param str color_scale The color scale to use for the plot
        \param int marker_size The size of the markers in the plot
        
        \return go.Figure The plotly figure of the mass function
        """

        masses = self.get_mass_values()

        # Plot the surface
        fig = go.Figure(data=[go.Scatter3d(z=masses, x=self._xvals, y=self._yvals,
                                        mode='markers', marker=dict(size=marker_size, color=masses, colorscale=color_scale))])

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
    


    def to_plotly_json(self, color_scale='viridis', marker_size=3):
        r"""!
        Function to convert the plotly figure to a JSON string. This is useful for rendering the plot
        in a web application.

        \param str color_scale The color scale to use for the plot
        \param int marker_size The size of the markers in the plot

        \return str The JSON string of the plotly figure
        """

        fig = self.to_plotly_figure(color_scale=color_scale, marker_size=marker_size)
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    


    
    def get_discretized_Laplacian_values(self):
        r"""!
        Compute the values of the discretized Laplacian for each point in the grid formed by the x and y values.
        This is done by creating an instance of the SphericalTwist class and then calling the mass method
        for each point in the grid.

        \return np.array The values of the discretized Laplacian for each point in the grid
        """


        sph = None
        
        if self._num_twists == 1:
            sph = SphericalTwist(LineBundle(self.line_bundle_1, catagory=self.catagory),
                          LineBundle(self.line_bundle_2, catagory=self.catagory),
                          degree=self.degree)

        elif self._num_twists == 2:
            sph = DoubleSphericalTwist(LineBundle(self.line_bundle_1, catagory=self.catagory),
                        LineBundle(self.line_bundle_2, catagory=self.catagory),
                        LineBundle(self.line_bundle_3, catagory=self.catagory),
                        degree=self.degree)
        else:
            raise ValueError("Invalid number of twists")
        

        def _disc_Lapl(x, y):
            if self.catagory == 'P2':
                return sph.mass(x + self.x_granularity, y) + sph.mass(x - self.x_granularity, y) + sph.mass(x, y + self.y_granularity) + sph.mass(x, y - self.y_granularity) - 4*sph.mass(x, y)
            elif self.catagory == 'P1':
                return sph.mass(complex(x + self.x_granularity, y)) + sph.mass(complex(x - self.x_granularity, y)) + sph.mass(complex(x, y + self.y_granularity)) + sph.mass(complex(x, y - self.y_granularity)) - 4*sph.mass(complex(x, y))
            elif self.catagory == 'K3':
                return sph.mass(x + self.x_granularity, y, self.degree) + sph.mass(x - self.x_granularity, y, self.degree) + sph.mass(x, y + self.y_granularity, self.degree) + sph.mass(x, y - self.y_granularity, self.degree) - 4*sph.mass(x, y, self.degree)
            else:
                raise ValueError("Invalid catagory")
            

        return np.array([_disc_Lapl(x, y) for x, y in zip(self._xvals, self._yvals)])
    

    def discretized_Laplacian_to_plotly_figure(self, color_scale='Plasma', marker_size=3):
        r"""!
        Generate a plotly figure of the discretized Laplacian for the given spherical twist in the given
        chamber of the stability manifold. The values of the discretized Laplacian are computed for each point
        in the grid formed by the x and y values, and then plotted as a 3D scatter plot.

        \param str color_scale The color scale to use for the plot
        \param int marker_size The size of the markers in the plot

        \return go.Figure The plotly figure of the discretized Laplacian
        """

        disc_vals = self.get_discretized_Laplacian_values()
        disc_vals_mean = np.mean(disc_vals)

        colors = np.abs(disc_vals - disc_vals_mean)  # Highlight stronger discontinuities
        fig = go.Figure(data=[go.Scatter3d(z=disc_vals, x=self._xvals, y=self._yvals,
                                        mode='markers', marker=dict(size=marker_size, color=colors, colorscale=color_scale))])

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




    
    def find_discontinuities(self, std_dev_threshold=1):

        disc_vals = self.get_discretized_Laplacian_values()
        disc_vals_mean = np.mean(disc_vals)
        disc_vals_std = np.std(disc_vals)


        ##############################
        #    Find Discontinuities    #
        ##############################

        discontinuities = []
        for x, y, disc_val in zip(self._xvals, self._yvals, disc_vals):
            # Our threshold will be 1 standard deviation above the mean
            if abs(disc_val - disc_vals_mean) > std_dev_threshold*disc_vals_std:
                discontinuities.append((x, y))

        return discontinuities
    



########################################
#               Main                   #
########################################


if __name__ == "__main__":
    
    mp = MassPlot(line_bundle_1=1, line_bundle_2=5,
                x_min=-10, x_max=10,
                y_min=0, y_max = 20,
                catagory='K3', degree=1)

    fig1 = mp.to_plotly_figure()
    fig1.show()

    fig2 = mp.discretized_Laplacian_to_plotly_figure()
    fig2.show()

