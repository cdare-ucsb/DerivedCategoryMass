from CoherentSheaf import LineBundle
from SphericalTwist import SphericalTwist
import numpy as np
import json
import plotly.graph_objects as go
import plotly.utils


def ints_to_mass_plot_P1_sing_twist(line_bundle_1, line_bundle_2, return_json = False):
    """
    Method which takes in two integers and returns a plot of the mass of the spherical twist
    of the two line bundles. This is used in the Flask app implementation to generate plots
    based on user form data.

    Parameters
    ----------
    line_bundle_1 : int
        The first line bundle degree in the spherical twist
    line_bundle_2 : int
        The second line bundle degree in the spherical twist
    return_json : bool
        A flag to indicate whether the plot should be returned as a JSON string or displayed
        in the browser. When passed to a Flask app, this should be True. Default is False.

    Returns
    -------
    str
        A JSON string representation of the plotly plot if return_json is True
    """

    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Input data must be integers")

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='P1'),
                          LineBundle(line_bundle_2, catagory='P1'), degree=1)


    # Define x values (spread around a region)
    x_vals = np.linspace(-5, 5, 200)  # X values from -2 to 2

    # Generate y values satisfying y > x^2
    y_vals = []
    for x in x_vals:
        y_range = np.linspace(0.1, 5, 160)  # 50 points per x value
        y_vals.append(y_range)

    # Convert to numpy array
    y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # Repeat x values to match the shape of y
    x_vals = np.repeat(x_vals, 160)  # Each x value repeats 10 times

    masses = np.array([sph.mass(complex(x, y)) for x, y in zip(x_vals, y_vals)])

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


def twist_triangle_to_json_P1(line_bundle_1, line_bundle_2):
    """
    Helper function to convert the data of a spherical twist of two line bundles to a JSON string
    for use in the Flask app. This is used to display the chain complex data of the spherical twist
    in the browser.

    Parameters
    ----------
    line_bundle_1 : int
        The first line bundle degree in the spherical twist
    line_bundle_2 : int
        The second line bundle degree in the spherical twist

    Returns
    -------
    str
        A JSON string representation of the chain complex

    Raises
    ------
    ValueError
        If the input data is not an integer
    """
    
    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Input data must be integers")

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='P1'),
                          LineBundle(line_bundle_2, catagory='P1'))
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


if __name__ == "__main__":
    ints_to_mass_plot_P1_sing_twist(1, 2, return_json=False)
