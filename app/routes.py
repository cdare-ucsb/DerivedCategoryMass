from flask import Blueprint, request, jsonify, render_template
from src.CY3Fold import generate_plotly_graph_1  # Import the function if it's defined in app.utils
from  app.models import load_simple_model,  preprocess_input, predict
from src import LineBundle, SphericalTwist, LePotier
import plotly.graph_objects as go
import plotly.utils
import numpy as np
import json





bp = Blueprint('routes', __name__)

# Sample route for prediction
@bp.route('/predict', methods=['POST'])
def predict_route():
    if request.is_json:
        input_data = request.json.get('input_data')
    else:
        input_data = request.form.get('input_data')
        if input_data:
            input_data = [float(x) for x in input_data.split(',')]

    if not input_data:
        return jsonify({"error": "No input data provided"}), 400

    model = load_simple_model('app/models/trained_model.pth')
    processed_input = preprocess_input(input_data)
    prediction = predict(model, processed_input)

    return render_template('result.html', prediction=prediction.item())



@bp.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')
    

    

@bp.route('/local-P2', methods=['GET', 'POST'])
def LocalP2():
    
    DLP = LePotier(width=5, granularity=5)
    DLP_2d_json = DLP.plot_drezet_le_potier(return_json=True, show_walls=True)
    chamber_struct = DLP.plot_continuing_chamber(return_json=True)

    return render_template('local-P2.html', DLP_2d_json = DLP_2d_json, chamber_struct = chamber_struct)
    # return render_template('index.html')

@bp.route('/local-P1', methods=['GET', 'POST'])
def LocalP1():

    k3_plot_json = generate_plotly_graph_1(return_json=True)
    return render_template('local-P1.html', k3_plot_json = k3_plot_json)




@bp.route('/plot_sph_twist_P2', methods=['POST'])
def plot_sph_twist_P2():
    line_bundle_1 = request.form['line_bundle_1']
    line_bundle_2 = request.form['line_bundle_2']

    try:
        line_bundle_1 = int(line_bundle_1)
        line_bundle_2 = int(line_bundle_2)

        plot_json = compute_plot_data_P2(line_bundle_1, line_bundle_2)
        chain_complex_data = compute_chain_complex_data_P2(line_bundle_1, line_bundle_2)
        
        
        return render_template('plot_P2.html', plot_json=plot_json, chain_complex=chain_complex_data)


    except ValueError:
        return jsonify({"error": "Invalid input data"}), 400
    
  


def compute_plot_data_P2(line_bundle_1, line_bundle_2):

    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Input data must be integers")

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='P2'),
                          LineBundle(line_bundle_2, catagory='P2'))

    DLP = LePotier(width=5, granularity=3)
    

    # Define x values (spread around a region)
    x_vals = np.linspace(-5, 5, 200)  # X values from -2 to 2

    # Generate y values satisfying y > x^2
    y_vals = []
    for x in x_vals:
        y_min = float(DLP.curve_estimate(x)) # Slightly above x^2
        y_max = 25  # Arbitrary upper bound
        y_range = np.linspace(y_min, y_max, 160)  # 50 points per x value
        y_vals.append(y_range)

    # Convert to numpy array
    y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # Repeat x values to match the shape of y
    x_vals = np.repeat(x_vals, 160)  # Each x value repeats 10 times

    masses = np.array([sph.mass(x, y) for x, y in zip(x_vals, y_vals)])

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

    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)



def compute_chain_complex_data_P2(line_bundle_1, line_bundle_2):
    
    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Input data must be integers")

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='P2'),
                          LineBundle(line_bundle_2, catagory='P2'))
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


    


@bp.route('/plot_sph_twist_P1', methods=['POST'])
def plot_sph_twist_P1():
    line_bundle_1 = request.form['line_bundle_1']
    line_bundle_2 = request.form['line_bundle_2']

    try:
        line_bundle_1 = int(line_bundle_1)
        line_bundle_2 = int(line_bundle_2)

        plot_json = compute_plot_data_P1(line_bundle_1, line_bundle_2)
        chain_complex_data = compute_chain_complex_data_P1(line_bundle_1, line_bundle_2)
        
        
        return render_template('plot_P1.html', plot_json=plot_json, chain_complex=chain_complex_data)


    except ValueError:
        return jsonify({"error": "Invalid input data"}), 400
    

def compute_chain_complex_data_P1(line_bundle_1, line_bundle_2):
    
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


def compute_plot_data_P1(line_bundle_1, line_bundle_2):

    if not isinstance(line_bundle_1, int) or not isinstance(line_bundle_2, int):
        raise ValueError("Input data must be integers")

    sph = SphericalTwist(LineBundle(line_bundle_1, catagory='P1'),
                          LineBundle(line_bundle_2, catagory='P1'))


    # Define x values (spread around a region)
    x_vals = np.linspace(-5, 5, 200)  # X values from -2 to 2

    # Generate y values satisfying y > x^2
    y_vals = []
    for x in x_vals:
        y_range = np.linspace(0, 5, 160)  # 50 points per x value
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

    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)