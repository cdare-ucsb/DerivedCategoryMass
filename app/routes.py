from flask import Blueprint, request, jsonify, render_template
from  app.models import load_simple_model,  preprocess_input, predict
from src import LineBundle, SphericalTwist
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
    prediction = None
    if request.method == 'POST':
        input_data = request.form['input_data']
        input_list = [float(x) for x in input_data.split(',')]
        model = load_simple_model()
        processed_input = preprocess_input(input_list)
        prediction = predict(model, processed_input).item()  # Convert prediction to scalar value

    return render_template('index.html', prediction=prediction)
    # return render_template('index.html')



@bp.route('/plot_sph_twist', methods=['POST'])
def plot_sph_twist():
    line_bundle_1 = request.form['line_bundle_1']
    line_bundle_2 = request.form['line_bundle_2']

    try:
        line_bundle_1 = int(line_bundle_1)
        line_bundle_2 = int(line_bundle_2)

        sph = SphericalTwist(LineBundle(line_bundle_1), LineBundle(line_bundle_2))

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

        plot_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        
        return render_template('plot.html', plot_json=plot_json)


    except ValueError:
        return jsonify({"error": "Invalid input data"}), 400
    
    
