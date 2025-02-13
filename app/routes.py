from flask import Blueprint, request, jsonify, render_template
import torch
from  app.models import load_simple_model, load_advanced_model, preprocess_input, predict

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
