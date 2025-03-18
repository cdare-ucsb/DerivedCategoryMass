from flask import Blueprint, request, jsonify, render_template, send_file
from app import socketio


# from app.models import load_simple_model,  preprocess_input, predict
from src.LocalP2 import LePotier, plot_successive_neighbors_ex1
from src.ProjectiveCY import K3GeometricChamber, complex_hypersurface_plotly_ex1
from src.SphericalTwist import SphericalTwist, DoubleSphericalTwist
from src.MassPlot import MassPlot
from src.CoherentSheaf import LineBundle
from src.model import SingleTwistModel

import traceback
import json
import plotly

import plotly.utils
import io
import os
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F  # Correct import for relu, softmax, etc.


# Temporary storage for model files (can be replaced with Redis or a database)
model_files = {}
stm = None

UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)  # Ensure upload directory exists


bp = Blueprint('routes', __name__)




@socketio.on("refresh")
def handle_refresh():
    print("\n\n\t\033[94mRefreshing page for all clients!\033[0m\n\n")
    socketio.emit("reload")




@bp.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')
    






##############################################
#              Local P1 Routes               #
##############################################     

@socketio.on('connect',  namespace='/local-P1')
def handle_connect():
    print("\n\n\t\033[94mClient connected to Local P1 page\033[0m\n\n")

@socketio.on('disconnect', namespace='/local-P1')
def handle_disconnect():
    print("\n\n\t\033[94mClient disconnected from Local P1 page\033[0m\n\n")




@bp.route('/local-P1', methods=['GET', 'POST'])
def LocalP1():

    k3_plot_json = complex_hypersurface_plotly_ex1(degree=4, return_json=True)
    return render_template('local-P1.html', k3_plot_json = k3_plot_json)


@bp.route('/plot_sph_twist_P1', methods=['POST'])
def plot_sph_twist_P1():
    line_bundle_1 = request.form['line_bundle_1']
    line_bundle_2 = request.form['line_bundle_2']

    try:
        line_bundle_1 = int(line_bundle_1)
        line_bundle_2 = int(line_bundle_2)

        sph = SphericalTwist(LineBundle(line_bundle_1, catagory='P1'),
                             LineBundle(line_bundle_2, catagory='P1'),
                             degree=1)
        
        mp = MassPlot(line_bundle_1=line_bundle_1,
                      line_bundle_2=line_bundle_2,
                      catagory='P1',
                      degree=1)

        plot_json = mp.to_plotly_json() 
        chain_complex_data = sph.defining_triangle_to_json()
        
        
        return render_template('plot_P1.html',
                                plot_json=plot_json,
                                chain_complex=chain_complex_data)

    except ValueError:
        print(traceback.format_exc())
        return jsonify({"error": "Invalid input data"}), 400
    

@socketio.on('train-model-P1', namespace='/local-P1')
def train_model_P1(data):
    """
    Handles the long-running task and emits progress updates.
    """
    line_bundle_1 = data.get("line_bundle_1")
    line_bundle_2 = data.get("line_bundle_2")
    filename = data.get("filename")

    print(f"\n\n\t\033[94mReceived request to process: {line_bundle_1}, {line_bundle_2}, {filename}\033[0m\n\n")

    stm = SingleTwistModel(line_bundle_1=line_bundle_1,
                           line_bundle_2=line_bundle_2,
                           catagory='P1',
                           degree=1, mode='disc')
    
    num_epochs = 5000
    
    # Define loss function and optimizer
    criterion = nn.MSELoss()  # Mean Squared Error loss
    optimizer = optim.Adam(stm.model.parameters(), lr=0.01)

    # Use learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                    factor=0.5, patience=20,
                                                    verbose=True)

    for epoch in range(num_epochs):
        # Training loop
        
        optimizer.zero_grad()
        
        # Forward pass
        pred = stm.model(stm.input_tensor)
        
        # Compute loss
        loss = criterion(pred, stm.output_tensor.view(-1, 1))
        
        # Backward pass
        loss.backward()
        optimizer.step()

        # Step the scheduler
        scheduler.step(loss)

        socketio.emit('progress', {'progress': int((100*epoch + 100) / num_epochs)},
                    namespace='/local-P1')  # Send progress update


    # Save model to memory instead of disk
    buffer = io.BytesIO()
    torch.save(stm.model.state_dict(), buffer)  # Save directly into memory
    buffer.seek(0)  # Move cursor to the beginning

    # Store the file in a temporary dictionary (you can use Redis instead)
    model_files[filename] = buffer.getvalue()  # Store as binary data

    print(f"\n\n\t\033[94mModel is ready for download: {filename}!\033[0m\n\n")

    # Notify the frontend that the file is ready
    socketio.emit("download_ready_P1", {"filename": filename}, namespace="/local-P1")

    
@bp.route('/download-model-P1/<filename>')
def download_model_P1(filename):
    """
    Serves the saved model file when requested.
    """

    if filename not in model_files:
        return "File not found!", 404

    buffer = io.BytesIO(model_files[filename])  # Retrieve the file from memory
    buffer.seek(0)

    return send_file(
        buffer,
        as_attachment=True,
        download_name=filename,
        mimetype="application/octet-stream"
    )

@bp.route('/upload_P1', methods=['POST'])
def upload_file_P1():
    """
    Handles file upload, reads x, y data, and triggers graph update.
    """
    if 'file' not in request.files:
        return "No file uploaded", 400

    file = request.files['file']
    
    if file.filename == '':
        return "No selected file", 400

    # Read uploaded file (assumes CSV format with 'x' and 'y' columns)
    try:
        stm = SingleTwistModel(line_bundle_1=1,
                               line_bundle_2=1,
                               catagory='P1',
                               degree=1, mode='disc')
        stm.model.load_state_dict(torch.load(file))
    except Exception as e:
        return f"Error reading file: {str(e)}", 400

    # Generate Plotly Graph
    fig = stm.predictions_to_plotly()

    #  Convert Figure to JSON
    disc_graph_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    # Send JSON to frontend using WebSockets
    socketio.emit("plot_update", {"disc_graph_json": disc_graph_json}, namespace="/local-P1")

    return "File uploaded successfully", 200











##############################################
#              Local P2 Routes               #
##############################################  


@socketio.on('connect',  namespace='/local-P2')
def handle_connect():
    print("\n\n\t\033[94mClient connected to Local P2 page\033[0m\n\n")

@socketio.on('disconnect', namespace='/local-P2')
def handle_disconnect():
    print("\n\n\t\033[94mClient disconnected from Local P2 page\033[0m\n\n")




    

@bp.route('/local-P2', methods=['GET', 'POST'])
def LocalP2():
    
    DLP = LePotier(width=5, granularity=5)
    DLP_2d_json = DLP.plot_region(return_json=True, show_walls=True)
    chamber_struct = plot_successive_neighbors_ex1(return_json=True)

    return render_template('local-P2.html',
                            DLP_2d_json = DLP_2d_json,
                            chamber_struct = chamber_struct)



@bp.route('/plot_sph_twist_P2', methods=['POST'])
def plot_sph_twist_P2():
    line_bundle_1 = request.form['line_bundle_1']
    line_bundle_2 = request.form['line_bundle_2']

    try:
        line_bundle_1 = int(line_bundle_1)
        line_bundle_2 = int(line_bundle_2)

        sph = SphericalTwist(LineBundle(line_bundle_1, catagory='P2'),
                             LineBundle(line_bundle_2, catagory='P2'),
                             degree=1)
        
        mp = MassPlot(line_bundle_1=line_bundle_1,
                      line_bundle_2=line_bundle_2,
                      catagory='P2',
                      degree=1)

        plot_json = mp.to_plotly_json()
        chain_complex_data = sph.defining_triangle_to_json()
        
        
        return render_template('plot_P2.html', plot_json=plot_json, chain_complex=chain_complex_data)


    except ValueError:
        print(traceback.format_exc())
        return jsonify({"error": "Invalid input data"}), 400
    

@socketio.on('train-model-P2', namespace='/local-P2')
def train_model_P2(data):
    """
    Handles the long-running task and emits progress updates.
    """
    line_bundle_1 = data.get("line_bundle_1")
    line_bundle_2 = data.get("line_bundle_2")
    filename = data.get("filename")

    print(f"\n\n\t\033[94mReceived request to process: {line_bundle_1}, {line_bundle_2}, {filename}\033[0m\n\n")

    stm = SingleTwistModel(line_bundle_1=line_bundle_1,
                           line_bundle_2=line_bundle_2,
                           catagory='P2',
                           degree=1, mode='disc')
    
    num_epochs = 5000
    
    # Define loss function and optimizer
    criterion = nn.MSELoss()  # Mean Squared Error loss
    optimizer = optim.Adam(stm.model.parameters(), lr=0.01)

    # Use learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                    factor=0.5, patience=20,
                                                    verbose=True)

    for epoch in range(num_epochs):
        # Training loop
        
        optimizer.zero_grad()
        
        # Forward pass
        pred = stm.model(stm.input_tensor)
        
        # Compute loss
        loss = criterion(pred, stm.output_tensor.view(-1, 1))
        
        # Backward pass
        loss.backward()
        optimizer.step()

        # Step the scheduler
        scheduler.step(loss)

        socketio.emit('progress', {'progress': int((100*epoch + 100) / num_epochs)},
                    namespace='/local-P2')  # Send progress update


    # Save model to memory instead of disk
    buffer = io.BytesIO()
    torch.save(stm.model.state_dict(), buffer)  # Save directly into memory
    buffer.seek(0)  # Move cursor to the beginning

    # Store the file in a temporary dictionary (you can use Redis instead)
    model_files[filename] = buffer.getvalue()  # Store as binary data

    print(f"\n\n\t\033[94mModel is ready for download: {filename}!\033[0m\n\n")

    # Notify the frontend that the file is ready
    socketio.emit("download_ready_P2", {"filename": filename}, namespace="/local-P2")

    
@bp.route('/download-model-P2/<filename>')
def download_model_P2(filename):
    """
    Serves the saved model file when requested.
    """

    if filename not in model_files:
        return "File not found!", 404

    buffer = io.BytesIO(model_files[filename])  # Retrieve the file from memory
    buffer.seek(0)

    return send_file(
        buffer,
        as_attachment=True,
        download_name=filename,
        mimetype="application/octet-stream"
    )

@bp.route('/upload_P2', methods=['POST'])
def upload_file_P2():
    """
    Handles file upload, reads x, y data, and triggers graph update.
    """
    if 'file' not in request.files:
        return "No file uploaded", 400

    file = request.files['file']
    
    if file.filename == '':
        return "No selected file", 400

    # Read uploaded file (assumes CSV format with 'x' and 'y' columns)
    try:
        stm = SingleTwistModel(line_bundle_1=1,
                               line_bundle_2=1,
                               catagory='P2',
                               degree=1, mode='disc')
        stm.model.load_state_dict(torch.load(file))
    except Exception as e:
        return f"Error reading file: {str(e)}", 400

    # Generate Plotly Graph
    fig = stm.predictions_to_plotly()

    #  Convert Figure to JSON
    disc_graph_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    # Send JSON to frontend using WebSockets
    socketio.emit("plot_update", {"disc_graph_json": disc_graph_json}, namespace="/local-P2")

    return "File uploaded successfully", 200













##############################################   
#              K3 Surface Routes             #
##############################################   


@socketio.on('connect',  namespace='/K3-surface')
def handle_connect():
    print("\n\n\t\033[94mClient connected to K3 Surface page\033[0m\n\n")

@socketio.on('disconnect', namespace='/K3-surface')
def handle_disconnect():
    print("\n\n\t\033[94mClient disconnected from K3 Surface page\033[0m\n\n")




@bp.route('/K3-surface', methods=['GET', 'POST'])
def K3Surface():
    K3 = K3GeometricChamber(degree=1) # Create double cover of P2.

    K3_alpha_beta_json = K3.plot_alpha_beta_plane(return_json=True)

    k3_plot_json = complex_hypersurface_plotly_ex1(degree=4, return_json=True)
    return render_template('K3-surface.html',
                            k3_plot_json = k3_plot_json,
                            K3_alpha_beta_json = K3_alpha_beta_json)
    

@bp.route('/plot_sing_sph_twist_K3', methods=['POST'])
def plot_sing_sph_twist_K3():
    line_bundle_1 = request.form['line_bundle_1']
    line_bundle_2 = request.form['line_bundle_2']

    try:
        line_bundle_1 = int(line_bundle_1)
        line_bundle_2 = int(line_bundle_2)

        sph = SphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                             LineBundle(line_bundle_2, catagory='K3'),
                             degree=1)
        
        mp = MassPlot(line_bundle_1=line_bundle_1,
                      line_bundle_2=line_bundle_2,
                      catagory='K3',
                      degree=1)

        plot_json = mp.to_plotly_json()
        chain_complex_data = sph.defining_triangle_to_json()
        
        return render_template('plot_K3_sing.html',
                                plot_json=plot_json,
                                chain_complex=chain_complex_data)
    

    except ValueError:
        print(traceback.format_exc())
        return jsonify({"error": "Invalid input data"}), 400
    
@bp.route('/plot_double_sph_twist_K3', methods=['POST'])
def plot_double_sph_twist_K3():
    line_bundle_1 = request.form['line_bundle_1']
    line_bundle_2 = request.form['line_bundle_2']
    line_bundle_3 = request.form['line_bundle_3']

    try:
        line_bundle_1 = int(line_bundle_1)
        line_bundle_2 = int(line_bundle_2)
        line_bundle_3 = int(line_bundle_3)

        sph = DoubleSphericalTwist(LineBundle(line_bundle_1, catagory='K3'),
                                    LineBundle(line_bundle_2, catagory='K3'),
                                    LineBundle(line_bundle_3, catagory='K3'),
                                    degree=1)

        mp = MassPlot(line_bundle_1=line_bundle_1,
                      line_bundle_2=line_bundle_2,
                      line_bundle_3=line_bundle_3,
                      catagory='K3',
                      degree=1)

        plot_json = mp.to_plotly_json()
        chain_complex_data = sph.defining_triangle_to_json()
        secondary_complex_data = sph.secondary_triangle_to_json()
        
        
        return render_template('plot_K3_double.html',
                                plot_json=plot_json,
                                chain_complex=chain_complex_data,
                                secondary_complex=secondary_complex_data)

    except ValueError:
        print(traceback.format_exc())
        return jsonify({"error": "Invalid input data"}), 400



@socketio.on('train-model-K3', namespace='/K3-surface')
def train_model_K3(data):
    """
    Handles the long-running task and emits progress updates.
    """
    line_bundle_1 = data.get("line_bundle_1")
    line_bundle_2 = data.get("line_bundle_2")
    filename = data.get("filename")

    print(f"\n\n\t\033[94mReceived request to process: {line_bundle_1}, {line_bundle_2}, {filename}\033[0m\n\n")

    stm = SingleTwistModel(line_bundle_1=line_bundle_1,
                           line_bundle_2=line_bundle_2,
                           catagory='K3',
                           degree=1, mode='disc')
    
    num_epochs = 5000
    
    # Define loss function and optimizer
    criterion = nn.MSELoss()  # Mean Squared Error loss
    optimizer = optim.Adam(stm.model.parameters(), lr=0.01)

    # Use learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                    factor=0.5, patience=20,
                                                    verbose=True)

    for epoch in range(num_epochs):
        # Training loop
        
        optimizer.zero_grad()
        
        # Forward pass
        pred = stm.model(stm.input_tensor)
        
        # Compute loss
        loss = criterion(pred, stm.output_tensor.view(-1, 1))
        
        # Backward pass
        loss.backward()
        optimizer.step()

        # Step the scheduler
        scheduler.step(loss)

        socketio.emit('progress', {'progress': int((100*epoch + 100) / num_epochs)},
                    namespace='/K3-surface')  # Send progress update


    # Save model to memory instead of disk
    buffer = io.BytesIO()
    torch.save(stm.model.state_dict(), buffer)  # Save directly into memory
    buffer.seek(0)  # Move cursor to the beginning

    # Store the file in a temporary dictionary (you can use Redis instead)
    model_files[filename] = buffer.getvalue()  # Store as binary data

    print(f"\n\n\t\033[94mModel is ready for download: {filename}!\033[0m\n\n")

    # Notify the frontend that the file is ready
    socketio.emit("download_ready_K3", {"filename": filename}, namespace="/K3-surface")


@bp.route('/download-model-K3/<filename>')
def download_model_K3(filename):
    """
    Serves the saved model file when requested.
    """

    if filename not in model_files:
        return "File not found!", 404

    buffer = io.BytesIO(model_files[filename])  # Retrieve the file from memory
    buffer.seek(0)

    return send_file(
        buffer,
        as_attachment=True,
        download_name=filename,
        mimetype="application/octet-stream"
    )


@bp.route('/upload_K3', methods=['POST'])
def upload_file_K3():
    """
    Handles file upload, reads x, y data, and triggers graph update.
    """
    if 'file' not in request.files:
        return "No file uploaded", 400

    file = request.files['file']
    
    if file.filename == '':
        return "No selected file", 400

    # Read uploaded file (assumes CSV format with 'x' and 'y' columns)
    try:
        stm = SingleTwistModel(line_bundle_1=1,
                               line_bundle_2=1,
                               catagory='K3',
                               degree=1, mode='disc')
        stm.model.load_state_dict(torch.load(file))
    except Exception as e:
        return f"Error reading file: {str(e)}", 400

    # Generate Plotly Graph
    fig = stm.predictions_to_plotly()

    #  Convert Figure to JSON
    disc_graph_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    # Send JSON to frontend using WebSockets
    socketio.emit("plot_update", {"disc_graph_json": disc_graph_json}, namespace="/K3-surface")

    return "File uploaded successfully", 200

    