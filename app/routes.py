from flask import Blueprint, request, jsonify, render_template
# from app.models import load_simple_model,  preprocess_input, predict
from src.LocalP2 import LePotier, plot_successive_neighbors_ex1
from src.ProjectiveCY import K3GeometricChamber, complex_hypersurface_plotly_ex1
from src.SphericalTwist import SphericalTwist, DoubleSphericalTwist
from src.MassPlot import MassPlot
from src.CoherentSheaf import LineBundle

import traceback


bp = Blueprint('routes', __name__)


@bp.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')
    

    

@bp.route('/local-P2', methods=['GET', 'POST'])
def LocalP2():
    
    DLP = LePotier(width=5, granularity=5)
    DLP_2d_json = DLP.plot_region(return_json=True, show_walls=True)
    chamber_struct = plot_successive_neighbors_ex1(return_json=True)

    return render_template('local-P2.html', DLP_2d_json = DLP_2d_json, chamber_struct = chamber_struct)
    # return render_template('index.html')

@bp.route('/local-P1', methods=['GET', 'POST'])
def LocalP1():

    k3_plot_json = complex_hypersurface_plotly_ex1(degree=4, return_json=True)
    return render_template('local-P1.html', k3_plot_json = k3_plot_json)

@bp.route('/K3-surface', methods=['GET', 'POST'])
def K3Surface():
    K3 = K3GeometricChamber(degree=1) # Create double cover of P2.

    K3_alpha_beta_json = K3.plot_alpha_beta_plane(return_json=True)

    k3_plot_json = complex_hypersurface_plotly_ex1(degree=4, return_json=True)
    return render_template('K3-surface.html',
                            k3_plot_json = k3_plot_json,
                            K3_alpha_beta_json = K3_alpha_beta_json)




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



