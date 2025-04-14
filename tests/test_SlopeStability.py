import os
import sys

from sympy import symbols, expand

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.StabilityCondition import HNFactor, HarderNarasimhanError, HarderNarasimhanFiltration, StabilityCondition, SlopeStability
from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.DerivedCategoryObject import GradedCoproductObject
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition



def test_SlopeStability_init():

    H, C = symbols('H C')
    divisor_data = DivisorData(basis=[H,C] ,
                            top_intersection_form={(H,H): 2,
                                                   (C,C):-2,
                                                   (H,C):0})
    geometry_context = GeometryContext(catagory="K3", divisor_data=divisor_data)
    slope_stability = SlopeStability(geometry_context)
    assert slope_stability.geometry_context == geometry_context
    assert slope_stability._hn_filt_cache == {}

def test_SlopeStability_line_bundle_destabilizers():
    H, C= symbols('H C')
    divisor_data = DivisorData(basis=[H,C] ,
                            top_intersection_form={(H,H): 2,
                                                   (C,C):-2,
                                                   (H,C):-1})
    geometry_context = GeometryContext(catagory="K3", divisor_data=divisor_data, polarization=H)
    slope_stability = SlopeStability(geometry_context)

    L = LineBundle(10*C+10*H, geometry_context)
    
    destabilizers = list(slope_stability.line_bundle_destabilizers(L))

    
    assert len(destabilizers) == 1
    assert destabilizers[0].divisor == 9*C + 10*H






    

def test_get_HN_factors():


    H, C, D = symbols("H C D")
    basis = [H, C, D]
    tensor_data = {
        (H, H): 4,
        (C,C):-2,
        (H,C):1,
        (D,D):-2,
        (H,D):2,
        (C,D):0
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data, polarization=H)

    C_coeff = 48
    D_coeff = 93
    H_coeff = 97
    lb1 = LineBundle(C_coeff*C + D_coeff*D + H_coeff*H, geometry_context=geometry_context)

    slope_stab = SlopeStability(geometry_context)
    lb1_slope_HN = slope_stab.get_HN_factors(lb1)
    print(lb1_slope_HN)

