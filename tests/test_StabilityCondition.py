import pytest
import os
import sys

from sympy import symbols, expand

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.ChernCharacter import ChernCharacter

from src.DerivedCategory.StabilityCondition import HNFactor, HarderNarasimhanError, HarderNarasimhanFiltration, StabilityCondition
from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.DerivedCategoryObject import GradedCoproductObject
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition



def test_StabilityCondition_init():
    H, C, D = symbols("H C D")
    basis3 = [H, C, D]
    tensor_data = {
        (H, H): 4,
        (C, C): -1,
        (D, D): -1,
        (H, C): 1,
        (C, D) : 0,
        (D, H) : 2
    }

    divisor_data = DivisorData(basis=basis3, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data)

    sigma1 = StabilityCondition(geometry_context, 0.5, 1.5)
    
    tensor_data_2 = {
        (H,): 1
    }
    divisor_data_2 = DivisorData(basis=[H], top_intersection_form=tensor_data_2)
    geometry_context2 = GeometryContext(catagory='P1', divisor_data=divisor_data_2)

    sigma2 = StabilityCondition(geometry_context2, complex(5,10))

    tensor_data_3 = {
        (H,H) : 1
    }
    divisor_data_3 = DivisorData(basis=[H], top_intersection_form=tensor_data_3)
    geometry_context3 = GeometryContext(catagory='P2', divisor_data=divisor_data_3)
    sigma3 = StabilityCondition(geometry_context3, 15.2,25)

    

    


    with pytest.raises(ValueError):
        StabilityCondition(geometry_context, 0.5, 1.5, 2.4)
    with pytest.raises(TypeError):
        StabilityCondition(geometry_context, "0.5", "2.4")
    with pytest.raises(TypeError):
        StabilityCondition(geometry_context2, 0.5)


def test_StabilityCondition_get_HN_factors():
    H, C, D = symbols("H C D")
    basis = [H]
    tensor_data = {
        (H, H): 4,
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data, polarization=H)

    sigma1 = StabilityCondition(geometry_context, 2, 2)
    
    lb1 = LineBundle(H, geometry_context=geometry_context)
    lb2 = LineBundle(2*H, geometry_context=geometry_context)
    lb3 = LineBundle(3*H, geometry_context=geometry_context)
    lb4 = LineBundle(4*H, geometry_context=geometry_context)
    lb5 = LineBundle(5*H, geometry_context=geometry_context)

    sph=SphericalTwistComposition([lb3, lb1])


    print("\n\n")
    print(sph.defining_triangle)
    print("\n\n")


    print(sigma1.get_HN_factors(sph))



    