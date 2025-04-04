import pytest
import os
import sys

from sympy import symbols, expand, cos


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData



def test_CoherentSheaf_init():
    # Test the initialization of the CoherentSheaf class
    
    # Test with valid inputs
    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H,) : 1
    }


    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='P1', divisor_data=divisor_data)

    lb1 = LineBundle(3*H, geometry_context=geometry_context)

    assert lb1.divisor == 3*H
    assert lb1.geometry_context.divisor_data.variety_dimension == 1

    assert lb1.chern_character[1] == expand(3*H)

    H2 = symbols("H2")
    basis2 = [H2]
    tensor_data_2 = {
        (H2,H2) : 1
    }

    divisor_data = DivisorData(basis=basis2, top_intersection_form=tensor_data_2)
    geometry_context2 = GeometryContext(catagory='P2', divisor_data=divisor_data)

    lb2 = LineBundle(4*H2, geometry_context=geometry_context2)

    assert lb2.divisor == 4*H2
    assert lb2.geometry_context.divisor_data.variety_dimension == 2
    assert lb2.chern_character[2] == expand(8*H2**2)

    H3, C, D = symbols("H3 C D")
    basis3 = [H3, C, D]
    tensor_data_3 = {
        (H3, H3): 4,
        (C, C): -1,
        (D, D): -1,
        (H3, C): 1,
        (C, D) : 0,
        (D, H3) : 2
    }

    divisor_data = DivisorData(basis=basis3, top_intersection_form=tensor_data_3)
    geometry_context3 = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb4 = LineBundle(3*H3 + 4*C, geometry_context=geometry_context3)
    assert lb4.divisor == 3*H3 + 4*C
    assert lb4.geometry_context.divisor_data.variety_dimension == 2
    assert lb4.chern_character[2] == expand(8*C**2 + 12*C*H3 + 9*H3**2/2)



    # Test improper initialization
    with pytest.raises(ValueError):
        LineBundle(3*H2, geometry_context=geometry_context)

    with pytest.raises(ValueError):
        LineBundle(cos(H), geometry_context=geometry_context)






# def test_CoherentSheaf_get_HN_factors():


    # H3, C, D = symbols("H3 C D")
    # basis3 = [H3, C, D]
    # tensor_data_3 = {
    #     (H3, H3): 4,
    #     (C, C): -1,
    #     (D, D): -1,
    #     (H3, C): 1,
    #     (C, D) : 0,
    #     (D, H3) : 2
    # }

    # divisor_data = DivisorData(basis=basis3, top_intersection_form=tensor_data_3)
    # geometry_context3 = GeometryContext(catagory='K3', divisor_data=divisor_data)

#     lb1 = LineBundle(3*H, geometry_context=geometry_context3)

#     assert lb1.get_HN_factors(1,2,3) == [lb1]

#     lb2 = LineBundle(3*H + 4*C, geometry_context=geometry_context3)

#     assert lb2.get_HN_factors(1,2,3) == [lb2]