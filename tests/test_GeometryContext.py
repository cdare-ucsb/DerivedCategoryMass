import pytest
import os
import sys

from sympy import symbols

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData


def test_GeometryContext_init():


    # Test with valid inputs
    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H,) : 1
    }


    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='P1', divisor_data=divisor_data)
    assert geometry_context.catagory == 'P1'
    assert geometry_context.variety_dimension() == 1

    H2 = symbols("H2")
    basis2 = [H2]
    tensor_data_2 = {
        (H2,H2) : 1
    }

    divisor_data = DivisorData(basis=basis2, top_intersection_form=tensor_data_2)
    geometry_context2 = GeometryContext(catagory='P2', divisor_data=divisor_data)
    assert geometry_context2.catagory == 'P2'
    assert geometry_context2.variety_dimension() == 2

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
    assert geometry_context3.catagory == 'K3'
    assert geometry_context3.variety_dimension() == 2
    assert geometry_context3.divisor_data.top_intersection_form[(H3, H3)] == 4
    assert geometry_context3.divisor_data.evaluate(3*H3**2) == 12


    with pytest.raises(NotImplementedError):
        # Test with invalid catagory
        GeometryContext(catagory='not_implemented', divisor_data=divisor_data)
    with pytest.raises(TypeError):
        # Test with invalid divisor_data
        GeometryContext(catagory='P1', divisor_data="not_a_divisor_data")
    with pytest.raises(ValueError):
        # Test with invalid divisor_data
        GeometryContext(catagory='P1', divisor_data=divisor_data, polarization="not_a_symbol")