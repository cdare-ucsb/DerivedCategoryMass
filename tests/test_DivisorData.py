import pytest
import os
import sys

from sympy import symbols, expand

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.GeometryContext import DivisorData


def test_DivisorData__init__():
    # Test with valid inputs
    x, y, z = symbols("x y z")
    basis = [x, y, z]
    tensor_data = {
        (x, x, x): 1,
        (x, x, y): 2,
        (x, y, z): -1,
        (x, y, y): 3,
        (y, y, y): 4
    }

    divisor_data = DivisorData(basis=basis,  top_intersection_form=tensor_data)

    assert divisor_data.basis == basis
    assert divisor_data.variety_dimension == 3
    assert divisor_data.top_intersection_form[(x, y, z)] == -1  

    # Test with symmetrization
    assert divisor_data.top_intersection_form[(z,y,x)] == -1  # Symmetrized value
    assert divisor_data.top_intersection_form[(x, y, x)] == 2  # Symmetrized value

    # Test with invalid inputs
    with pytest.raises(TypeError):
        # Test with invalid basis
        DivisorData(basis="not_a_list", top_intersection_form=tensor_data)
    with pytest.raises(TypeError):
        # Test with invalid tensor data
        DivisorData(basis=basis, top_intersection_form="not_a_dict")



def test_DivisorData__evaluate():
    # Test with valid inputs
    x, y, z = symbols("x y z")
    basis = [x, y, z]
    tensor_data = {
        (x, x, x): 1,
        (x, x, y): 2,
        (x, y, z): -1,
        (x, y, y): 3,
        (y, y, y): 4
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)

    expr1 = x**2 + 2*x*y + 3*y**2
    expr2 = x + y + z
    result = divisor_data.evaluate(expr1, expr2)

    assert result == 1 + 2*2 + 3*3 + 2 + 2*3 + 3*4  + 0 -2  + 0



    H = symbols("H")
    basis_2 = [H]
    tensor_data_2 = {
        (H, H): 2,
    }

    term1 = 1 + 3*H + 4.5*H**2
    term2 = -1 + 1j*H+ (10 - 20*1j) * H**2

    print(expand(term1*term2))

    divisor_data_2 = DivisorData(basis=basis_2, top_intersection_form=tensor_data_2)
    # result = divisor_data_2.evaluate(1 + 3*H + 4.5*H**2, -1 + 1j*H+ (10 - 20*1j) * H**2)
    # print(result)


    # Test with invalid inputs
    with pytest.raises(TypeError):
        # Test with invalid expression
        divisor_data.evaluate("not_a_sympy_expression")

    with pytest.raises(ValueError):
        # Test with invalid expression
        divisor_data.evaluate(x**2 + y**2, expr1)