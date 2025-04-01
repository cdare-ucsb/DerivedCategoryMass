import pytest
import os
import sys

from sympy import symbols

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.IntersectionForm import IntersectionForm


def test_IntersectionForm__init__():
    # Test with valid inputs
    x, y, z = symbols("x y z")
    basis = [x, y, z]
    dimension = 3
    tensor_data = {
        (x, x, x): 1,
        (x, x, y): 2,
        (x, y, z): -1,
        (x, y, y): 3,
        (y, y, y): 4
    }

    intersection_form = IntersectionForm(basis=basis, dimension=dimension, tensor_data=tensor_data)

    assert intersection_form.basis == basis
    assert intersection_form.dimension == dimension
    assert intersection_form.tensor[(x, y, z)] == -1  

    # Test with symmetrization
    assert intersection_form.tensor[(z,y,x)] == -1  # Symmetrized value
    assert intersection_form.tensor[(x, y, x)] == 2  # Symmetrized value

    # Test with invalid inputs
    with pytest.raises(TypeError):
        # Test with invalid basis
        IntersectionForm(basis="not_a_list", dimension=dimension, tensor_data=tensor_data)
    with pytest.raises(ValueError):
        # Test with invalid dimension
        IntersectionForm(basis=basis, dimension=-1, tensor_data=tensor_data)
    with pytest.raises(TypeError):
        # Test with invalid tensor data
        IntersectionForm(basis=basis, dimension=dimension, tensor_data="not_a_dict")



def test_IntersectionForm__evaluate():
    # Test with valid inputs
    x, y, z = symbols("x y z")
    basis = [x, y, z]
    dimension = 3
    tensor_data = {
        (x, x, x): 1,
        (x, x, y): 2,
        (x, y, z): -1,
        (x, y, y): 3,
        (y, y, y): 4
    }

    intersection_form = IntersectionForm(basis=basis, dimension=dimension, tensor_data=tensor_data)

    expr1 = x**2 + 2*x*y + 3*y**2
    expr2 = x + y + z
    result = intersection_form.evaluate(expr1, expr2)

    assert result == 1 + 2*2 + 3*3 + 2 + 2*3 + 3*4  + 0 -2  + 0

    # Test with invalid inputs
    with pytest.raises(TypeError):
        # Test with invalid expression
        intersection_form.evaluate("not_a_sympy_expression")

    with pytest.raises(ValueError):
        # Test with invalid expression
        intersection_form.evaluate(x**2 + y**2, expr1)