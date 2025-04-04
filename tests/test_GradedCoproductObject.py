import pytest
import os
import sys

from sympy import symbols, expand

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.DerivedCategoryObject import GradedCoproductObject
from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData
from src.DerivedCategory.CoherentSheaf import LineBundle


@pytest.fixture
def coprod_K3_ex1():

    H3, C, D = symbols("H C D")
    basis3 = [H3, C, D]
    tensor_data = {
        (H3, H3): 4,
        (C, C): -1,
        (D, D): -1,
        (H3, C): 1,
        (C, D) : 0,
        (D, H3) : 2
    }

    divisor_data = DivisorData(basis=basis3, top_intersection_form=tensor_data)
    geometry_context3 = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(3*H3, geometry_context=geometry_context3)
    lb2 = LineBundle(4*C, geometry_context=geometry_context3)
    lb3 = LineBundle(5*D, geometry_context=geometry_context3)
    lb4 = LineBundle(6*H3, geometry_context=geometry_context3)
    lb5 = LineBundle(7*C, geometry_context=geometry_context3)

    return [lb1, lb2, lb3, lb4, lb5]

@pytest.fixture
def coprod_K3_ex2():

    H, C = symbols("H C")
    basis3 = [H, C]
    tensor_data = {
        (H, H) : 6,
        (C, C) : -1,
        (H, C) : 0,
    }

    divisor_data = DivisorData(basis=basis3, top_intersection_form=tensor_data)
    geometry_context3 = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(H, geometry_context=geometry_context3)
    lb2 = LineBundle(2*C, geometry_context=geometry_context3)

    return [lb1, lb2], H, C


def test_GradedCoproductObject_init(coprod_K3_ex1):
    

    # Test with valid inputs
    graded_coproduct = GradedCoproductObject(coprod_K3_ex1)


def test_GradedCoproductObject_cherncharacter(coprod_K3_ex2):


    [lb1, lb2], H, C = coprod_K3_ex2


    graded_coproduct = GradedCoproductObject([lb1, lb2])
    assert graded_coproduct.chernCharacter().expr == expand(2*C**2 + 2*C + H**2/2 + H + 2)


