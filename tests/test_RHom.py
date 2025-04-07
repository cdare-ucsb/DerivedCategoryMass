import pytest
import os
import sys

from sympy import symbols

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))



from src.DerivedCategory.DerivedCategoryObject import GradedCoproductObject
from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.RHom import RHom
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition


@pytest.fixture
def coprod_K3_ex1():

    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H, H) : 6
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(3*H, geometry_context=geometry_context)
    lb2 = LineBundle(2*H, geometry_context=geometry_context)

    coprod_lb1s = GradedCoproductObject(object_vector=[ lb1, lb1 ],
                                 shift_vector=[0, 3],
                                 dimension_vector=[3, 6])
    
    coprod_lb2s = GradedCoproductObject(object_vector=[ lb2, lb2 ],
                                 shift_vector=[0, 3],
                                 dimension_vector=[3, 6])
    
    

    return [lb1, lb2, coprod_lb1s, coprod_lb2s]


def test_RHom_lb_to_lb(coprod_K3_ex1):
    lb1, lb2, coprod_lb1s, coprod_lb2s = coprod_K3_ex1
    rhom1 = RHom(lb1, lb2)
    
    assert rhom1 is not None

    assert rhom1[-2] == 3*(-1)**2 + 2
    assert rhom1.get(0) is None
    
    rhom2 = RHom(lb2, lb1)
    assert rhom2 is not None
    assert rhom2[0] == 3*(1)**2 + 2
    assert rhom2.get(-2) is None    

def test_RHom_lb_to_graded_coproduct_object(coprod_K3_ex1):
    lb1, lb2, coprod_lb1s, coprod_lb2s = coprod_K3_ex1
    rhom = RHom(lb2, coprod_lb1s)
    
    assert rhom is not None
    
    assert rhom[0] == coprod_lb1s.dimension_vector[0]*RHom(lb2, lb1)[0]
    assert rhom[3] == coprod_lb1s.dimension_vector[1]*RHom(lb2, lb1)[0]

    rhom2 = RHom(lb1, coprod_lb2s)

    assert rhom2 is not None

    assert rhom2[-2] == coprod_lb2s.dimension_vector[0]*RHom(lb1, lb2)[-2]
    assert rhom2[1] == coprod_lb2s.dimension_vector[1]*RHom(lb1, lb2)[-2]


def test_RHom_lb_to_sph_twist():

    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H, H) : 6
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(7*H, geometry_context=geometry_context)
    lb2 = LineBundle(5*H, geometry_context=geometry_context)
    lb3 = LineBundle(3*H, geometry_context=geometry_context)

    sph = SphericalTwistComposition(line_bundle_vector=[lb1, lb2, lb3])

    assert sph.defining_triangle[0].dimension_vector[0] ==  RHom(lb3, lb2)[0] * 14 - RHom(lb3, lb1)[0]

def test_RHom_lb_to_double_sph_twist():

    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H, H) : 6
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(H, geometry_context=geometry_context)
    lb2 = LineBundle(3*H, geometry_context=geometry_context)
    lb3 = LineBundle(5*H, geometry_context=geometry_context)
    lb4 = LineBundle(7*H, geometry_context=geometry_context)
    

    assert RHom(lb1, SphericalTwistComposition(line_bundle_vector=[lb4, lb3, lb2]))[2] == RHom(lb1, lb2)[0] * 146 - RHom(lb1, SphericalTwistComposition([lb4, lb3]) )[1]

