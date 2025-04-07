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
def single_twist_K3_ex1():

    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H, H) : 6
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(3*H, geometry_context=geometry_context)
    lb2 = LineBundle(5*H, geometry_context=geometry_context)


    return [lb1, lb2]


@pytest.fixture
def double_twist_K3_ex1():
    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H, H) : 6
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(3*H, geometry_context=geometry_context)
    lb2 = LineBundle(5*H, geometry_context=geometry_context)
    lb3 = LineBundle(7*H, geometry_context=geometry_context)

    return [lb1, lb2, lb3]

@pytest.fixture
def triple_twist_K3_ex1():

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

    return [lb1, lb2, lb3, lb4]




def test_SphericalTwistComposition_init_ex1(single_twist_K3_ex1):

    lb1, lb2 = single_twist_K3_ex1

    sph1 = SphericalTwistComposition(line_bundle_vector=[lb2, lb1])

    H = sph1.geometry_context.divisor_data.basis[0]

    # There should only be one instance of the SphericalTwistComposition
    # with the same line bundles
    assert sph1 is SphericalTwistComposition(line_bundle_vector=[ LineBundle(5*H, sph1.geometry_context),
                                                                 LineBundle(3*H, sph1.geometry_context)])
    
    assert sph1 is not SphericalTwistComposition(line_bundle_vector=[ LineBundle(3*H, sph1.geometry_context),
                                                                    LineBundle(5*H, sph1.geometry_context)])
    

def test_SphericalTwistComposition_init_ex2(double_twist_K3_ex1):

    lb1, lb2, lb3 = double_twist_K3_ex1

    sph2 = SphericalTwistComposition(line_bundle_vector=[lb3, lb2, lb1])

    H = sph2.geometry_context.divisor_data.basis[0]

    assert sph2 is SphericalTwistComposition(line_bundle_vector=[ LineBundle(7*H, sph2.geometry_context),
                                                                LineBundle(5*H, sph2.geometry_context),
                                                                 LineBundle(3*H, sph2.geometry_context)])

    assert sph2 is not SphericalTwistComposition(line_bundle_vector=[ LineBundle(5*H, sph2.geometry_context),
                                                                LineBundle(7*H, sph2.geometry_context),
                                                                 LineBundle(3*H, sph2.geometry_context)])
    
    assert sph2 is not SphericalTwistComposition(line_bundle_vector=[ LineBundle(7*H, sph2.geometry_context),
                                                                LineBundle(5*H, sph2.geometry_context)])


    
def test_SphericalTwistComposition_get_defining_triangle_ex1(single_twist_K3_ex1):
    lb1, lb2 = single_twist_K3_ex1

    sph1 = SphericalTwistComposition(line_bundle_vector=[lb2, lb1])

    dim_hom = RHom(lb1, lb2)

    assert isinstance(sph1.defining_triangle[0], GradedCoproductObject)
    assert isinstance(sph1.defining_triangle[1], LineBundle)
    assert isinstance(sph1.defining_triangle[2], SphericalTwistComposition)

    assert sph1.defining_triangle[0].dimension_vector[0] == dim_hom[0]
    assert sph1.defining_triangle[0].object_vector[0] is lb1
    assert sph1.defining_triangle[0].shift_vector[0] == 0


def test_SphericalTwistComposition_get_defining_triangle_ex2(double_twist_K3_ex1):

    lb1, lb2, lb3 = double_twist_K3_ex1

    sph2 = SphericalTwistComposition(line_bundle_vector=[lb3, lb2, lb1])

    assert isinstance(sph2.defining_triangle[0], GradedCoproductObject)
    assert isinstance(sph2.defining_triangle[1], SphericalTwistComposition)
    assert isinstance(sph2.defining_triangle[2], SphericalTwistComposition)
    
    assert isinstance(sph2.defining_triangle[0].object_vector[0], LineBundle)
    assert sph2.defining_triangle[0].dimension_vector[0] == 146
    assert sph2.defining_triangle[0].shift_vector[0] == 1


def test_SphericalTwistComposition_get_defining_triangle_ex3(triple_twist_K3_ex1):

    lb1, lb2, lb3, lb4 = triple_twist_K3_ex1

    sph3 = SphericalTwistComposition(line_bundle_vector=[lb4, lb3, lb2, lb1])

    print("\n" + str(sph3.defining_triangle))




def test_SphericalTwistComposition_get_canonical_triangles_ex1(single_twist_K3_ex1):
    lb1, lb2 = single_twist_K3_ex1

    sph1 = SphericalTwistComposition(line_bundle_vector=[lb2, lb1])

    assert sph1.canonical_triangles() == [sph1.defining_triangle]


def test_SphericalTwistComposition_get_canonical_triangles_ex2(double_twist_K3_ex1):
    
    lb1, lb2, lb3 = double_twist_K3_ex1

    sph2 = SphericalTwistComposition(line_bundle_vector=[lb3, lb2, lb1])

    canonical_triangles = sph2.canonical_triangles()
    assert len(canonical_triangles) == 2
    assert canonical_triangles[0].object1.dimension_vector[0] == 146
    assert canonical_triangles[0].object1.object_vector[0] is lb1
    assert canonical_triangles[0].object1.shift_vector[0] == 1

    assert canonical_triangles[0].object2.line_bundle_vector[1] is lb2
    assert canonical_triangles[0].object2.line_bundle_vector[0] is lb3

    assert canonical_triangles[1].object1.dimension_vector[0] == 14

    expected_first_obj = GradedCoproductObject(object_vector=[SphericalTwistComposition(line_bundle_vector=[lb2, lb1])], shift_vector=[0], dimension_vector=[14])

    assert canonical_triangles[1].object1 is expected_first_obj
