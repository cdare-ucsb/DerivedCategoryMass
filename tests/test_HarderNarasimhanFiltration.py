import pytest
import os
import sys

from sympy import symbols, expand

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.ChernCharacter import ChernCharacter

from src.DerivedCategory.StabilityCondition import HNFactor, HarderNarasimhanError, HarderNarasimhanFiltration
from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.DerivedCategoryObject import GradedCoproductObject



def test_HNFilt_init():
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
    geometry_context3 = GeometryContext(catagory='K3', divisor_data=divisor_data)

    lb1 = LineBundle(3*H, geometry_context=geometry_context3)
    lb2 = LineBundle(4*C, geometry_context=geometry_context3)
    lb3 = LineBundle(5*D, geometry_context=geometry_context3)
    lb4 = LineBundle(6*H, geometry_context=geometry_context3)
    lb5 = LineBundle(7*C, geometry_context=geometry_context3)

    HN_filt = HarderNarasimhanFiltration(stable_objects=[lb1,lb2,lb3,lb4,lb5],
                                        phase_vector=[0.2, 0.4, -1.1, 3.2, 5.9],
                                        dimension_vector=[1, 1, 1, 4, 1],
                                        shift_vector=[0, 0, 2, 4, 5])
    print("\n")
    print(HN_filt)

    print("\n\nMin:")
    print(min(HN_filt))

    print("\n\nMax:")
    print(max(HN_filt))

    assert isinstance(min(HN_filt), HNFactor)
    assert isinstance(max(HN_filt), HNFactor)


    with pytest.raises(TypeError):
        HN_filt = HarderNarasimhanFiltration(stable_objects=lb1, phase_vector=[0.2])


def test_GradedCoprod_to_HN_filt():
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

    graded_coproduct = GradedCoproductObject(object_vector=[lb1, lb2], shift_vector=[-1,0], dimension_vector=[10, 1])

    phases = [2.1, 3.2]

    print("\n\n")
    print(graded_coproduct)

    hn_filt = HarderNarasimhanFiltration( stable_objects=[], phase_vector=[] )

    print("\n\n")
    print(hn_filt)

    hn_filt += HarderNarasimhanFiltration(stable_objects=[lb1], shift_vector=[-1], phase_vector=[2.1], dimension_vector=[10])
    
    print("\n\n")
    print(hn_filt)

    hn_filt += HarderNarasimhanFiltration(stable_objects=[lb2], shift_vector=[0], phase_vector=[3.2], dimension_vector=[1])
    print("\n\n")
    print(hn_filt)

    hn_filt *= 2
    print("\n\n")
    print(hn_filt)
    
    