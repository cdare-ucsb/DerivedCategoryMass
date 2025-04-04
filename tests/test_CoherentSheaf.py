import pytest
import os
import sys

from sympy import symbols, expand, cos


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from src.DerivedCategory.CoherentSheaf import LineBundle



def test_CoherentSheaf_init():
    # Test the initialization of the CoherentSheaf class
    
    H, C = symbols("H C")

    lb1 = LineBundle(3*H, catagory='K3', basis=[H,C])

    assert lb1.divisor == 3*H
    assert lb1.catagory == 'K3'

    assert lb1.chern_character[2] == expand(9*H**2 / 2)

    H2 = symbols("H2")
    lb2 = LineBundle(3*H2, catagory='P1')

    assert lb2.divisor == 3*H2
    assert lb2.catagory == 'P1'

    lb3 = LineBundle(3*H, catagory='K3', basis=[H, C])

    assert lb3 is lb1
    assert lb3 is not lb2

    # Test improper initialization
    with pytest.raises(ValueError):
        lb4 = LineBundle(3*H2, catagory='K3', basis=[H, C])

    with pytest.raises(ValueError):
        lb5 = LineBundle(cos(H), catagory='K3', basis=[H, C])

    with pytest.raises(NotImplementedError):
        lb6 = LineBundle(3*H, catagory='K4', basis=[H, C])




def test_CoherentSheaf_get_HN_factors():


    # Test the get_HN_factors method
    H, C = symbols("H C")

    lb1 = LineBundle(3*H, catagory='K3', basis=[H,C])

    assert lb1.get_HN_factors(1,2,3) == [lb1]

    lb2 = LineBundle(3*H + 4*C, catagory='K3', basis=[H,C])

    assert lb2.get_HN_factors(1,2,3) == [lb2]