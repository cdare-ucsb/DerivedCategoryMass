import pytest
import os
import sys

from sympy import symbols, expand


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.CoherentSheaf import LineBundle



def test_CoherentSheaf_init():
    # Test the initialization of the CoherentSheaf class
    
    H, C = symbols("H C")

    lb1 = LineBundle(3*H, catagory='K3', allowed_basis=[H,C])

    assert lb1.divisor == 3*H
    assert lb1.catagory == 'K3'
    assert lb1.allowed_basis == [H, C]

    assert lb1.chern_character[2] == expand(9*H**2 / 2)

    H2 = symbols("H2")
    lb2 = LineBundle(3*H2, catagory='P1')

    assert lb2.divisor == 3*H2
    assert lb2.catagory == 'P1'
    assert lb2.allowed_basis == [H2]