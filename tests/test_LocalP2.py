import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from DeepGanModel.scripts import LocalP2 as LP2
import numpy as np

def test_chain_complex():
    # Test the chain complex
    linebundle1 = LP2.LineBundle(-3)
    linebundle2 = LP2.LineBundle(-2)
    linebundle3 = LP2.LineBundle(-1)
    linebundle4 = LP2.LineBundle(0)

    chaincomplex = LP2.ChainComplex([linebundle1, linebundle2, linebundle3, linebundle4], [3, 1, 4, 2], [1,2,3,4])


    assert len(str(chaincomplex)) == 129




def test_chern_character():
    # Test the Chern class
    linebundle1 = LP2.LineBundle(3)
    chern = linebundle1.chernCharacter()

    assert chern.ch2 == 4.5
    assert chern.ch1 == 3
    assert chern.ch0 == 1


    linebundle2 = LP2.LineBundle(-2)
    chern = linebundle2.chernCharacter()

    assert chern.ch2 == 2.0
    assert chern.ch1 == -2
    assert chern.ch0 == 1

    linebundle3 = LP2.LineBundle(0)
    chern = linebundle3.chernCharacter()

    assert chern.ch2 == 0.0
    assert chern.ch1 == 0
    assert chern.ch0 == 1

    chaincomplex = LP2.ChainComplex([linebundle1, linebundle2, linebundle3], [-5,-4,-3])

    chern = chaincomplex.chernCharacter()

    assert chern.ch2 == -2.5
    assert chern.ch1 == -5.0
    assert chern.ch0 == -1.0




def test_spherical_twist():

    sph_twist_comp = LP2.SphericalTwist(LP2.LineBundle(-3), LP2.LineBundle(-4))
    assert sph_twist_comp.chernCharacter() == LP2.ChernCharacter(4, -13, 21.5)

    


def test_is_cotangent_bundle_sum():

    cot = LP2.CotangentBundle(0)
    assert cot.isCotangentBundleSum() == True
    assert cot.chernCharacter() == LP2.ChernCharacter(2, -3, 1.5)


    cot1 = LP2.CotangentBundle(1)
    assert cot1.isCotangentBundleSum() == True
    assert cot1.chernCharacter() == LP2.ChernCharacter(2, -1, -0.5)

    # Show that the Euler sequence holds
    O0 = LP2.LineBundle(0)
    Oneg1 = LP2.LineBundle(-1)
    assert cot.chernCharacter() + O0.chernCharacter() ==  3 * Oneg1.chernCharacter()


    sum_of_cots = 3 * cot1.chernCharacter() 

    # Make sure this correctly identifies sums of chern characters
    assert sum_of_cots.isCotangentBundleSum()

    incorrect_sum_of_cots = 3 * cot1.chernCharacter() + cot.chernCharacter()
    assert not incorrect_sum_of_cots.isCotangentBundleSum()
    

def test_complex_central_charge():

    linebundle1 = LP2.LineBundle(-3)
    linebundle2 = LP2.LineBundle(-2)

    s = 0.5
    q = 0.9

    assert linebundle1.central_charge(s,q) == complex(-3.6,-3.5)
    assert linebundle2.central_charge(s,q) == complex(-1.1,-2.5)


    chaincomplex = LP2.ChainComplex([linebundle1, linebundle2], [1, 2])
    assert chaincomplex.central_charge(s,q) == complex(2.5,1)

