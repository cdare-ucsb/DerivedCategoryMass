import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import LocalP2_stab as LP2
import numpy as np

def test_chain_complex():
    # Test the chain complex
    linebundle1 = LP2.LineBundle(-3)
    linebundle2 = LP2.LineBundle(-2)
    linebundle3 = LP2.LineBundle(-1)
    linebundle4 = LP2.LineBundle(0)

    chaincomplex = LP2.ChainComplex([linebundle1, linebundle2, linebundle3, linebundle4], -5)

    assert str(chaincomplex) == '-5         -4         -3         -2\nO(-3) ---> O(-2) ---> O(-1) ---> O(0)'


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

    chaincomplex = LP2.ChainComplex([linebundle1, linebundle2, linebundle3], -5)

    chern = chaincomplex.chernCharacter()

    assert chern.ch2 == -2.5
    assert chern.ch1 == -5.0
    assert chern.ch0 == -1.0




# def test_spherical_twist():

    


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
    