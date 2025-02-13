import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import LocalP2_stab
import numpy as np

def test_chain_complex():
    # Test the chain complex
    linebundle1 = LocalP2_stab.LineBundle(-3)
    linebundle2 = LocalP2_stab.LineBundle(-2)
    linebundle3 = LocalP2_stab.LineBundle(-1)
    linebundle4 = LocalP2_stab.LineBundle(0)

    chaincomplex = LocalP2_stab.ChainComplex([linebundle1, linebundle2, linebundle3, linebundle4], -5)

    assert str(chaincomplex) == '-5         -4         -3         -2\nO(-3) ---> O(-2) ---> O(-1) ---> O(0)'


def test_chern_class():
    # Test the Chern class
    linebundle1 = LocalP2_stab.LineBundle(3)
    chern = linebundle1.chernClass()

    assert chern.ch2 == 4.5
    assert chern.ch1 == 3
    assert chern.ch0 == 1


    linebundle2 = LocalP2_stab.LineBundle(-2)
    chern = linebundle2.chernClass()

    assert chern.ch2 == 2.0
    assert chern.ch1 == -2
    assert chern.ch0 == 1

    linebundle3 = LocalP2_stab.LineBundle(0)
    chern = linebundle3.chernClass()

    assert chern.ch2 == 0.0
    assert chern.ch1 == 0
    assert chern.ch0 == 1

    chaincomplex = LocalP2_stab.ChainComplex([linebundle1, linebundle2, linebundle3], -5)

    chern = chaincomplex.chernClass()

    assert chern.ch2 == -2.5
    assert chern.ch1 == -5.0
    assert chern.ch0 == -1.0


def test_spherical_twist():

    