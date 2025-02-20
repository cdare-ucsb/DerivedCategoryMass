import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import ChernCharacter, CoherentSheaf, LineBundle, ChainComplex



def test_chern_character():
    # Test the Chern class
    linebundle1 = LineBundle(3)
    chern = linebundle1.chernCharacter()

    assert chern.ch2 == 4.5
    assert chern.ch1 == 3
    assert chern.ch0 == 1


    linebundle2 = LineBundle(-2)
    chern = linebundle2.chernCharacter()

    assert chern.ch2 == 2.0
    assert chern.ch1 == -2
    assert chern.ch0 == 1

    linebundle3 = LineBundle(0)
    chern = linebundle3.chernCharacter()

    assert chern.ch2 == 0.0
    assert chern.ch1 == 0
    assert chern.ch0 == 1

    chaincomplex = ChainComplex([linebundle1, linebundle2, linebundle3], [-5,-4,-3])

    chern = chaincomplex.chernCharacter()

    assert chern.ch2 == -2.5
    assert chern.ch1 == -5.0
    assert chern.ch0 == -1.0
