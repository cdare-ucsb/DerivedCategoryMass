import pytest
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src import LineBundle, ChainComplex




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


def test_adding_characters():

    linebundle1 = LineBundle(3)
    linebundle2 = LineBundle(-2)

    chern1 = linebundle1.chernCharacter()
    chern2 = linebundle2.chernCharacter()

    chern = chern1 + chern2

    assert chern.ch2 == 6.5
    assert chern.ch1 == 1
    assert chern.ch0 == 2

    chern = chern1 - chern2

    assert chern.ch2 == 2.5
    assert chern.ch1 == 5
    assert chern.ch0 == 0