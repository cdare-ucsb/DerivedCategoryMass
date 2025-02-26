import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.ChernCharacter import ChernCharacterP2
from src.CoherentSheaf import CotangentBundleP2, LineBundle 
from src.DerivedCategoryObject import DerivedCategoryObject, ChainComplex



def test_chain_complex():
    # Test the chain complex
    linebundle1 = LineBundle(-3)
    linebundle2 = LineBundle(-2)
    linebundle3 = LineBundle(-1)
    linebundle4 = LineBundle(0)

    chaincomplex = ChainComplex([linebundle1, linebundle2, linebundle3, linebundle4], [3, 1, 4, 2], [1,2,3,4])


    assert len(str(chaincomplex)) == 129





def test_complex_central_charge():

    linebundle1 = LineBundle(-3)
    linebundle2 = LineBundle(-2)

    s = 0.5
    q = 0.9

    assert linebundle1.central_charge(s,q) == complex(-3.6,-3.5)
    assert linebundle2.central_charge(s,q) == complex(-1.1,-2.5)


    chaincomplex = ChainComplex([linebundle1, linebundle2], [1, 2])
    assert chaincomplex.central_charge(s,q) == complex(2.5,1)

    assert chaincomplex.get_largest_phase(s,q) == 1.368058363928518
    assert chaincomplex.get_smallest_phase(s,q) == 0.2455170585827645

