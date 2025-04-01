import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.DerivedCategory.ChernCharacter import ChernCharacterP2
from src.DerivedCategory.CoherentSheaf.CoherentSheaf import CotangentBundleP2, LineBundle
from src.DerivedCategory.DerivedCategoryObject import ChainComplex


def test_is_cotangent_bundle_sum():

    cot = CotangentBundleP2(0)
    assert cot.isCotangentBundleSum() == True
    assert cot.chernCharacter() == ChernCharacterP2(2, -3, 1.5)


    cot1 = CotangentBundleP2(1)
    assert cot1.isCotangentBundleSum() == True
    assert cot1.chernCharacter() == ChernCharacterP2(2, -1, -0.5)

    # Show that the Euler sequence holds
    O0 = LineBundle(0)
    Oneg1 = LineBundle(-1)
    assert cot.chernCharacter() + O0.chernCharacter() ==  3 * Oneg1.chernCharacter()


    sum_of_cots = 3 * cot1.chernCharacter() 

    # Make sure this correctly identifies sums of chern characters
    assert sum_of_cots.isCotangentBundleSum()

    incorrect_sum_of_cots = 3 * cot1.chernCharacter() + cot.chernCharacter()
    assert not incorrect_sum_of_cots.isCotangentBundleSum()