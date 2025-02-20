import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import ChernCharacter, CotangentBundle, LineBundle, ChainComplex


def test_is_cotangent_bundle_sum():

    cot = CotangentBundle(0)
    assert cot.isCotangentBundleSum() == True
    assert cot.chernCharacter() == ChernCharacter(2, -3, 1.5)


    cot1 = CotangentBundle(1)
    assert cot1.isCotangentBundleSum() == True
    assert cot1.chernCharacter() == ChernCharacter(2, -1, -0.5)

    # Show that the Euler sequence holds
    O0 = LineBundle(0)
    Oneg1 = LineBundle(-1)
    assert cot.chernCharacter() + O0.chernCharacter() ==  3 * Oneg1.chernCharacter()


    sum_of_cots = 3 * cot1.chernCharacter() 

    # Make sure this correctly identifies sums of chern characters
    assert sum_of_cots.isCotangentBundleSum()

    incorrect_sum_of_cots = 3 * cot1.chernCharacter() + cot.chernCharacter()
    assert not incorrect_sum_of_cots.isCotangentBundleSum()