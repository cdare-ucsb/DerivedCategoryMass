import pytest
import os
import sys

from sympy import symbols, expand

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.ChernCharacter import ChernCharacter




def test_ChernCharacter__init__():

    # Test with valid inputs

    x, y, z = symbols("x y z")
    expr1 = 1 + 2*x + 3*y + 4*x**2 + 5*y**2 + 6*x*y
    print(expr1.as_poly(x, y))

    ch1 = ChernCharacter(expr1, basis=[x,y], dimension=2)

    assert ch1.expr == expr1
    assert ch1.basis == [x, y]
    assert ch1.dimension == 2

    with pytest.raises(ValueError):
        # Test with invalid basis
        ChernCharacter(expr=3+z, basis=[x, y], dimension=3)
    with pytest.raises(TypeError):
        # Test with invalid dimension
        ChernCharacter(expr=x+y, basis=[x, y], dimension="a")
    with pytest.raises(ValueError):
        # Test with invalid dimension
        ChernCharacter(expr=x+y, basis=[x, y], dimension=-3)
    with pytest.raises(TypeError):
        # Test with invalid expression
        ChernCharacter(expr="x+y", basis=[x, y], dimension=0)



def test_ChernCharacter_exp():

    # Test with valid inputs
    x, y, z = symbols("x y z")
    basis = [x, y]
    linear_expression = 3*x + y

    ch_exp = ChernCharacter.exp(linear_expression, basis=basis, dimension=2)
    assert ch_exp[0] == 1
    assert ch_exp[1] == linear_expression
    assert ch_exp[2] == expand((linear_expression)**2 / 2)


    # Test with invalid inputs
    with pytest.raises(ValueError):
        # Test with invalid basis
        ChernCharacter.exp(3+z, basis=[x, y], dimension=3)
    with pytest.raises(TypeError):
        # Test with invalid dimension
        ChernCharacter.exp(x+y, basis=[x, y], dimension="a")
    with pytest.raises(ValueError):
        # Test with invalid dimension
        ChernCharacter.exp(x+y, basis=[x, y], dimension=-3)
    with pytest.raises(TypeError):
        # Test with invalid expression
        ChernCharacter.exp("x+y", basis=[x, y], dimension=0)
    with pytest.raises(ValueError):
        # Test with invalid expression
        ChernCharacter.exp(x**2 + y**2, basis=[x, y], dimension=2)


def test_ChernCharacter_multiplication():


    # Test with multiplying two Chern characters

    x, y = symbols("x y")
    basis = [x, y]
    expr1 = 1 + 2*x + 3*y
    expr2 = 4 + 5*x + 6*y
    ch1 = ChernCharacter(expr1, basis=basis, dimension=2)
    ch2 = ChernCharacter(expr2, basis=basis, dimension=2)
    mult = ch1 * ch2

    assert mult[0] == 4
    assert mult[1] ==  13*x + 18*y
    assert mult[2] ==  10*x**2 + 18*y**2 + 27*x*y

    # Test with multiplying by a scalar

    scalar = 5
    mult_scalar = ch1 * scalar
    assert mult_scalar[0] == 5
    assert mult_scalar[1] == 10*x + 15*y
    assert mult_scalar[2] == 0


    with pytest.raises(ValueError):
        # Test with different bases
        z, w = symbols("z w")
        ch3 = ChernCharacter(expr1, basis=[z,w], dimension=2)

        ch2*ch3
    with pytest.raises(ValueError):
        # Test with different dimensions
        ch3 = ChernCharacter(expr1, basis=[x,y], dimension=3)
        ch1*ch3
    with pytest.raises(TypeError):
        # Test with not a ChernCharacter nor a scalar
        ch1 * "x+y"



def test_ChernCharacter_iter():

    # Test with valid inputs
    x, y = symbols("x y")
    expr = 1 + 2*x + 3*y + 4*x**2 + 5*y**2 + 6*x*y
    ch = ChernCharacter(expr, dimension=2)

    # Test the iterator
    for i, term in enumerate(ch):
        if i == 0:
            assert term == 1
        elif i == 1:
            assert term == 2*x + 3*y
        else:
            assert term == 4*x**2 + 5*y**2 + 6*x*y



    assert sum(1 for _ in ChernCharacter.exp(x + y, basis=[x,y], dimension=6)) == 7
