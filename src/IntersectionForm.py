from sympy import Symbol, expand, S, prod
from itertools import permutations, product
from collections import defaultdict
from collections.abc import Mapping



class IntersectionForm:
    r"""!
    A class which represents the matrix for an intersection form on the cohomology ring. Specifically, this class only
    tracks the values obtained when intersecting the top-dimensional number of divisors. For example, if this class is meant
    to model the intersection form on a smooth projective threefold and the Picard group is generated by [X, Y, Z], then the intersection form should only keep tract of the intersection values of X.X.X or X.Y.Z, etc. 

    The intersection form is a symmetric bilinear form on the cohomology ring, which means that one only needs to provide the information of the intersection numbers on a single permutation of the generators. The intersection form is then symmetrized by summing over all permutations of the generators.

    The intersection form is represented as a dictionary, where the keys are tuples of basis elements and the values are the
    intersection numbers. For example, if the intersection form is given by the dictionary { (X, Y): 1, (Y, Z): 2 }, then the
    intersection number of X and Y is 1, and the intersection number of Y and Z is 2. The intersection form is symmetric, so
    the intersection number of Y and X is also 1. The intersection number of X and Z is not defined, so it is not included in the
    dictionary.
    """

    def __init__(self, basis, dimension, tensor_data):
        r"""!
        Initialize the intersection form with a specified basis consisting of SymPy symbols for divisor classes, the dimension of the ambient variety, and a dictionary / Mapping of intersection numbers. As the intersection form is symmetric, the user only needs to provide the intersection numbers for one permutation of the generators. The intersection form is then symmetrized by summing over all permutations of the generators.

        \param basis list A list of SymPy symbols representing the basis elements of the cohomology ring.
        \param dimension int The dimension of the ambient variety.
        \param tensor_data Mapping A dictionary / Mapping of intersection numbers, where the keys are tuples of basis elements and the values are the intersection numbers. The keys must be tuples of the same length as the dimension and the values must be numeric.

        \throws TypeError If the basis is not a list, if the basis elements are not SymPy symbols, if the dimension is not an integer, if the tensor data is not a dictionary, if the tensor keys are not tuples, or if the tensor values are not numeric.
        \throws ValueError If the dimension is negative, if the tensor keys do not have the same length as the dimension, or if the tensor values are not numeric.
        """

        if not isinstance(basis, list):
            raise TypeError("Basis must be a list")
        if not all(isinstance(b, Symbol) for b in basis):
            raise TypeError("All basis elements must be SymPy symbols")
        if not isinstance(dimension, int):
            raise TypeError("Dimension must be an integer")
        if dimension < 0:
            raise ValueError("Dimension must be non-negative")
        if not isinstance(tensor_data, Mapping):
            raise TypeError("Tensor data must be a dictionary")
        if not all(isinstance(k, tuple) for k in tensor_data.keys()):
            raise TypeError("Tensor keys must be tuples")
        if not all(len(k) == dimension for k in tensor_data.keys()):
            raise ValueError("Tensor keys must have the same length as the dimension")
        if not all(isinstance(v, (int, float, S)) for v in tensor_data.values()):
            raise TypeError("Tensor values must be numeric")

        
        self.basis = basis # The list of valid SymPy symbols that the intersection form is defined for

        self.dimension = dimension # The dimension of the ambient variety

        self.tensor = defaultdict(lambda: S(0)) # The dictionary or mapping object used to store the intersection numbers. The default value is 0, which is the identity element for addition.

        # Symmetrize input
        for key, val in tensor_data.items():
            for perm in set(permutations(key)):
                self.tensor[perm] += S(val) 


    def _expand_expr(self, expr):
        r"""!
        Helper function which takes a SymPy expression representing a polynomial in some formal variables over the basis, and breaks the expression into a list of tuples, where each tuple contains the coefficient and a list of the basis elements in the monomial. For example, if the input expression is 2*H**3 + 3*D**2*H, then the output will be [(2, [H, H, H]), (3, [D, D, H])]. The function also checks that all symbols in the expression are in the basis.

        \param expr SymPy expression A SymPy expression representing a polynomial in some formal variables over the basis. The expression must be a linear combination of the basis elements.

        \throws ValueError If the expression contains symbols not in the basis, or if the expression contains powers of symbols that are not integers or are negative.

        \return list A list of tuples, where each tuple contains the coefficient and a list of the basis elements in the monomial. For example, if the input expression is 2*H**3 + 3*D**2*H, then the output will be [(2, [H, H, H]), (3, [D, D, H])].
        """

        # Verify that the expression is a linear combination of the basis elements
        used_symbols = expr.free_symbols
        invalid = used_symbols - set(self.basis)
        if invalid:
            raise ValueError(f"Expression uses symbols not in the basis: {', '.join(str(s) for s in invalid)}")
    

        expr = expand(expr)
        terms = []
        # Split the expression into its terms (separated by + or -)
        for term in expr.as_ordered_terms():
            # Further split into coefficient and monomial, e.g. 3*H**2*D -> (3, H**2*D)
            coeff, monomial = term.as_coeff_Mul()

            # Break monomial into base ** exp factors
            factors = []
            for factor in monomial.as_ordered_factors():
                if factor.is_Pow:
                    base, exp = factor.args
                    if not exp.is_Integer or exp < 0:
                        raise ValueError(f"Unsupported power: {factor}")
                    # represent H**3 as H*H*H and D**2 as D*D
                    factors.extend([base] * int(exp))
                else:
                    # Separates out individual factors --- e.g. H*D into H and D
                    factors.append(factor)

            terms.append((coeff, factors))
        return terms
    


    def evaluate(self, *exprs):
        r"""!
        Function which computes the intersection number of a given SymPy expression using the stored intersection form / dictionary provided at initialization. The function takes a variable number of SymPy expressions as input, and computes the intersection number obtained by intersecting the monomials in the expressions. The function returns the intersection number as a SymPy expression.

        \param exprs list A variable number of SymPy expressions representing polynomials in some formal variables over the basis. The expressions must be linear combinations of the basis elements.

        \throws ValueError If the degree of some product of terms in the input is not equal to the dimension of the intersection form.

        \return SymPy expression A SymPy expression representing the intersection number obtained by intersecting the monomials in the expressions. The expression is a linear combination of the basis elements.
        """

        expr_terms = [self._expand_expr(e) for e in exprs]
        total = S(0)


        for combo in product(*expr_terms):
            coeffs = [c for c, _ in combo]
            mon_lists = [m for _, m in combo]
            flat_mons = [m for sublist in mon_lists for m in sublist]

            if len(flat_mons) != self.dimension:
                raise ValueError(f"The degree of some product of terms in the input is {len(flat_mons)}, but the intersection form is only defined for degree {self.dimension}.")

            # If the tuple of monomials is not in the dictionary, return 0
            val = self.tensor.get(tuple(flat_mons), 0)
            total += S(val) * prod(coeffs)

        return total
    

