from sympy import symbols, expand, S, prod
from itertools import permutations, product
from collections import defaultdict

class IntersectionStructure:

    def __init__(self, basis, dimension, tensor_data):
        
        self.basis = basis
        self.dimension = dimension
        self.tensor = defaultdict(lambda: S(0))

        # Symmetrize input
        for key, val in tensor_data.items():
            for perm in set(permutations(key)):
                self.tensor[perm] += S(val) 


    def _expand_expr(self, expr):
        """
        Return a list of (coefficient, [monomial factors]) for a SymPy expression.
        For example, D**2 * H returns (1, [D, D, H])
        """
        expr = expand(expr)
        terms = []
        for term in expr.as_ordered_terms():
            coeff, monomial = term.as_coeff_Mul()

            # Break monomial into base ** exp factors
            factors = []
            for factor in monomial.as_ordered_factors():
                if factor.is_Pow:
                    base, exp = factor.args
                    if not exp.is_Integer or exp < 0:
                        raise ValueError(f"Unsupported power: {factor}")
                    factors.extend([base] * int(exp))
                else:
                    factors.append(factor)

            terms.append((coeff, factors))
        return terms

    def evaluate(self, *exprs):
        """
        Compute the n-fold intersection number T(expr1, ..., exprn)
        where each expr is a SymPy linear combination of basis elements.
        """

        expr_terms = [self._expand_expr(e) for e in exprs]
        total = S(0)

        for combo in product(*expr_terms):
            coeffs = [c for c, _ in combo]
            mon_lists = [m for _, m in combo]
            flat_mons = [m for sublist in mon_lists for m in sublist]

            if len(flat_mons) != self.dimension:
                continue  # skip invalid combinations (e.g. D^2 * H in surface with dim=2)

            val = self.tensor.get(tuple(flat_mons), 0)
            total += S(val) * prod(coeffs)

        return total
    




class ChernCharacter:
    def __init__(self, expr, basis, dimension):
        self.expr = expand(expr)
        self.basis = basis
        self.dimension = dimension

    def degree_part(self, degree):
        """
        Extracts the part of the expression that lies in total degree = `degree`
        with respect to the basis.
        """
        degree_expr = 0
        for term in self.expr.as_ordered_terms():
            coeff, mon = term.as_coeff_Mul()
            try:
                poly = mon.as_poly(*self.basis)
                total_deg = poly.total_degree()
            except Exception:
                total_deg = 0  # mon = 1 or constant

            if total_deg == degree:
                degree_expr += coeff * mon

        return degree_expr
    
    

    def __repr__(self):
        return f"ChernCharacter({self.expr})"

    def __add__(self, other):
        return ChernCharacter(self.expr + other.expr, self.basis)

    def __mul__(self, other):
        if isinstance(other, ChernCharacter):
            if self.basis != other.basis:
                raise ValueError("Cannot multiply Chern characters with different bases.")
            if self.dimension != other.dimension:
                raise ValueError("Cannot multiply Chern characters with different dimensions.")
            

            new_expr = expand(self.expr * other.expr)

            truncated_expr = 0
            for term in new_expr.as_ordered_terms():
                coeff, mon = term.as_coeff_Mul()
                try:
                    poly = mon.as_poly(*self.basis)
                    total_deg = poly.total_degree()

                except Exception:
                    total_deg = 0

                if total_deg <= self.dimension:
                    truncated_expr += coeff * mon

            return ChernCharacter(expand(truncated_expr), self.basis, self.dimension)
        
        elif isinstance(other, (float,int)):
            return ChernCharacter(expand(self.expr * S(other)), self.basis, self.dimension)
        else:
            raise TypeError("Can only multiply ChernCharacter by another ChernCharacter or a scalar.")


    def __rmul__(self, other):
        return self.__mul__(other)
    

    def total_expression(self):
        return self.expr
    

if __name__ == "__main__":

    # Basis
    D, H = symbols("D H")

    # Coefficients for ch(E)
    a, b, c, d, e, f = symbols("a b c d e f")

    # Chern character expression
    expr = 1 + 3*D + 5*H + 6*D**2 + 7*D*H + 9*H**2
    ch1 = ChernCharacter(expr, basis=[D, H], dimension=2)

    print("Chern character:", ch1)
    print("Degree 0 part:", ch1.degree_part(0))  # a
    print("Degree 2 part:", ch1.degree_part(1))  # b*D + c*H
    print("Degree 4 part:", ch1.degree_part(2))  # d*D**2 + e*D*H + f*H**2

    ch2 = ChernCharacter(2*D + 3*H, basis=[D, H], dimension=2)

    ch3 = ch1 * ch2
    print("\n\n\n")
    print(f"Multiplication of \n\n{ch1}\n{ch2}\n-----------------\n{ch1*ch2}")

    print("\n\n\n")
    print("Extraction of last term of Chern character")
    print(ch3.degree_part(2)) # 6*D**2 + 19*D*H + 15*H**2


    tensor_data = {
        (D, D): 1,
        (D, H): 2,
        (H, H): 3,
    }

    print("\n\n\nEvaluate using the tensor data")
    intersection_structure = IntersectionStructure([D, H], 2, tensor_data)
    print(intersection_structure.evaluate(ch3.degree_part(2))) # 6 + 2(19) + 3(15) = 6 + 38 + 45 = 89




