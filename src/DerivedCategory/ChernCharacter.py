from sympy import expand, S, Symbol, factorial, Expr, Add, Mul
from collections import defaultdict
from typing import List





class ChernCharacter():
    r"""!
    A class representing the Chern Character of an object in the derived category of some smooth projective variety.
    The Chern Character is a graded object in the cohomology ring of our variety —— ASSUMING THAT ALL CLASSES ARE ALGEBRAIC,
    we can represent the Chern Character as a polynomial in the basis of divisor classes. This is known to hold for cellular
    varieties (such as Grassmannians), toric varieties, and complete intersections. Ultimately, the representation of a Chern
    Character as a polynomial over some formal variables allows us flexibility in defining things such as intersection numbers,
    exponential characters, and so on.
    """
        

    def __init__(self, expr : Expr, dimension : int, basis : List[Symbol]=None):
        r"""!
        Initialize a Chern Character object from a SymPy expression representing a polynomial in the divisor classes of the variety. The list of viable divisor classes must be passed to the object as well so that we can check that the expression is valid.

        \param expr SymPy expression The expression representing the Chern character. This should be a polynomial in the basis of divisors.

        \param basis list of SymPy symbols The basis of divisors that the Chern character is defined over. This should be a list of SymPy symbols. By default, this is None, which means that the basis will be inferred from the expression.

        \param dimension int The dimension of the underlying variety that the Chern Character is defined over. This is the expected degree of the Chern character's polynomial representation, which allows us to truncate the product of Characters.

        \throws TypeError If the expr is not a polynomial in the given basis of divisors
        \throws TypeError If the basis is not a list of SymPy symbols
        \throws TypeError If the dimension is not an integer
        \throws ValueError If the dimension is negative
        """

        # Input validation
        if basis is not None and not isinstance(basis, list):
            raise TypeError("Basis must be a list of SymPy symbols.")
        if basis is not None and not all(isinstance(b, Symbol) for b in basis):
            raise TypeError("Basis must be a list of SymPy symbols.")
        

        # If the basis is not provided, we will use the free symbols in the expression as the basis
        if basis is not None:
            self.basis = basis ## This is the list of divisors that the Chern character is defined over. This is a list of SymPy symbols.
        else:
            self.basis = list(expr.free_symbols)


        if not isinstance(dimension, (int, float)):
            raise TypeError("Dimension object must be a number that can be cast to an integer.")
        if int(dimension) < 0:
            raise ValueError("Dimension must be a non-negative integer.")
        if not isinstance(expr, Expr):
            raise TypeError("Expression must be a SymPy expression.")
        if not expr.free_symbols.issubset(set(self.basis)) or not expr.is_polynomial(*self.basis):
            raise ValueError("Expression must be a polynomial in the basis of divisors.")


        self.expr = expand(expr) ## This is the SymPy expression for the Chern character as a polynomial in the formal variables corresponding to the basis of divisors.
        
        self.dimension = int(dimension) ## This is the dimension of the underlying variety that the Chern Character is defined over. The dimension is useful since it tells us the expected degree of the Chern polynomial, which allows us to truncate the product of Characters

        self._monomial_coeffs = {} 

        self._degree_components = {}


        #####
        # Set up lookup-dictionary for monomials -> coefficients
        #####
        for term in self.expr.as_ordered_terms():
            coeff, mon = term.as_coeff_Mul()
            
            if mon == 1:
                self._degree_components[0] = self._degree_components.get(0, 0) + coeff
            else:
                try:
                    poly = mon.as_poly(*self.basis)
                    deg = poly.total_degree()
                    self._monomial_coeffs[mon] = coeff
                    self._degree_components[deg] = self._degree_components.get(deg, 0) + coeff * mon
                except Exception:
                    continue




    def __add__(self, other):
        r"""!
        Add two Chern Character objects together. The Chern Character objects must be defined over the same basis of divisors and have the same dimension in order to add them. The resulting Chern Character object is simply the element-wise sum of the two Chern Character objects. This is done by simply adding the two expressions together and simplifying the result.

        \param ChernCharacter other The Chern Character object to add to the current Chern Character object

        \return ChernCharacter The sum of the two Chern Character objects, which is simply the element-wise sum of the graded elements

        \throws TypeError If other is not a ChernCharacter object

        \throws ValueError If the two Chern Character objects are incompatible in the sense that they are defined over different bases or have different dimensions
        """

        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only add ChernCharacter objects together.")
        if not self.basis == other.basis:
            raise ValueError("Chern Characters must have the same basis to add them together.")
        if not self.dimension == other.dimension:
            raise ValueError("Chern Characters must have the same dimension to add them together.")
        
        return ChernCharacter(expand(self.expr + other.expr), dimension=self.dimension, basis=self.basis)
    
    def __sub__(self, other):
        r"""!
        Subtract two Chern Character objects together. The Chern Character objects must be defined over the same basis of divisors and have the same dimension in order to subtract them. The resulting Chern Character object is simply the element-wise difference of the two Chern Character objects. This is done by simply subtracting the two expressions together and simplifying the result.

        \param ChernCharacter other The Chern Character object to subtract from the current Chern Character object

        \return ChernCharacter The difference of the two Chern Character objects, which is simply the element-wise difference of the graded elements

        \throws TypeError If other is not a ChernCharacter object

        \throws ValueError If the two Chern Character objects are incompatible in the sense that they are defined over different bases or have different dimensions
        """

        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only add ChernCharacter objects together.")
        if not self.basis == other.basis:
            raise ValueError("Chern Characters must have the same basis to add them together.")
        if not self.dimension == other.dimension:
            raise ValueError("Chern Characters must have the same dimension to add them together.")
        
        return ChernCharacter(expand(self.expr - other.expr), basis=self.basis, dimension=self.dimension)
    
    def __mul__(self, other):
        r"""!
        Method to define how a ChernCharacter object gets multiplied by another object. As the cohomology ring is a C-algebra, we should naturally be able to multiply ChernCharacter objects by scalars (i.e. floating point numbers or integers). We also want to be able to multiply ChernCharacter objects by other ChernCharacter objects. This is done in a graded fashion, meaning that we will multiply the two ChernCharacter objects like polynomials. For example, if we have two ChernCharacter objects A and B, then the product of A and B will be a new ChernCharacter object whose degree 2 component will consist of all i + j = 2 terms from the two ChernCharacters, and so on.

        \param other (int, float, ChernCharacter) The object that we are multiplying the ChernCharacter by. This can be a scalar or another ChernCharacter object. If it is a ChernCharacter object, we wish to treat the two ChernCharacters like polynomials by multiplying them in a graded fashion.

        \return ChernCharacter The Chern Character object obtained by multiplying the original Chern Character with either a scalar or another Chern Character object. If the other object is a scalar (int or float), this will simply scale all the entries by that number. If the other object is a ChernCharacter, then the resulting ChernCharacter will resemble the product of the two as polynomials --- that is, the degree 2 component will consist of all i + j = 2 terms from the two ChernCharacters, and so on.

        \throws TypeError If scalar is not an integer, floating point number, or ChernCharacter object
        \throws ValueError If the two Chern Character objects are incompatible in the sense that they are defined over different bases or have different dimensions
        """

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

            return ChernCharacter(expr=expand(truncated_expr), basis=self.basis, dimension=self.dimension)
        
        elif isinstance(other, (float,int)):
            return ChernCharacter(expr=expand(self.expr * S(other)), basis=self.basis, dimension=self.dimension)
        else:
            raise TypeError("Can only multiply ChernCharacter by another ChernCharacter or a scalar.")

    

    def __rmul__(self, other):
        r"""!
        This is effectively the same as the __mul__ method, though allows for the scalar to be on the left side of the multiplication.

        \param other (int, float, ChernCharacter) The object that we are multiplying the ChernCharacter by. This can be a scalar or another ChernCharacter object. If it is a ChernCharacter object, we wish to treat the two ChernCharacters like polynomials by multiplying them in a graded fashion.

        \return ChernCharacter The Chern Character object obtained by multiplying the original Chern Character with either a scalar or another Chern Character object. If the other object is a scalar (int or float), this will simply scale all the entries by that number. If the other object is a ChernCharacter, then the resulting ChernCharacter will resemble the product of the two as polynomials --- that is, the degree 2 component will consist of all i + j = 2 terms from the two ChernCharacters, and so on.


        \throws TypeError If scalar is not an integer, floating point number, or ChernCharacter object
        """

        return self.__mul__(other)
    

    
    def __eq__(self, other):
        r"""!
        Check if two Chern Character objects are equal. The Chern Character objects are equal if their graded elements are equal.

        \param ChernCharacter other The Chern Character object to compare to the current Chern Character object

        \return bool True if the Chern Character objects are equal, False otherwise

        """

        if not isinstance(other, ChernCharacter):
            return False
        
        return expand(self.expr) == expand(other.expr) and self.basis == other.basis and self.dimension == other.dimension
    

    def __str__(self):
        r"""!
        String representation of the Chern Character object. 

        \return str A string representation of the Chern Character object
        """

        return str(self.expr)
    
    def __getitem__(self, key):
        r"""!
        Allow for indexing of the Chern Character object, such as ChernCharacter[0] or ChernCharacter[1]. This will return a SymPy expression containing all the monomials of the specified degree in the Chern polynomial corresponding to the character.

        \param 

        \return SymPy expression The SymPy expression fo the monomials corresponding to the specified degree in the Chern polynomial

        \throws TypeError If key is not an integer
        \throws IndexError If key is not in the range of the Chern Character object
        """

        if key in self._monomial_coeffs:
            return self._monomial_coeffs[key]
        
        elif isinstance(key, int):
            return self._degree_components.get(key, 0)
        
        else:
            return 0
        

    
    def __iter__(self):
        r"""!
        Allow for iteration over the Chern Character object by degree.

        \return iter An iterator over the Chern Character object
        """

        return (deg for deg in sorted(self._degree_components.keys()))
    
    def __reversed__(self):
        r"""!
        Return the reversed Chern Character object by degree; this will start from the highest degree component.

        \return reversed The reversed Chern Character object
        """
        return (deg for deg in sorted(self._degree_components.keys(), reverse=True))
    

    def __len__(self):
        r"""!
        Return the length of the Chern Character object; since every ChernCharacter is represented by its ChernPolynomial, the length is simply the highest degree + 1 (accounting for the constant term / degree 0 term)

        \return int The length of the Chern Character object which is just the dimension + 1
        """

        return self.dimension + 1
    
    
    def __contains__(self, expr):
        r"""!
        Check if an item is in the Chern Character object. This is simply a pass-through to the numpy array.

        \param 

        \return bool True if the expr is in the Chern Character object, False otherwise
        """

        if not isinstance(expr, Expr):
            raise TypeError("Only SymPy expressions can be used with 'in'.")

        expr = expand(expr)

        for term in expr.as_ordered_terms():
            coeff, mon = term.as_coeff_Mul()
            if mon == 1:
                if self._degree_components.get(0, 0) != coeff:
                    return False
            else:
                if mon not in self._monomial_coeffs:
                    return False
                if self._monomial_coeffs[mon] != coeff:
                    return False

        return True
        
    
    
    def __setitem__(self, key, value):
        """
        Sets either:
        - the coefficient of a specific monomial (e.g., x, x*y), or
        - the entire symbolic expression of a given degree (e.g., 0, 1, 2).
        """
        from sympy import expand

        # Case 1: monomial key (e.g., x or x*y)
        if isinstance(key, Expr):
            old_coeff = self._monomial_coeffs.get(key, 0)
            self.expr -= old_coeff * key
            self.expr += value * key
            self._monomial_coeffs[key] = value

            # Update the corresponding degree component
            try:
                deg = key.as_poly(*self.basis).total_degree()
                # Recompute degree component afresh
                self._degree_components[deg] = sum(
                    c * m for m, c in self._monomial_coeffs.items()
                    if m.as_poly(*self.basis).total_degree() == deg
                )
            except Exception:
                pass

        # Case 2: scalar (degree 0) component
        elif isinstance(key, int) and key == 0:
            # Remove the old degree-0 scalar from expr
            old_scalar = self._degree_components.get(0, 0)
            self.expr -= old_scalar
            self.expr += value
            self._degree_components[0] = value

        # Case 3: full replacement of a degree component (e.g., ch[2] = 5*x**2 + y**2)
        elif isinstance(key, int):
            # Remove all existing terms of that degree
            to_remove = []
            for mon, coeff in self._monomial_coeffs.items():
                try:
                    if mon.as_poly(*self.basis).total_degree() == key:
                        self.expr -= coeff * mon
                        to_remove.append(mon)
                except Exception:
                    continue
            for mon in to_remove:
                del self._monomial_coeffs[mon]

            old_deg_expr = self._degree_components.get(key, 0)
            self.expr -= old_deg_expr
            self.expr += value
            self._degree_components[key] = value

            # Parse new expression into monomial → coeff
            for term in expand(value).as_ordered_terms():
                coeff, mon = term.as_coeff_Mul()
                if mon != 1:
                    self._monomial_coeffs[mon] = coeff

        else:
            raise KeyError(f"Unsupported key type for __setitem__: {key}")



    
    def __hash__(self) -> int:
        r"""!
        Hash function for the Chern Character object. This is useful for using the Chern Character object as a key in a dictionary or set.

        \return int The hash of the Chern Character object
        """

        return hash((
            tuple(self.expr.as_ordered_terms()),
            tuple(self.basis),
            self.dimension
        ))


    @staticmethod
    def exp(linear_expr, dimension, basis=None):
        r"""!
        This is a factory-method for creating ChernCharacter objects that are the exponential of a linear expression; by exponential of a linear expression, we mean a Chern Character arising from a formal power series of the form
        \[
        \sum_{k=0}^{\infty} \frac{(linear\_expr)^k}{k!}
        \]
        where the linear expression is a polynomial in the basis of divisors. This is useful for computing the Chern character of a coherent sheaf or vector bundle.

        \param linear_expr SymPy expression The linear expression to exponentiate. This should be a polynomial in the basis of divisors.

        \param basis list of SymPy symbols The basis of divisors that the Chern character is defined over. This should be a list of SymPy symbols. If None, the basis will be inferred from the linear expression.

        \param dimension int The dimension of the underlying variety that the Chern Character is defined over. This is the expected degree of the Chern polynomial, which allows us to truncate the product of Characters.

        \return ChernCharacter The Chern Character object obtained by exponentiating the linear expression. This is a Chern Character object that is the exponential of the linear expression.

        \throws TypeError If linear_expr is not a SymPy expression
        \throws TypeError If basis is not a list of SymPy symbols
        \throws TypeError If dimension is not an integer
        \throws ValueError If dimension is negative
        \throws ValueError If the linear expression is not a polynomial in the basis of divisors
        """

        if not isinstance(linear_expr, Expr):
            raise TypeError("Linear expression must be a SymPy expression.")
        if basis is not None and not isinstance(basis, list):
            raise TypeError("Basis must be a list of SymPy symbols.")
        if basis is not None and not all(isinstance(b, Symbol) for b in basis):
            raise TypeError("Basis must be a list of SymPy symbols.")
        
        # If the basis is not provided, we will use the free symbols in the expression as the basis
        if not basis:
            basis = list(linear_expr.free_symbols)


        if not isinstance(dimension, (int, float)):
            raise TypeError("Dimension must be an integer.")
        if dimension < 0:
            raise ValueError("Dimension must be a non-negative integer.")
        if not linear_expr.free_symbols.issubset(set(basis)) or not linear_expr.is_polynomial(*basis):
            raise ValueError("Linear expression must be a polynomial in the basis of divisors.")
        




        ret_expr = S(0)
        for k in range(0, dimension + 1):
            term = expand(linear_expr**k / factorial(k))
            # Truncate any term exceeding total degree = dimension
            truncated_term = 0
            for mon in term.as_ordered_terms():
                coeff, part = mon.as_coeff_Mul()
                try:
                    deg = part.as_poly(*basis).total_degree()
                except Exception:
                    deg = 0
                if deg <= dimension:
                    truncated_term += coeff * part
                else:
                    raise ValueError(f"Input was not a linear expression; found term of degree {deg} in term: {mon}")
            ret_expr += truncated_term

        return ChernCharacter(ret_expr, dimension, basis=basis)





        

