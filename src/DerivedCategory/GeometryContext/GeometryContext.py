from .DivisorData import DivisorData
from sympy import Symbol, Add, Expr, Mul


from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']


class GeometryContext():

    def __init__(self, catagory : str, divisor_data : DivisorData, polarization : Symbol = None):

        if catagory not in IMPLEMENTED_CATAGORIES:
            raise NotImplementedError(f"Catagory {catagory} is not implemented.")
        
        if not isinstance(divisor_data, DivisorData):
            raise TypeError("Intersection form must be an instance of DivisorData object.")
        
        if polarization is not None and not isinstance(polarization, Symbol):
            raise TypeError("Polarization must be a SymPy symbol.")
        if polarization is not None and not _is_linear_combination(polarization, divisor_data.basis):
            raise ValueError("Polarization must be a linear combination of elements in the divisor basis.")
        
        self.polarization = polarization
    

        ## Verify that divisor_data is valid
        if catagory == "P1":
            # Verify that intersection_form data is valid
            if divisor_data.variety_dimension != 1:
                raise ValueError("P1 variety must have dimension 1.")
            if len(divisor_data.basis) != 1:
                raise ValueError("P1's Picard group is generated by a single element")
            if divisor_data.top_intersection_form != {(divisor_data.basis[0],): 1}:
                raise ValueError("P1 variety's top intersection form must assign the single divisor to 1.")
            if polarization is not None and polarization != divisor_data.basis[0]:
                raise ValueError("P1 variety's polarization must be the single ample generator.")
            
        elif catagory == "P2":
            
            if divisor_data.variety_dimension != 2:
                raise ValueError("P2 variety must have dimension 2.")
            if len(divisor_data.basis) != 1:
                raise ValueError("P2's Picard group is generated by a single element")
            if divisor_data.top_intersection_form != {(divisor_data.basis[0],divisor_data.basis[0]): 1}:
                raise ValueError("P2 variety's top intersection form must assign the squared ample generator to 1.")
            if polarization is not None and polarization != divisor_data.basis[0]:
                raise ValueError("P2 variety's polarization must be the single ample generator.")
            
        elif catagory == 'K3':
            if divisor_data.variety_dimension != 2:
                raise ValueError("K3 variety must have dimension 2.")
            
        else:
            # TODO: Maybe adjust this to allow for custom classes
            raise NotImplementedError(f"Catagory {catagory} is not implemented.")
            
        self.catagory = catagory
        self.divisor_data = divisor_data

    def __hash__(self):
        return hash(self.catagory)
    
    def __eq__(self, other):
        if not isinstance(other, GeometryContext):
            return False
        return self.catagory == other.catagory and self.divisor_data == other.divisor_data and self.polarization == other.polarization
    
    def variety_dimension(self):
        return self.divisor_data.variety_dimension
    



def _is_linear_combination(expr: Expr, basis_symbols: list[Symbol]) -> bool:
    """
    Checks if expr is a linear combination of basis_symbols.
    That is, expr = sum_i (a_i * b_i), where b_i in basis_symbols and a_i ∈ ℚ or ℝ.
    No constants, no nonlinear terms, and no unknown symbols allowed.
    """
    expr = expr.expand()
    basis_set = set(basis_symbols)

    terms = expr.as_ordered_terms()

    for term in terms:
        if term.is_Atom:
            # A bare symbol or number is not allowed unless it's in basis
            if term in basis_set:
                continue
            return False

        coeff, rest = term.as_coeff_Mul()

        if rest.is_Atom:
            # Must be a single basis symbol
            if rest not in basis_set:
                return False
        else:
            # If rest is a product or something more complex: invalid
            # E.g. H*C, H**2, etc.
            symbols_in_term = rest.free_symbols
            if not symbols_in_term <= basis_set:
                return False
            return False

    return True