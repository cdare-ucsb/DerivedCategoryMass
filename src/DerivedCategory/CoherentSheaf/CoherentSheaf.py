from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject
from sympy import Expr, Symbol, PolynomialError
from typing import List


###############################################################################
#                                                                             #
#                            Single-Degree Objects                            #
# ----------------------------------------------------------------------------#
#  These objects are used to represent coherent sheaves, vector bundles, and  #
#  line bundles. Line bundles (i.e. vector bundles of rank 1) will be our     #
#  building blocks for the sake of this project, since every coherent sheaf   #
#  in projective space has a resolution by line bundles. Moreover, all        #         
#  spherical twists that we will consider can be represented as a composition #
#  of twists around (the pushforward of) line bundles.                        #    
#                                                                             #
###############################################################################


from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']



class CoherentSheaf(DerivedCategoryObject):
    r"""!
    Generic class for coherent sheaves on a projective variety. This class is
    intended to be subclassed by more specific classes like LineBundle, which
    represent line bundles on projective varieties. By itself, it does not
    encode much more data than the Chern Character of the sheaf, since 
    numerical stability conditions technicially only depend on the Chern
    Character. 

    The parent class to this is DerivedCategoryObject, which is a more general
    class that theoretically does not even require a ChernCharacter --- just
    a string label.
    """

    _instances = {}

    def __new__(cls, chern_character : Expr, catagory : str, allowed_basis : List[Symbol] =None):

        key = (tuple(chern_character), catagory, tuple(allowed_basis) if allowed_basis else None)
        if key not in cls._instances:
            instance = super().__new__(cls)
            cls._instances[key] = instance
    
        return cls._instances[key]

    
    def __init__(self, chern_character : Expr, catagory : str, allowed_basis : List[Symbol] =None):
        r"""!
        Initializes an instance of CoherentSheaf with the specified Chern Character
        and catagory.

        \param chern_character: ChernCharacter
            The Chern Character of the coherent sheaf
        \param catagory: str
            The catagory of the coherent sheaf. Currently implemented catagories
            are 'P1', 'P2', and 'K3'

        \raises ValueError: If the catagory is not implemented, or if the Chern Character is not of the correct length
        \raises TypeError: If the Chern Character is not an instance of Chern
        """

        if hasattr(self, '_initialized') and self._initialized:
            return

        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {catagory} is not implemented.")
        
        
        if not isinstance(chern_character, ChernCharacter):
            raise TypeError("Chern Character must be a ChernCharacter object")


        self.catagory = catagory ## The catagory of the coherent sheaf (e.g. 'P1', 'P2', 'K3')

        self.chern_character = chern_character ## The Chern Character of the coherent sheaf, passed as a ChernCharacter object.

        self.allowed_basis = allowed_basis ## The allowed basis of divisors for the coherent sheaf

        self._initialized = True ## Mark the instance as initialized

        

    def chernCharacter(self):
        r"""!
        Simply accessor method that is intended to be overridden by subclasses, but is 
        also implemented for more general parent / container classes like ChainComplex
        and DerivedCategoryObject. Using the same name allows more modularity in the
        code, and allows for more general functions to be written that can be applied
        to a variety of objects.

        \return ChernCharacter The Chern Character of the coherent sheaf
        """

        return self.chern_character

        
    def shift(self, n):
        r"""!
        Override of the DerivedChatagoryObject shift method. This method shifts the coherents sheaf,
        considered as a complex concentrated in degree 0, by n units. The implementation of this 
        method is crucial to allow including the sheaf in a distinguished triangle, since any triangle
        can be rotated right or left. Since a Coherent sheaf does not keep track of its cohomological 
        information, the method must return a CoherentSheafCoproduct concentrated in a single degree.

        \param int n
            The number of units to shift the coherent sheaf by

        \return ChainComplex A ChainComplex concentrated in a single degree, shifted by n units
        """

        from src.DerivedCategory.GradedCoproductObject import LineBundleCoproduct # include in the method to avoid circular import
        return LineBundleCoproduct( sheaf_vector=[self], shift_vector=[n], dimension_vector=[1])
        

    def __str__(self):
        r"""!
        String representation of the coherent sheaf. This is intended to be overridden by subclasses
        to provide a more informative string representation. Currently this method only returns the
        following string:

        'CoherentSheaf with Chern Character <ch0, ch1, ch2>'

        assuming that the Chern character has 3 entries.

        \return str A string representation of the coherent sheaf
        """

        return f'CoherentSheaf with Chern Character {self.chern_character}'
    

    def get_HN_factors(self, *args):
        return super().get_HN_factors(*args)
    
    def is_semistable(self, *args):
        return super().is_semistable(*args)





class LineBundle(CoherentSheaf):
    """!
    Main class for line bundles on a projective variety. Line bundles are specifically locally free
    sheaves (i.e. vector bundles) of rank 1. In the cases of Local P1 and Local P2, the line bundles
    will serve as the building blocks of the derived category, since every coherent sheaf admits a 
    resolution by line bundles (coming from the `canonical' exceptional collection). However, on K3 
    surfaces, it is not generally true that even vector bundles can be decomposed into sums of line
    bundels.

    Since LineBundles are specifically rank 1, their second chern class c_2(L) always vanishes. Since
    the second Chern Character is (c_1^2(E) - c_2(E))/2, the second Chern Character of a line bundle
    is simply c_1^2(E)/2; therefore, this class will conveniently only store the degree of the line
    bundle, which is the first Chern Character.
    """



    def __init__(self, divisor : Expr, catagory : str, allowed_basis : List[Symbol] = None):
        r"""!
        Initializes an instance of LineBundle with the specified degree and catagory. The Chern Character
        of the line bundle is automatically computed based on the degree and catagory.

        \param int degree The degree of the line bundle
        \param str catagory The catagory of the line bundle. Currently implemented catagories are 'P1', 'P2', and 'K3'

        \raises ValueError If the degree is not an integer
        \raises NotImplementedError If the catagory is not implemented


        \var catagory str The catagory of the line bundle
        \var chern_character ChernCharacter The Chern Character of the line bundle
        """
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise NotImplementedError(f"Catagory {catagory} is not implemented.")
        if not isinstance(divisor, Expr):
            raise ValueError(f"Divisor must be a SymPy Expr: currently passed {type(divisor)}")
        
        if allowed_basis is None:
            allowed_basis = list(divisor.free_symbols)
    
        else:
            used_vars = divisor.free_symbols
            allowed_vars = set(allowed_basis)
            if not used_vars <= allowed_vars:
                raise ValueError(f"divisor uses symbols {used_vars - allowed_vars} which are not in allowed basis {allowed_vars}")

        # Check it's a linear polynomial (degree 1 or less)
        try:
            poly = divisor.as_poly(*allowed_basis)
            if poly.total_degree() > 1:
                raise ValueError("divisor must be linear in the basis elements")
        except PolynomialError:
            raise ValueError("divisor is not a valid polynomial in the basis")


        dimension = 0

        if catagory == 'P1':
            dimension = 2
        elif catagory == 'P2' or catagory == 'K3':
            dimension = 3
        else:
            raise NotImplementedError(f"Catagory {catagory} is not implemented.")

        super().__init__(ChernCharacter.exp(divisor, dimension=dimension), catagory, allowed_basis=allowed_basis)

        self.divisor = divisor ## The divisor of the line bundle, which is a polynomial in the allowed basis

        
        


    def is_semistable(self, *_):
        r"""!
        A result of Macrì-Schmidt (Lectures on Bridgeland Stability, 2016) is that whenever a surface
        has Picard rank 1, line bundles are stable everywhere. This will specifically be used for the
        case of Local P2 (in which case the pushforward i_* preserves this fact), and K3 surfaces. For
        local P1, the line bundles are stable everywhere by definition of the tilt.

        \param tuple args
            The parameters of the stability condition. The number of parameters should be equal to the
            number of parameters required by the central charge for the given catagory. For example, a P1
            object requires a single complex number parameter, while a P2 object requires two real number
            parameters. These are not in fact used, but included to match the format of other classes.

        \return bool True since line bundles are stable in our currently implemented examples.
        """

        
        return True


    def __str__(self):
        r"""!
        String representation of the line bundle. Since all of our implemented catagories come 
        from objects of Picard Rank 1, the line bundles are all derived from the structure sheaf.
        In particular, we can represent any line bundle as O(d) for some integer d.

        \return str A string representation of the line bundle, with the format 'O(d)'
        """

        return f'O({self.divisor})'
    
    def __eq__(self, other):
        r"""!
        Equality comparison for line bundles. Two line bundles are considered equal if they have the
        same degree and catagory.

        \param LineBundle other The line bundle to compare to

        \return bool True if the line bundles have the same catagory and degree, False otherwise
        """
        if not isinstance(other, LineBundle):
            return False
        if other.catagory != self.catagory:
            return False
        return self.divisor == other.divisor
    
    