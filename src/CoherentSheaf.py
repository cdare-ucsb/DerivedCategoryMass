from ChernCharacter import ChernCharacter
from DerivedCategoryObject import DerivedCategoryObject
import math
import cmath


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


IMPLEMENTED_CATAGORIES = ['P1', 'P2', 'K3']



class CoherentSheaf(DerivedCategoryObject):
    """
    Generic class for coherent sheaves on a projective variety. This class is
    intended to be subclassed by more specific classes like LineBundle, which
    represent line bundles on projective varieties. By itself, it does not
    encode much more data than the Chern Character of the sheaf, since 
    numerical stability conditions technicially only depend on the Chern
    Character. 

    The parent class to this is DerivedCategoryObject, which is a more general
    class that theoretically does not even require a ChernCharacter --- just
    a string label.

    Attributes:
    ----------
    chern_character : ChernCharacter
        The Chern Character of the coherent sheaf
    catagory : str
        The catagory of the coherent sheaf. Currently implemented catagories
        are 'P1', 'P2', and 'K3'
    
    """
    
    def __init__(self, chern_character, catagory):
        """
        Initializes an instance of CoherentSheaf with the specified Chern Character
        and catagory.

        Parameters:
        ----------
        chern_character : ChernCharacter
            The Chern Character of the coherent sheaf
        catagory : str
            The catagory of the coherent sheaf. Currently implemented catagories
            are 'P1', 'P2', and 'K3'

        Raises:
        -------
        ValueError
            If the catagory is not implemented, or if the Chern Character is not of the correct length
        TypeError
            If the Chern Character is not an instance of Chern
        """
        
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {catagory} is not implemented.")
        
        if catagory == 'P1' and len(chern_character) != 2:
            raise ValueError("P1 objects should have a Chern Character of length 2")
        if catagory == 'P2' and len(chern_character) != 3:
            raise ValueError("P2 objects should have a Chern Character of length 3")
        if catagory == 'K3' and len(chern_character) != 3:
            raise ValueError("K3 objects should have a Chern Character of length 3")
        
        if not isinstance(chern_character, ChernCharacter):
            raise TypeError("Chern Character must be a ChernCharacter object")


        self.catagory = catagory
        self.chern_character = chern_character

        

    def chernCharacter(self):
        """
        Simply accessor method that is intended to be overridden by subclasses, but is 
        also implemented for more general parent / container classes like ChainComplex
        and DerivedCategoryObject. Using the same name allows more modularity in the
        code, and allows for more general functions to be written that can be applied
        to a variety of objects.

        Returns:
        -------
        ChernCharacter
            The Chern Character of the coherent sheaf

        """
        return self.chern_character

    def phase(self, *args):
        """
        Computes the phase of the central charge of the coherent sheaf. The central charge
        is an element of the dual of the numerical Grothendieck group; in other words, a 
        funtction

        Z : K -> C

        where K is the numerical Grothendieck group, and C is the complex numbers. The phase
        of the central charge is the argument of this complex number.

        Parameters:
        ----------
        *args : float or int
            The parameters of the central charge. The number of parameters should be equal
            to the number of parameters required by the central charge for the given catagory.
            For example, a P1 object requires a single complex number parameter, while a P2
            object requires two real number parameters.

        Returns:
        -------
        float
            The phase of the central charge of the coherent sheaf, in units of pi

        """
        return cmath.phase(self.central_charge(*args)) / math.pi

    def central_charge(self, *args):
        """
        Computes the central charge of the coherent sheaf. The central charge is a function
        that takes in the parameters of the stability condition, and returns a complex number.
        The central charge is a function of the Chern Character of the sheaf that is additive
        on exact sequences, with coefficients that depend on the stability condition.

        Parameters:
        ----------
        *args : float or int
            The parameters of the central charge. The number of parameters should be equal
            to the number of parameters required by the central charge for the given catagory.
            For example, a P1 object requires a single complex number parameter, while a P2
            object requires two real number parameters.

        Returns:
        -------
        complex
            The central charge of the coherent sheaf as a complex number

        Raises:
        -------
        ValueError
            If the number of arguments is incorrect
        TypeError
            If the arguments are not of the correct type
        NotImplementedError
            If the catagory is not implemented
        """

        if self.catagory == 'P1':
            # check that args is a single complex number
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 central charges should have a single complex parameter")

            return -1*self.chern_character[1] + args[0]*self.chern_character[0]
        
        elif self.catagory == 'P2':
            # check that args is a pair of real numbers
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 central charges should have two real number parameters")

            return complex(-1*self.chern_character[2] +
                            args[1] * self.chern_character[0],
                              self.chern_character[1] - args[0] * self.chern_character[0])
    
        elif self.catagory == 'K3':
            # check that args is a pair of real numbers, as well as an integer representing the degre
            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")

            alpha = args[0]
            beta = args[1]
            d = args[2]
            
            return complex(2*d*alpha * self.chern_character[1] - self.chern_character[2] - self.chern_character[0] + (beta**2 - alpha**2)*d*self.chern_character[0], 
                           2*d*self.chern_character[1] - 2*d*alpha*beta*self.chern_character[0])
    
        else:
            # Catagory is not currently implemented
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        
    def shift(self, n):
        """
        Override of the DerivedChatagoryObject shift method. This method shifts the coherents sheaf,
        considered as a complex concentrated in degree 0, by n units. The implementation of this 
        method is crucial to allow including the sheaf in a distinguished triangle, since any triangle
        can be rotated right or left. Since a Coherent sheaf does not keep track of its cohomological 
        information, the method must return a ChainComplex concentrated in a single degree.

        Parameters:
        ----------
        n : int
            The number of units to shift the coherent sheaf by

        Returns:
        -------
        ChainComplex
            A ChainComplex concentrated in a single degree, shifted by n units
        
        """
        from ChainComplex import ChainComplex # include in the method to avoid circular import
        return ChainComplex( sheaf_vector=[self], shift_vector=[n], dimension_vector=[1])
        

    def __str__(self):
        """
        String representation of the coherent sheaf. This is intended to be overridden by subclasses
        to provide a more informative string representation. Currently this method only returns the
        following string:

        'CoherentSheaf with Chern Character <ch0, ch1, ch2>'

        assuming that the Chern character has 3 entries.

        Returns:
        -------
        str
            A string representation of the coherent sheaf
        """
        return f'CoherentSheaf with Chern Character {self.chern_character}'

    def __hash__(self):
        """
        Hash function for the coherent sheaf. This is implemented to allow for the coherent sheaf to be
        used as a key in a dictionary. This functionality is primarily implemented in the ChainComplex
        class, where a dictionary of Coherent sheaves is used to account for duplicate sheaves of the same
        type.
        """
        return hash(self.chern_character)




class LineBundle(CoherentSheaf):
    """
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

    Attributes:
    ----------
    degree : int
        The degree of the line bundle
    catagory : str
        The catagory of the line bundle. Currently implemented catagories are 'P1', 'P2', and 'K3'
    chern_character : ChernCharacter
        The Chern Character of the line bundle; it will always be stored as [1, degree] for P1 objects,
        and [1, degree, degree^2/2] for P2 and K3 objects.

    """

    def __init__(self, degree, catagory):
        """
        Initializes an instance of LineBundle with the specified degree and catagory. The Chern Character
        of the line bundle is automatically computed based on the degree and catagory.

        Parameters:
        ----------
        degree : int
            The degree of the line bundle
        catagory : str
            The catagory of the line bundle. Currently implemented catagories are 'P1', 'P2', and 'K3'

        Raises:
        -------
        ValueError
            If the degree is not an integer
        NotImplementedError
            If the catagory is not implemented
        """
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise NotImplementedError(f"Catagory {catagory} is not implemented.")
        if not isinstance(degree, int):
            raise ValueError(f"degree must be an integer: currently passed {type(degree)}")
        
        self.degree = degree
        self.catagory = catagory
        if self.catagory == 'K3' or self.catagory == 'P2':
            # Since K3 surfaces and P2 have a second Chern Character, we need to store ch_2
            self.chern_character = ChernCharacter([1, self.degree, float(self.degree**2)/2])
        elif self.catagory == 'P1':
            # For local P1 the Chern character is only the rank and degree
            self.chern_character = ChernCharacter([1, self.degree])
        


    def is_semistable(self, *args):
        """
        A result of Macrì-Schmidt (Lectures on Bridgeland Stability, 2016) is that whenever a surface
        has Picard rank 1, line bundles are stable everywhere. This will specifically be used for the
        case of Local P2 (in which case the pushforward i_* preserves this fact), and K3 surfaces. For
        local P1, the line bundles are stable everywhere by definition of the tilt.

        Parameters:
        ----------
        *args : float or int
            The parameters of the stability condition. The number of parameters should be equal to the
            number of parameters required by the central charge for the given catagory. For example, a P1
            object requires a single complex number parameter, while a P2 object requires two real number
            parameters. These are not in fact used, but included to match the format of other classes.

        
        Returns:
        -------
        bool
            True since line bundles are stable in our currently implemented examples.
        """
        return True


    def __str__(self):
        """
        String representation of the line bundle. Since all of our implemented catagories come 
        from objects of Picard Rank 1, the line bundles are all derived from the structure sheaf.
        In particular, we can represent any line bundle as O(d) for some integer d.

        Returns:
        -------
        str
            A string representation of the line bundle, with the format 'O(d)'

        """
        return f'O({self.degree})'
    
    def __eq__(self, other):
        """
        Equality comparison for line bundles. Two line bundles are considered equal if they have the
        same degree and catagory.

        Parameters:
        ----------
        other : LineBundle
            The line bundle to compare to

        Returns:
        -------
        bool
            True if the line bundles have the same catagory and degree, False otherwise
        """
        if not isinstance(other, LineBundle):
            return False
        if other.catagory != self.catagory:
            return False
        return self.degree == other.degree
    
    def __hash__(self):
        """
        Hash function for the line bundle. This is implemented to allow for the line bundle to be
        used as a key in a dictionary. This functionality is primarily implemented in the ChainComplex
        class, where a dictionary of Coherent sheaves is used to account for duplicate sheaves of the same
        type.

        Returns:
        -------
        int
            The hash of the line bundle
        """
        return hash(self.degree)

