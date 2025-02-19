import math
import re


###############################################################################
#                                                                             #
#                           Numerical Characters                              #
# ----------------------------------------------------------------------------#
#  These objects are used to represent Chern classes of coherent sheaves,     #
#  vector bundles, and line bundles. They are used to determine if a given    #
#  object is a sum of line bundles or cotangent bundles.                      #
#                                                                             #
###############################################################################


class ChernCharacter():
    """
    The Chern character is a topological invariant that will heavily be used for computing
    homological information about the coherent sheaves, mainly due to theoretical 
    results such as Hirzebruch-Riemann-Roch. Specifically, for an invertible sheaf L = O(D),
    one has 

    χ(L) = h^0(L) - h^1(L) + h^2(L) - ...
         = ∫_X ch(L) td(X)
         = c_1(L)^2/2 + c_1(L) c_1(T_X)/2 + c_2(L) + ...

    In particular, the dimension of the Hom-spaces will be calculated from the chern characters
    and vice-versa. The Chern character will be represented as a tuple (ch0, ch1, ch2) where
    ch0 is the rank, ch1 is the degree, and ch2 is the second Chern class.

    Attributes:
    ----------
    ch0 : int
        The rank of the coherent sheaf
    ch1 : int
        The degree of the coherent sheaf
    ch2 : float
        The second Chern class of the coherent sheaf
        
    """
    def __init__(self, ch0, ch1, ch2): 
        '''
        Initialize an instance of ChernCharacter with the specified characteristic classes
        Notice that ch0 is the same as the rank and ch1 is the same as the degree. However,
        it is not true that ch2 = c2

        Parameters:
        ----------
        ch0 : int
            The chern character in degree 0 (i.e. the rank for sheaves)
        ch1 : int
            The chern character in degree 1 (i.e. the degree for sheaves)
        ch2 : float
            The chern character in degree 2 (i.e. (c1^2 - 2c2)/2 for sheaves)

        '''
        self.ch0 = ch0
        self.ch1 = ch1
        self.ch2 = ch2

    def __str__(self):
        """
        String representation of the Chern Character class

        Returns:
        -------
        str
            A string representation of the Chern Character
        """
        return f'<{self.ch0}, {self.ch1}, {self.ch2}>'
    
    def __add__(self, other):
        """
        Method to add two Chern Characters together. This is done by adding the corresponding
        components of the Chern Character.

        Parameters:
        ----------
        other : ChernCharacter
            The Chern Character to add to the current Chern Character

        Returns:
        -------
        ChernCharacter
            The sum of the two Chern Characters
        """

        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only add ChernCharacter objects together.")

        return ChernCharacter(self.ch0 + other.ch0, self.ch1 + other.ch1, self.ch2 + other.ch2)
    
    def __sub__(self, other):
        """
        Method to subtract two Chern Characters. This is done by subtracting the corresponding
        components of the Chern Character.

        Parameters:
        ----------
        other : ChernCharacter
            The Chern Character to subtract from the current Chern Character
        
        Returns:
        -------
        ChernCharacter
            The difference of the two Chern Characters
        """

        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only subtract ChernCharacter objects together.")

        return ChernCharacter(self.ch0 - other.ch0, self.ch1 - other.ch1, self.ch2 - other.ch2)
    
    def __mul__(self, scalar):
        """
        Method to multiply a Chern Character by a scalar. This is done by multiplying each
        component of the Chern Character by the scalar.

        Parameters:
        ----------
        scalar : int
            The scalar to multiply the Chern Character by

        Returns:
        -------
        ChernCharacter
            The Chern Character multiplied by the scalar
        """

        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")

        return ChernCharacter(self.ch0 * scalar, self.ch1 * scalar, self.ch2 * scalar)
    
    def __rmul__(self, scalar):
        """
        Method to multiply a Chern Character by a scalar. This is done by multiplying each
        component of the Chern Character by the scalar.

        Parameters:
        ----------
        scalar : int
            The scalar to multiply the Chern Character by

        Returns:
        -------
        ChernCharacter
            The Chern Character multiplied by the scalar
        """

        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")

        return ChernCharacter(self.ch0 * scalar, self.ch1 * scalar, self.ch2 * scalar)
    
    def __eq__(self, other):
        """
        Method to determine if two Chern Characters are equal. This is done by checking if
        the corresponding components of the Chern Character are equal.

        Parameters:
        ----------
        other : ChernCharacter
            The Chern Character to compare to the current Chern Character

        Returns:
        -------
        bool
            True if the two Chern Characters are equal, False otherwise
        """

        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only compare ChernCharacter objects.")

        return self.ch0 == other.ch0 and self.ch1 == other.ch1 and self.ch2 == other.ch2
    
    def isLineBundleSum(self):
        """
        Method to determine if the Chern Character is the same as a certain sum of line bundles,
        up to a shift. The Chern Character of a line bundle is given by (1, deg, deg^2 / 2); since
        the Chern Character is additive, we can determine if a given Chern Character is a sum
        of line bundles by checking if all of the characteristic classes are divisible by
        the rank (loosely speaking, since the second character has a factor of 1/2).

        Returns:
        -------
        bool
            True if the Chern Character is a sum of line bundles, False otherwise
        """
        if self.ch0 == 0:
            # At the very least represents a torsion sheaf. However, any sum of line
            # bundles must have positive rank.
            return False
        
        if self.ch1 % self.ch0 != 0:
            # If the degree is not divisible by the rank, the vector <ch0, ch1, ch2> cannot
            # be a multiple of <1, d, d^2 / 2> for any integer d.
            return False

        # since ch1 = d * ch0, we can determine the possible degree
        possible_deg = int(self.ch1 / self.ch0)
        # check if the second Chern class is of the correct format 
        return self.ch2 == self.ch0 * float(possible_deg**2) / 2
    
    def isCotangentBundleSum(self):
        """
        Method to determine if the Chern Character is the same of a certain sum of cotangent bundles.
        of a FIXED DEGREE. Using the Euler exact sequence

        0 -> Ω(d) -> O(d-1)^3 -> O(d) -> 0

        we can easily verify that the Chern character of the (twisted) cotangent bundle is
        (2, 2d - 3, d^2 - 3d + 3/2). Since the Chern Character is additive, we can determine
        whether a given Chern character represents the sum of cotangent bundles (of a single degree)
        by checking if the characteristic classes are divisible by 2 (loosely speaking, since the
        second character has a factor of 1/2), and if the degree is of the correct form.

        Returns:
        -------
        bool
            True if the Chern Character is a sum of cotangent bundles, False otherwise

        """

        if self.ch0 % 2 != 0:
            # Any sum of cotangent bundles must have even class
            return False

        num_copies = int(self.ch0 / 2) # Cotangent bundle is rank 2, so see how many copies of line bundle there are

        # if this is a cotangent bundles, the degree must be of the form 2d - 3. In particular, 
        # when it is a sum of cotangent bundles, it will be of the form (2d - 3) * num_copies.
        # The following ine simply recovers d.
        original_deg = ( (float(self.ch1) / num_copies) + 3) / 2

        if original_deg != int(original_deg):
            # If the degree is not of the form 2d - 3, then it cannot be a sum of cotangent bundles
            return False

        # Next check that the second Chern class is of the correct form
        # we may in fact use the additivity of the Chern Character on exact sequences
        # as well as the Euler exact sequence above to reduce this problem to making sure
        # the middle term is indeed of the form O(d-1)^3.
        new_bundle = ChernCharacter(self.ch0 + num_copies,
                                self.ch1 + num_copies * original_deg,
                                  self.ch2 + (num_copies * original_deg**2 / 2))
        # Account for possible shift which makes everything negative
        shifted_bundle = ChernCharacter(-self.ch0 - num_copies, 
                                    -self.ch1 - num_copies * original_deg,
                                    -self.ch2 - (num_copies * original_deg**2 / 2))    

        
        
        new_bundle_is_cot =  (new_bundle.isLineBundleSum() and new_bundle.ch0 % 3 == 0)
        shifted_bundle_is_cot = (shifted_bundle.isLineBundleSum() and shifted_bundle.ch0 % 3 == 0)
        return new_bundle_is_cot or shifted_bundle_is_cot
    






###############################################################################
#                                                                             #
#                            Single-Degree Objects                            #
# ----------------------------------------------------------------------------#
#  These objects are used to represent coherent sheaves, vector bundles, and  #
#  line bundles. Line bundles (i.e. vector bundles of rank 1) will be our     #
#  building blocks for the sake of this project, since every coherent sheaf   #
#  on P^2 has a resolution by line bundles. Moreover, all spherical twists    #         
#  that we will consider can be represented as a composition of twists around #
#  (the pushforward of) line bundles.                                         #    
#                                                                             #
###############################################################################

class CoherentSheaf():
    '''
    This class represents a general coherent sheaf on a projective plane. 
    Coherent sheaves are generalizations of vector bundles since they can have torsion.
    Since we will assume that the base field is the complex numbers, we can assume that
    all coherent sheaves are simple represented by C[x,y,z]-modules that are localized 
    in the sense that they agree on the overlap of distinguished affine open sets (e.g.
    U_x = {x != 0}, U_y = {y != 0}, U_z = {z != 0}).

    While coherent sheaves are not uniquely determined by their characteristic classes, 
    we will use Chern classes as the minimal amount of information to represent a coherent sheaf.


    Attributes:
    ----------
    rank : int
        The rank of the coherent sheaf
    deg : int
        The degree of the coherent sheaf
    c2 : float
        The second Chern class of the coherent sheaf
    '''

    def __init__(self, rank, c1, c2 = 0):
        """
        Initialize an instance of CoherentSheaf with the specified chern
        classes.

        Parameters:
        ----------
        rank : int
            The rank of the coherent sheaf
        deg : int
            The degree of the coherent sheaf
        c2 : float
            The second Chern class of the coherent sheaf
        """
        self.c0 = int(rank)
        self.c1 = int(c1)
        self.c2 = c2
    
    def chernCharacter(self):
        """
        Method which returns a Chern Characteristic class of the coherent sheaf
        The chern character is defined on line bundles as follows:

        ch(L) = exp(c1(L)) = 1 + c1(L) + c1^2(L) / 2 + ...

        However, for vector bundles the chern character will generally be written as

        ch(E) = rank(E) + c1(E) + (c_1^2(E) - c2(E)) / 2 + (c_1^3(E) - 3c1(E)c2(E) + c3(E)) / 6 + ...

        The numerics of the Chern Character are handled by a separate class, so that this method
        only returns a specific chern character.

        Returns:
        -------
        ChernCharacter
            The Chern character of the coherent sheaf
        
        """
        ch2 = float( self.c1**2 - 2 * self.c2 ) / 2 # The second Chern class of the sheaf
        return ChernCharacter(self.c0, self.c1, ch2)

    def __str__(self):
        """
        String representation of the coherent sheaf class

        Returns:
        -------
        str
            A string representation of the coherent sheaf
        """
        return f'CoherentSheaf of rank {self.c0}, degree {self.c1}, and c2 {self.c2}'
    


class VectorBundle(CoherentSheaf):
    """
    A vector bundle is a special type of coherent sheaf; in particular, it is a locally free
    sheaf of finite rank. These will ultimately be better generalizations of line bundles, and
    any direct sum of line bundles is an example of a vector bundle. Other notable vector bundles
    that we will use are the tangent and cotangent bundles, which are (somewhat) canonical rank
    2 vector bundles. There is not much additional information to add to the class, since the
    rank, degree, and second Chern class are sufficient to determine the Chern Character. Moreover,
    any additional information would be redundant since we are not keeping track of torsion.

    Attributes:
    ----------
    rank : int
        The rank of the vector bundle
    deg : int
        The degree of the vector bundle
    c2 : float
        The second Chern class of the vector bundle
    """

    def __init__(self, rank, deg, c2):
        """
        Initialize an instance of VectorBundle with the specified chern
        classes.

        Parameters:
        ----------
        rank : int
            The rank of the vector bundle
        deg : int
            The degree of the vector bundle
        c2 : float
            The second Chern class of the vector bundle
        """
        self.rank = int(rank)
        self.deg = int(deg)
        self.c2 = float(c2)

    def __str__(self):
        """
        String representation of the vector bundle class

        Returns:
        -------
        str
            A string representation of the vector bundle
        """
        return f'<{self.rank},{self.deg},{self.c2}>'
    
    def chernCharacter(self):
        """
        Method which returns a Chern Characteristic class of the vector bundle
        The chern character is defined on line bundles as follows:

        ch(L) = exp(c1(L)) = 1 + c1(L) + c1^2(L) / 2 + ...

        However, for vector bundles the chern character will generally be written as

        ch(E) = rank(E) + c1(E) + (c_1^2(E) - c2(E)) / 2 + (c_1^3(E) - 3c1(E)c2(E) + c3(E)) / 6 + ...

        The numerics of the Chern Character are handled by a separate class, so that this method
        only returns a specific chern character.

        Returns:
        -------
        ChernCharacter
            The Chern character of the vector bundle
        
        """
        return super().chernCharacter()
    
    def isLineBundleSum(self):
        """
        Helper function to determine if the vector bundle is a sum of line bundles. The function 
        simpily calls the isLineBundleSum() method of the Chern Character class.

        Returns:
        -------
        bool
            True if the vector bundle is a sum of line bundles, False otherwise
        """
        return self.chernCharacter().isLineBundleSum()
    

    def isCotangentBundleSum(self):
        """
        Helper function to determine if the vector bundle is a sum of cotangent bundles. The function
        simply calls the isCotangentBundleSum() method of the Chern Character class.

        Returns:
        -------
        bool
            True if the vector bundle is a sum of cotangent bundles, False otherwise
        """
        return self.chernCharacter().isCotangentBundleSum()



class LineBundle(VectorBundle):
    """
    A line bundle is a special type of vector bundle; in particular, it is a rank 1 vector bundle.
    These will ultimately be our building blocks for the sake of this project, since every coherent
    sheaf on P^2 has a resolution by line bundles. Moreover, all spherical twists that we will consider
    can be represented as a composition of twists around (the pushforward of) line bundles.

    As the Picard group of P^2 is isomorphic to Z, we can represent any line bundle by its degree.

    Attributes:
    ----------
    c0 : int
        The rank of the line bundle will always be 1
    c1 : int
        The degree of the line bundle
    c2 : float
        The second Chern class of the line bundle, which will always be of the form d^2 / 2
    """
    
    def __init__(self, deg):
        """
        Initialize an instance of LineBundle with the specified degree.

        Parameters:
        ----------
        deg : int
            The degree of the line bundle
        """
        self.c0 = 1
        self.c1 = int(deg)
        self.c2 = 0

    def __str__(self):
        """
        String representation of the line bundle class; the standard notation of a line bundle
        of degree d is O(d), since every line bundle is a twist of the structure sheaf.

        Returns:
        -------
        str
            A string representation of the line bundle
        """
        return f'O({self.c1})'
    


    
class CotangentBundle(VectorBundle):
    """
    An implementation of the canonical bundle of the complex projective plane. The canonical bundle
    is a rank 2 vector bundle, and is the dual of the tangent bundle; it is generated by the
    differentials of the homogeneous coordinates x, y, and z. Locally, these may be represented as
    d(y/z), d(x/z), and d(x/y) (respectively), which simplify to 

                    d(y/z) = dy/z - y dz/z^2    
                    d(x/z) = dx/z - x dz/z^2

    This gives an embedding of sheaves Ω --> O(-1)^3 whose cokernel is the structure sheaf. 

    Attributes:
    ----------
    c0 : int
        The rank of the cotangent bundle, which is 2
    c1 : int
        The first chern class of the cotangent bundle, which is of the form 2d-3
    c2 : float
        The second chern class of the cotangent bundle, which is of the form d^2 -3d + 3
    """
    def __init__(self, deg=0):
        """
        Initialize an instance of CotangentBundle with the specified degree.
        By default, the degree is 0 so that the class fits into the standard Euler Sequence

        Parameters:
        ----------
        deg : int
            The degree of the cotangent bundle
        """
        self.c0 = 2
        self.c1 = 2*int(deg) - 3
        self.c2 = deg**2 - 3*deg + 3

    def __str__(self):
        """
        String representation of the cotangent bundle class; the standard notation of a cotangent bundle
        of degree d is Ω(d).

        Returns:
        -------
        str
            A string representation of the cotangent bundle
        """
        return f'\u03a9({int((self.c1 + 3) / 2)})'
    
    def chernCharacter(self):
        """
        Method which returns a Chern Characteristic class of the cotangent bundle. The chern
        Character of the cotangent bundle is always of the form (2, 2d - 3, d^2 - 3d + 3/2).
        """
        ch2 = float(self.c1**2 - 2 * self.c2) / 2
        return ChernCharacter(self.c0, self.c1, ch2)






###############################################################################
#                                                                             #
#                             Multi-Degree Objects                            #
# ----------------------------------------------------------------------------#
#  These objects are used to represent chain complexes over a the category of #              
#  coherent sheaves; these objects are stored primarily as an array of degree #
#  , ranks, and homological shifts. Despite the fact that ChainComplexes are  #
#  partially implemented, LineBundleComplex objects are much more useful for  #
#  several methods in this project (note in particular that every coherent    #
#  sheaf on P^2 has a resolution by line bundles). Consequently, it is        #
#  important to make sure every function in ChainComplex is overridden in     #
#  LineBundleComplex.                                                         #
#                                                                             #
#  The Distinguished triangle further generalizes the notion of ordered       #                                                     #
#  exact sequences in the Derived category. As our main object of interest,   #
#  SphericalTwistComposition, exists in the Derived category, it will be      #
#  crucial to understand the passage from complexes of coherent sheaves to    #
#  distinguished triangles.                                                   #
#                                                                             #
###############################################################################



class DerivedCategoryObject():
    """
    General parent class for both ChainComplex and LineBundleComplex objects. This class
    is used to represent objects in the derived category of coherent sheaves on P^2. Since
    objects in the derived category are technically equivalence classes of complexes of
    coherent sheaves (up to quasi-isomorphism), we will use an abstract parent class to 
    capture only string information and numerical information.

    Attributes:
    ----------
    string : str
        A string representation of the derived category
    chern_character : ChernCharacter
        The Chern Character of the derived category object
    """


    def __init__(self, string = "0", chern_character = ChernCharacter(0,0,0)):
        """
        Initialize an instance of DerivedCategoryObject with the specified string and Chern Character.

        Parameters:
        ----------
        string : str
            The string representation of the derived category object
        chern_character : ChernCharacter
            The Chern Character of the derived category object
        """
        self.string = string
        self.chern_character = chern_character
    
    def __str__(self):
        """
        String representation of the derived category object

        Returns:
        -------
        str
            A string representation of the derived category
        """
        if self.string is not None:
            return self.string
        else:
            return "0"


    def chernCharacter(self):
        """
        Method to return the Chern Character of the derived category object

        Returns:
        -------
        ChernCharacter
            The Chern Character of the derived category object
        """
        return self.chern_character


    def shiftComplex(self, shift):
        """
        Method to shift the derived category object by a given homological shift

        Parameters:
        ----------
        shift : int
            The homological shift to apply to the derived category object

        Returns:
        -------
        DerivedCategoryObject
            The derived category object shifted by the homological shift
        """
        new_string =  self._update_string_by_shift(self.string, shift)
        new_chern = self.chern_character * (-1)**shift
        return DerivedCategoryObject(new_string, new_chern)


    def _update_string_by_shift(my_str, n):
        """
        Static helper function to update the possible string representations of the abstract
        DerivedCategoryObject by a homological shift. For example, if an object is called A[3],
        then a shift of 2 will yield A[5] (as a string). If the object is unshifted to start with,
        then shifting the object should return A[n] (as a string). The only object that should not
        be shifted is the zero object, which is represented by the string "0".

        Parameters:
        ----------
        my_str : str
            The string representation of the derived category object
        n : int
            The homological shift to apply to the derived category object

        Returns:
        -------
        str
            The string representation of the derived category object shifted by the homological shift

        Raises:
        -------
        TypeError
            If my_str is not a string
            If n is not an integer

        """

        if not isinstance(my_str, str):
            raise TypeError("my_str must be a string.")
        if not isinstance(n, int):
            raise TypeError("n must be an integer.")

        # If the object is zero, it should remain zero
        if my_str == "0":
            return "0"
        
        # This regex checks for a pattern at the end of the string that looks like "[number]"
        pattern = r'\[(\d+)\]$'
        match = re.search(pattern, my_str)
        
        if match:
            # Extract the current number k, add n, and format the new number.
            k = int(match.group(1))
            new_k = k + n
            # Replace the [k] at the end with [new_k]
            updated_str = re.sub(pattern, f'[{new_k}]', my_str)
            return updated_str
        else:
            # If there's no [number] at the end, append [n]
            return my_str + f'[{n}]'







class ChainComplex(DerivedCategoryObject):
    """
    For any abelian category, a chain complex is a sequence of objects and morphisms between them
    such that the composition of any two consecutive morphisms is the zero morphism. In the derived
    category of coherent sheaves on P^2, we can represent a chain complex as a sequence of coherent
    sheaves with a shift. For instance, a general complex will be of the form

              i=-n       i=-n+1      i=-n+2    ...
    0 ------> E1 -------> E2 --------> E3 ---> ...

    (A priori, there is no reason the complexes cant also descend infinitely in the other direction). 
    For the purposes of this project, only finite complexes will be considered. Such a complex can be
    stored in a similar way to a DenseVector object --- namely, since the majority of entries in the
    complex will be zero, we can store the complex as a list of coherent sheaves and a shift vector.

    Attributes:
    ----------
    sheaf_vector : list
        A list of coherent sheaves in the complex
    shift_vector : list
        A list of homological shifts in the complex
    dimension_vector : list
        A list of the number of direct sums of each coherent sheaf in the complex
    
    """

    def __init__(self, sheaf_vector, shift_vector, dimension_vector = None):
        """
        Initialize an instance of ChainComplex with the specified sheaf vector and shift vector.
        The initialization performs several checks to ensure that the input is valid.

        Parameters:
        ----------
        sheaf_vector : list
            A list of coherent sheaves in the complex
        shift_vector : list
            A list of homological shifts in the complex
        dimension_vector : list
            A list of the number of direct sums of each coherent sheaf in the complex. By default,
            it will simply be set to be [1, 1, ..., 1] so that there is only one copy of each sheaf.

        Raises:
        -------
        ValueError
            If the sheaf vector is empty
            If the sheaf vector and shift vector have different lengths
        TypeError
            If any element of the sheaf vector is not a CoherentSheaf object
            If any element of the shift vector is not an integer
        """

        if not sheaf_vector:
            raise ValueError("sheaf_vector cannot be empty.")

        if not all(isinstance(obj, CoherentSheaf) for obj in sheaf_vector):
            raise TypeError("All elements of complex_vector must be instances of CoherentSheaf.")

        if not all(isinstance(shift, int) for shift in shift_vector):
            raise TypeError("All elements of shift_vector must be integers.")

        if len(sheaf_vector) != len(shift_vector):
            raise ValueError("sheaf_vector and shift_vector must have the same length.") 

        if dimension_vector is None:
            dimension_vector = [1] * len(sheaf_vector)
        elif not all(isinstance(dim, int) for dim in dimension_vector):
            raise TypeError("All elements of dimension_vector must be integers.")
        elif len(sheaf_vector) != len(dimension_vector):
            raise ValueError("sheaf_vector and dimension_vector must have the same length.")
        # Dimension cannot be non-negative
        elif not all(dim >= 0 for dim in dimension_vector):
            raise ValueError("All elements of dimension_vector must be non-negative integers.") 
        
        self.sheaf_vector = sheaf_vector
        self.dimension_vector = dimension_vector
        self.shift_vector = shift_vector

        # If an element of the complex has dimension 0, we can get rid of it using helper method
        self._remove_zeros_from_dimension_vector()

        

    def __str__(self):
        """
        String representation of the chain complex. The complex is represented in cohomological order 
        (which technically would be descending order of the shifts, since IR[-2] means the complex with
        a copy of IR in index 2). The individual coherent sheaves in the complex are represented by their
        own respective print functinos --- this will generally be cumbersome for arbitrary vector bundles,
        but more clean for named instances like O(1) or Ω(3). 

        Example
        -------
        For a complex with O(3)

        """
         # Zip the three lists together and sort by descending shift
        bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
        bundles_sorted = sorted(bundles, key=lambda x: x[2], reverse=True)
        
        # Define the arrow string used to join bottom elements
        arrow = " --------> "
        # We'll use the same spacing (without arrow characters) to join top columns
        join_space = " " * len(arrow)
        
        top_columns = []
        bottom_columns = []
        for sheaf, dim, shift in bundles_sorted:
            # Create the bottom string: e.g., "O(5)[3]"
            bottom_elem = f"{sheaf}[{shift}]"
            col_width = len(bottom_elem)
            # Create the top string: e.g., "⊕7" and right-align it to match the bottom column width
            top_elem = f"⊕{dim}"
            top_columns.append(top_elem.rjust(col_width))
            bottom_columns.append(bottom_elem)
        
        top_line = join_space.join(top_columns)
        bottom_line = arrow.join(bottom_columns)
        return top_line + "\n" + bottom_line

        
    

    def chernCharacter(self):
        """
        Helper function to compute the Chern Character of the chain complex. The Chern Character of
        a chain complex is the alternating sum of the Chern Characters of the individual sheaves in
        the complex. Since the Chern character is additive, we may multiply the Chern Characters by
        the dimension of the sheaf to represent direct sums of sheaves.
        """
        cherns = [sheaf.chernClass() for sheaf in self.sheaf_vector]

        ch0 = 0
        ch1 = 0
        ch2 = 0

        for i in range(len(cherns)):
            ch0 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch0
            ch1 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch1
            ch2 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch2

        return ChernCharacter(ch0, ch1, ch2)


    
    def shiftComplex(self, shift):
        """
        Method to shift the chain complex by a given homological shift

        Parameters:
        ----------
        shift : int
            The homological shift to apply to the chain complex

        Returns:
        -------
        ChainComplex
            The chain complex shifted by the homological shift
        """
        new_shift_vector = [shift + s for s in self.shift_vector]

        return ChainComplex(self.sheaf_vector, new_shift_vector, self.dimension_vector)

    def isShiftOfSheaf(self):
        """
        Simple helper function which checks if the complex is a shift of a single sheaf

        Returns:
        -------
        bool
            True if the complex is a shift of a single sheaf, False otherwise
        """
        return len(self.complex) == 1


    def _remove_zeros_from_dimension_vector(self):
        """
        Helper function which iterates through the dimension vector, and if a certain Coherent sheaf
        is only included 0 times, we may effectively erase it.
        """
        for i in range(len(self.dimension_vector)):
            if self.dimension_vector[i] == 0:
                del self.sheaf_vector[i]
                del self.shift_vector[i]
                del self.dimension_vector[i]
                



    


class DistinguishedTriangle():
    """
    A distinguished triangle is a sequence of three objects in the derived category of coherent sheaves
    (in our case) with morphisms of complexes between them such that the third complex is always 
    quasi-isomorphic to the mapping cone of the morphism between the first two complexes. In the derived
    category, distinguished triangles generalize the notion of exact sequences since they also give rise
    to long exact sequences on the level of cohomology.

    In our implementation, we will use the classes DerivedCategoryObject and ChainComplex to represent
    the objects in the distinguished triangle. 

    Attributes:
    ----------
    object1 : DerivedCategoryObject
        The first object in the distinguished triangle
    object2 : DerivedCategoryObject
        The second object in the distinguished triangle
    object3 : DerivedCategoryObject
        The third object in the distinguished triangle

    """

    def __init__(self, derived_object1, derived_object2, derived_object3):
        """
        Initialize an instance of DistinguishedTriangle with the specified objects, and
        then update the Chern Characters of the objects if they are not already set.

        Parameters:
        ----------
        derived_object1 : DerivedCategoryObject
            The first object in the distinguished triangle
        derived_object2 : DerivedCategoryObject
            The second object in the distinguished triangle
        derived_object3 : DerivedCategoryObject
            The third object in the distinguished triangle

        Raises:
        -------
        TypeError
            If any of the objects are not instances of DerivedCategoryObject
        
        """

        if not isinstance(derived_object1, DerivedCategoryObject):
            raise TypeError("derived_object1 must be an instance of DerivedCategoryObject.")
        if not isinstance(derived_object2, DerivedCategoryObject):
            raise TypeError("derived_object2 must be an instance of DerivedCategoryObject.")
        if not isinstance(derived_object3, DerivedCategoryObject):
            raise TypeError("derived_object3 must be an instance of DerivedCategoryObject.")

        self.object1 = derived_object1
        self.object2 = derived_object2
        self.object3 = derived_object3

        # call on helper method to potentially update chern characters of any
        # DerivedCategoryObjects that have chern character (0,0,0)
        self._update_chern_characters()


    def __str__(self):
        """
        String representation of the distinguished triangle. Since the chain complexes are printed
        horizontally, the distinguished triangle will be printed vertically.

        Returns:
        -------
        str
            A string representation of the distinguished triangle in the form
                A
                |
                |
                |
                V
                B
                |
                |
                |
                V
                C
        """
        ret_str = str(self.object1) 
        ret_str += '\n|\n|\n|\nV\n'
        ret_str += str(self.object2)
        ret_str += '\n|\n|\n|\nV\n'
        ret_str += str(self.object3) 
        return ret_str

    
    def shiftRight(self):
        """
        Method to shift the distinguished triangle to the right. This is done by wrapping the last
        object in the previous triangle around to the first object in the new triangle (homologically 
        shifted by -1) and then letting the previous first object be the new second object. The new
        third object is the previous second object.

        Returns:
        -------
        DistinguishedTriangle
            The distinguished triangle shifted to the right
        """
        return DistinguishedTriangle(self.object3.shiftComplex(-1), self.object1, self.object2)
    def shiftLeft(self):
        """
        Method to shift the distinguished triangle to the left. This is done by wrapping the first
        object in the previous triangle around to the last object in the new triangle (homologically
        shifted by 1) and then letting the previous second object be the new first object. The new
        third object is the previous second object.

        Returns:
        -------
        DistinguishedTriangle
            The distinguished triangle shifted to the left
        """

        return DistinguishedTriangle(self.object2, self.object3, self.object1.shiftComplex(1))

   

    def _update_chern_characters(self):
        """
        Helper function to update the Chern Characters of the objects in the distinguished triangle
        if they are not already set. This is done by using the additivity of the Chern Character on
        exact sequences. In particular, if the Chern Character of the first object is (0,0,0), then
        the Chern Character of the first object is updated to be the difference of the Chern Characters
        of the second and third objects. Similarly, if the Chern Character of the second object is (0,0,0),
        then the Chern Character of the second object is updated to be the sum of the Chern Characters of
        the first and third objects. If the Chern Character of the third object is (0,0,0), then the Chern
        Character of the third object is updated to be the difference of the Chern Characters of the first
        and second objects.

        """

        if self.object1.chernCharacter() == ChernCharacter(0,0,0) and not isinstance(self.object1, ChainComplex):
            # Update the Chern character of the first complex

            chern2 = self.object2.chernCharacter()
            chern3 = self.object3.chernCharacter()

            self.object1.chern_character = chern2 - chern3
        elif self.object2.chernCharacter() == ChernCharacter(0,0,0) and not isinstance(self.object2, ChainComplex):
            # Update the Chern character of the second complex

            chern1 = self.object1.chernCharacter()
            chern3 = self.object3.chernCharacter()

            self.object2.chern_character = chern1 + chern3
        elif self.object3.chernCharacter() == ChernCharacter(0,0,0) and not isinstance(self.object3, ChainComplex):
            # Update the Chern character of the third complex

            chern1 = self.object1.chernCharacter()
            chern2 = self.object2.chernCharacter()
            
            self.object3.chern_character = chern2 - chern1

            
           
              



    







###############################################################################
#                                                                             #
#                       Spherical Twist Composition                           #
# ----------------------------------------------------------------------------#
#  This object is used to represent the composition of spherical twists in    #
#  in the derived category of local P2 consisting of sheaves supported on the #
#  zero divisor P2. In particular, the only compositions of twists we consier #
#  are twists around (pushforwards) of line bundles. Theoretically, since the #
#  derived category of coherent sheaves on P2 is constructible, the braid     #
#  relations between spherical twists allow any object obtained as a series   #
#  of spherical twists to be represented as a sequence of spherical twists    #
#  specifically around line bundles.                                          #
#                                                                             #
#  In order to determine possible Harder-Narasimhan filtrations of the        #
#  spherical twists, we need to iteratively keep track of previous Harder-    #
#  Narasimhan filtrations. This is done by expanding the leaves of a binary   #
#  tree, where each node is a spherical twist, and the children of each node  #
#  are the spherical twists obtained by applying a successive twist.          #
#                                                                             #   
###############################################################################      




class Node:
    def __init__(self, key):
        self.key = key
        self.left = None
        self.right = None


class BinaryFiltrationTree:
    def __init__(self, root_node=None):
        if root_node is not None and not isinstance(root_node, Node):
            raise TypeError("root_node must be an instance of Node.")
        # Private member variable for the root node
        self._root_node = root_node

    def expand_leaves(self, funct1, funct2):
        """
        Traverses the tree and adds two children to every leaf node.
        The left child is generated by funct1() and the right child by funct2().
        """
        self._expand_leaves_recursive(self._root_node, funct1, funct2)

    def _expand_leaves_recursive(self, node, funct1, funct2):
        if node is None:
            return
        # Check if the node is a leaf node (no children)
        if node.left is None and node.right is None:
            node.left = Node(funct1(node.key))
            node.right = Node(funct2(node.key))
        else:
            self._expand_leaves_recursive(node.left, funct1, funct2)
            self._expand_leaves_recursive(node.right, funct1, funct2)



    

class SphericalTwistComposition():

    def __init__(self, integer_array):

        if not integer_array:
            raise ValueError("Array is empty")
        
        self.integer_array = integer_array
            
        self.base_bundle = LineBundle(integer_array[0])
        self.defining_triangle_tree = None

        if len(integer_array) > 1:
            line_bundle_2 = LineBundle(integer_array[1])

            first_twist = self.__sph_twist_LineBundles(line_bundle_2, self.base_bundle)
            shifted_twist = first_twist.shiftLeft()
            print(shifted_twist)
            root_node = Node(first_twist)
            self.defining_triangle_tree = BinaryFiltrationTree(root_node=root_node)


    
    def _dimHom_LineBundles(self, line_bundle_1, line_bundle_2):
        degree_dif = line_bundle_2.deg - line_bundle_1.deg

        if degree_dif == 0:
            return (1, 0, 0, 1)
        elif degree_dif > -3 and degree_dif < 0:
            rank3 = math.comb(line_bundle_1.deg - line_bundle_2.deg + 2, 2)
            return (0, 0, 0, rank3)
        elif degree_dif > 0 and degree_dif < 3:
            rank0 = math.comb(degree_dif + 2, 2)
            return (rank0, 0, 0, 0)
        elif degree_dif >= 3:
            rank0 = math.comb(degree_dif + 2, 2)
            rank1 = math.comb(degree_dif - 1, 2)
            return (rank0, rank1, 0, 0)
        else:
            rank2 = math.comb(line_bundle_1.deg - line_bundle_2.deg -1, 2)
            rank3 = math.comb(line_bundle_1.deg - line_bundle_2.deg + 2, 2)    
            return (0, 0, rank2, rank3)



    def __sph_twist_LineBundles(self, line_bundle_1, line_bundle_2):
        '''
        use triangle  i^* i_* E -> E -> E x O(-3)[2]
        '''

        homDims = self._dimHom_LineBundles(line_bundle_1, line_bundle_2)

        dimension_vector = [] 
        shift_vector = []

        for i in range(len(homDims)):
            
            if homDims[i] == 0:
                continue
            dimension_vector.append(homDims[i])
            shift_vector.append(-1*i)

        object1 = ChainComplex([line_bundle_1, line_bundle_1], dimension_vector, shift_vector)
        object2 = ChainComplex([line_bundle_2], [1], [0])
        object3 = DerivedCategoryObject(string=f"Tw_{line_bundle_1.deg} O({line_bundle_2.deg})")

        return DistinguishedTriangle(object1, object2, object3)




        

if __name__ == "__main__":
    linebundle1 = LineBundle(-3)
    linebundle2 = LineBundle(-2)
    linebundle3 = LineBundle(-1)
    linebundle4 = LineBundle(0)
    linebundle5 = LineBundle(1)


    cot = CotangentBundle(0)
    print(cot)
    print(cot.chernCharacter())
    print(cot.isCotangentBundleSum())

    chern1 = cot.chernCharacter() + linebundle4.chernCharacter() 
    chern2 = 3 * linebundle3.chernCharacter()

    print("Chern 1: {}".format(chern1))
    print("Chern 2: {}".format(chern2))

    print("Chern1 == chern2: {}".format(chern1 == chern2))




