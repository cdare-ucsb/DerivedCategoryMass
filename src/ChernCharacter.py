###############################################################################
#                                                                             #
#                           Numerical Characters                              #
# ----------------------------------------------------------------------------#
#  These objects are used to represent Chern classes of coherent sheaves,     #
#  vector bundles, and line bundles. They are used to determine if a given    #
#  object is a sum of line bundles or cotangent bundles.                      #
#                                                                             #
###############################################################################



############################
#  Generic Parent Class    #
############################


class ChernCharacter():
    """
    Generic parent class for the chern character. Chern characters should operate like
    lists of floating point numbers, since they are graded objects concentrated in degrees 
    going from 0 to the dimension of the variety. 
    """

    def __init__(self):
        pass




############################
#     P1 Implementation    #
############################

class ChernCharacterP1(ChernCharacter):

    def __init__(self, ch0, ch1):
        self.ch0 = ch0
        self.ch1 = ch1

    def __str__(self):
        """
        String representation of the Chern Character class

        Returns:
        -------
        str
            A string representation of the Chern Character
        """
        return f'<{self.ch0}, {self.ch1}>'
    
    def __add__(self, other):
        """
        Method to add two Chern Characters together. This is done by adding the corresponding
        components of the Chern Character.

        Parameters:
        ----------
        other : ChernCharacterP1
            The Chern Character to add to the current Chern Character

        Returns:
        -------
        ChernCharacterP1
            The sum of the two Chern Characters
        """

        if not isinstance(other, ChernCharacterP1):
            raise TypeError("Can only add ChernCharacterP1 objects together.")

        return ChernCharacterP2(self.ch0 + other.ch0, self.ch1 + other.ch1)
    

    def __sub__(self, other):
        """
        Method to subtract two Chern Characters. This is done by subtracting the corresponding
        components of the Chern Character.

        Parameters:
        ----------
        other : ChernCharacterP1
            The Chern Character to subtract from the current Chern Character
        
        Returns:
        -------
        ChernCharacterP1
            The difference of the two Chern Characters
        """

        if not isinstance(other, ChernCharacterP1):
            raise TypeError("Can only subtract ChernCharacterP1 objects together.")

        return ChernCharacterP2(self.ch0 - other.ch0, self.ch1 - other.ch1)
    
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
        ChernCharacterP1
            The Chern Character multiplied by the scalar
        """

        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacterP1 objects by integers.")

        return ChernCharacterP2(self.ch0 * scalar, self.ch1 * scalar)
    
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
        ChernCharacterP1
            The Chern Character multiplied by the scalar
        """

        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacterP1 objects by integers.")

        return ChernCharacterP2(self.ch0 * scalar, self.ch1 * scalar)
    
    def __eq__(self, other):
        """
        Method to determine if two Chern Characters are equal. This is done by checking if
        the corresponding components of the Chern Character are equal.

        Parameters:
        ----------
        other : ChernCharacterP1
            The Chern Character to compare to the current Chern Character

        Returns:
        -------
        bool
            True if the two Chern Characters are equal, False otherwise
        """

        if not isinstance(other, ChernCharacterP1):
            return False

        return self.ch0 == other.ch0 and self.ch1 == other.ch1
    



############################
#     P2 Implementation    #
############################

class ChernCharacterP2(ChernCharacter):
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
        Initialize an instance of ChernCharacterP2 with the specified characteristic classes
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
        other : ChernCharacterP2
            The Chern Character to add to the current Chern Character

        Returns:
        -------
        ChernCharacterP2
            The sum of the two Chern Characters
        """

        if not isinstance(other, ChernCharacterP2):
            raise TypeError("Can only add ChernCharacterP2 objects together.")

        return ChernCharacterP2(self.ch0 + other.ch0, self.ch1 + other.ch1, self.ch2 + other.ch2)
    
    def __sub__(self, other):
        """
        Method to subtract two Chern Characters. This is done by subtracting the corresponding
        components of the Chern Character.

        Parameters:
        ----------
        other : ChernCharacterP2
            The Chern Character to subtract from the current Chern Character
        
        Returns:
        -------
        ChernCharacterP2
            The difference of the two Chern Characters
        """

        if not isinstance(other, ChernCharacterP2):
            raise TypeError("Can only subtract ChernCharacterP2 objects together.")

        return ChernCharacterP2(self.ch0 - other.ch0, self.ch1 - other.ch1, self.ch2 - other.ch2)
    
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
        ChernCharacterP2
            The Chern Character multiplied by the scalar
        """

        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacterP2 objects by integers.")

        return ChernCharacterP2(self.ch0 * scalar, self.ch1 * scalar, self.ch2 * scalar)
    
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
        ChernCharacterP2
            The Chern Character multiplied by the scalar
        """

        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacterP2 objects by integers.")

        return ChernCharacterP2(self.ch0 * scalar, self.ch1 * scalar, self.ch2 * scalar)
    
    def __eq__(self, other):
        """
        Method to determine if two Chern Characters are equal. This is done by checking if
        the corresponding components of the Chern Character are equal.

        Parameters:
        ----------
        other : ChernCharacterP2
            The Chern Character to compare to the current Chern Character

        Returns:
        -------
        bool
            True if the two Chern Characters are equal, False otherwise
        """

        if not isinstance(other, ChernCharacterP2):
            return False

        return self.ch0 == other.ch0 and self.ch1 == other.ch1 and self.ch2 == other.ch2
    
    def isLineBundleP2Sum(self):
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
    
    def isCotangentBundleP2Sum(self):
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
        new_bundle = ChernCharacterP2(self.ch0 + num_copies,
                                self.ch1 + num_copies * original_deg,
                                  self.ch2 + (num_copies * original_deg**2 / 2))
        # Account for possible shift which makes everything negative
        shifted_bundle = ChernCharacterP2(-self.ch0 - num_copies, 
                                    -self.ch1 - num_copies * original_deg,
                                    -self.ch2 - (num_copies * original_deg**2 / 2))    

        
        
        new_bundle_is_cot =  (new_bundle.isLineBundleP2Sum() and new_bundle.ch0 % 3 == 0)
        shifted_bundle_is_cot = (shifted_bundle.isLineBundleP2Sum() and shifted_bundle.ch0 % 3 == 0)
        return new_bundle_is_cot or shifted_bundle_is_cot
    

