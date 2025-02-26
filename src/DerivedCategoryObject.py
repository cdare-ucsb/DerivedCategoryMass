from ChernCharacter import ChernCharacterP2, ChernCharacterP1, ChernCharacter
import re



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
    def __init__(self, string = "0", chern_character = None):
        self.string = string
        if chern_character and not isinstance(chern_character, ChernCharacter):
            raise TypeError("chern_character must be an instance of ChernCharacterP1 or ChernCharacterP2.")
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
        ChernCharacterP1
            The Chern Character of the derived category object
        """
        return self.chern_character
    
    def shiftObject(self, shift):
        """
        Method to shift the derived category object by a given homological shift

        Parameters:
        ----------
        shift : int
            The homological shift to apply to the derived category object

        Returns:
        -------
        DerivedCategoryObjectP1
            The derived category object shifted by the homological shift
        """
        new_string = _update_string_by_shift(self.string, shift)
        if self.chern_character:
            new_chern = int((-1)**shift) * self.chern_character 
        else:
            new_chern = None

        return DerivedCategoryObject(new_string, new_chern)
    



class DerivedCategoryObjectP1():

    def __init__(self, string = "0", chern_character = None):
        self.string = string
        if chern_character and not isinstance(chern_character, ChernCharacterP1):
            raise TypeError("chern_character must be an instance of ChernCharacterP1.")
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
        ChernCharacterP1
            The Chern Character of the derived category object
        """
        return self.chern_character
    

    def shiftObject(self, shift):
        """
        Method to shift the derived category object by a given homological shift

        Parameters:
        ----------
        shift : int
            The homological shift to apply to the derived category object

        Returns:
        -------
        DerivedCategoryObjectP1
            The derived category object shifted by the homological shift
        """
        new_string = _update_string_by_shift(self.string, shift)
        if self.chern_character:
            new_chern = int((-1)**shift) * self.chern_character 
        else:
            new_chern = None

        return DerivedCategoryObjectP1(new_string, new_chern)
    
    def central_charge(self, w):
        """
        Method to compute the central charge of the chain complex. The central charge of a chain complex
        is the alternating sum of the central charges of the individual sheaves in the complex. Since the
        central charge is additive, we may multiply the central charges by the dimension of the sheaf to
        represent direct sums of sheaves.

        Parameters:
        ----------
        w : complex
            The complex parameter for a geometric stability condition in P1
        

        Returns:
        -------
        complex
            The central charge of the chain complex as a complex number

        Raises:
        -------
        TypeError
            If w is not a complex number

        """
        
        if not isinstance(w, complex):
            raise TypeError("w must be a complex number.")
        
        if self.chern_character:
            raise ValueError("DerivedCategoryObject not initialized, cannot compute central charge.")
    
        return complex(-self.chern_character.ch1 + w * self.chern_character.ch0)




class DerivedCategoryObjectP2():
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
    chern_character : ChernCharacterP2
        The Chern Character of the derived category object
    """


    def __init__(self, string = "0", chern_character = None):
        """
        Initialize an instance of DerivedCategoryObject with the specified string and Chern Character.

        Parameters:
        ----------
        string : str
            The string representation of the derived category object
        chern_character : ChernCharacterP2
            The Chern Character of the derived category object

        Raises:
        -------
        TypeError
            If chern_character is not an instance of ChernCharacterP2
        """
        self.string = string

        if chern_character and not isinstance(chern_character, ChernCharacterP2):
            raise TypeError("chern_character must be an instance of ChernCharacterP2.")
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
        ChernCharacterP2
            The Chern Character of the derived category object
        """
        return self.chern_character
    

    def central_charge(self, s, q):
        """
        Method to compute the central charge of the chain complex. The central charge of a chain complex
        is the alternating sum of the central charges of the individual sheaves in the complex. Since the
        central charge is additive, we may multiply the central charges by the dimension of the sheaf to
        represent direct sums of sheaves.

        Parameters:
        ----------
        s : float
            The parameter controlling the imaginary part of the central charge
        q : float
            The parameter controlling the real part of the central charge

        Returns:
        -------
        complex
            The central charge of the chain complex as a complex number

        Raises:
        -------
        TypeError
            If s or q are not floating-point decimals

        """
        
        if not isinstance(s, float) or not isinstance(q, float):
            raise TypeError("s and q must be floating-point decimals.")
        
        if self.chern_character:
            raise ValueError("DerivedCategoryObject not initialized, cannot compute central charge.")
    
        return complex(-self.chern_character.ch2 + q * self.chern_character.ch0, self.chern_character.ch1 - s * self.chern_character.ch0)


    def shiftObject(self, shift):
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
        new_string = _update_string_by_shift(self.string, shift)
        if self.chern_character:
            new_chern = int((-1)**shift) * self.chern_character 
        else:
            new_chern = None

        return DerivedCategoryObjectP2(new_string, new_chern)





###############################################################################
#                         Static Helper Functions                             #
###############################################################################
    

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