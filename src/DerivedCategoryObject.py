from ChernCharacter import ChernCharacter
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


IMPLEMENTED_CATAGORIES = ['P1', 'P2', 'K3']



class DerivedCategoryObject():
    def __init__(self, catagory, string = "0", chern_character = None):
        
        if catagory not in IMPLEMENTED_CATAGORIES:    
            raise ValueError(f"Catagory {catagory} is not implemented.")

        if chern_character and not isinstance(chern_character, ChernCharacter):
            raise TypeError("chern_character must be an instance of ChernCharacterP1 or ChernCharacterP2.")
        
        self.catagory = catagory
        self.chern_character = chern_character
        self.string = string

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
    
    def shift(self, shift):
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

        return DerivedCategoryObject(catagory=self.catagory, string=new_string, chern_character=new_chern)
    


    
    def central_charge(self, *args):
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

        if self.chern_character is None:
            raise ValueError("DerivedCategoryObject not initialized, cannot compute central charge.")
        
        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("w must be a complex number.")
            
            chern_char = self.chern_character
            return complex(-chern_char[1] + args[0]*chern_char[0])
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("s and q must be floating-point decimals.")
            
            chern_char = self.chern_character
            return complex(-chern_char[2] + args[1] * chern_char[0], chern_char[1] - args[0] * chern_char[0])
        elif self.catagory == 'K3':

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
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        
    def is_semistable(self, *args):
        pass
        





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