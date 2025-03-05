from .ChernCharacter import ChernCharacter
import re



IMPLEMENTED_CATAGORIES = ['P1', 'P2', 'K3']



class DerivedCategoryObject():
    """
    This class acts as a general parent class for objects in the derived category of coherent sheaves.
    The derived category is a triangulated category (though in many cases is a dg-category) that is
    constructed over an abelian category --- in geometric contexts, this is usually the category of
    coherent sheaves on a variety. The derived category is a way to encode the information of the
    homological data of the abelian category in a more structured way, though it is considerably abstract
    since an element of the derived complex technically represents all possible resolutions of the object.

    One can typically think of objects in the derived catagory as chain complexes of sheaves; however, it
    will be true in general that multiple chain complexes of varying lengths can represent the same object
    in the derived category. For example, one can consider the cotangent bundle Ω^1 on a variety IP^2; there
    is a minimal resolution of the cotangent bundle by:

                               0 -> O(-3) -> O(-2)^3 -> Ω^1 -> 0

    obtained by replacing the last two terms in the standard Koszul resolution by Ω^1 = ker(O(-1)^3 --> O). It then
    follows that the derived category object representing Ω^1 is the same as the two-term complex O(-3) -> O(-2)^3. 


    As this is the most general class used in the current program, it will assume the least amount of information
    given since we can only assume the object fits into at least one distinguished triangle which may help
    determine the objects Chern character / numerical data. Thus, the bare minimum information this
    should encode is the catagory of the object and a string representation of the object. 

    
    Attributes:
    ----------
    catagory : str
        The catagory of the derived category object. This should be one of the implemented catagories
        ['P1', 'P2', 'K3'].
    string : str
        A string representation / label of the derived category object.
    chern_character : ChernCharacter
        The Chern Character of the derived category object. This is optional and can be set later, especially
        when the object is put into a distinguished triangle with two objects whose Chern characters are known.

    """
    def __init__(self, catagory, string = "0", chern_character = None):
        """
        Initialize an instance of the DerivedCategoryObject with the specified catagory, string representation,
        and Chern Character. The Chern Character is optional and can be set later. 

        Parameters:
        ----------
        catagory : str
            The catagory of the derived category object. This should be one of the implemented catagories
            ['P1', 'P2', 'K3'].
        string : str
            A string representation / label of the derived category object. Default is "0".
        chern_character : ChernCharacter
            The Chern Character of the derived category object. Default is None.

        Raises:
        -------
        NotImplementedError
            If the catagory is not implemented
        TypeError
            If the Chern Character is not an instance of Chern

        """
        
        if catagory not in IMPLEMENTED_CATAGORIES:    
            raise NotImplementedError(f"Catagory {catagory} is not implemented.")

        if chern_character and not isinstance(chern_character, ChernCharacter):
            raise TypeError("chern_character must be an instance of ChernCharacter.")
        
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
        new_string = self.update_string_by_shift(self.string, shift)
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
        """
        Method to determine if the derived category object is semistable with respect to a given stability condition. 
        This will simply act as a wrapper for the central charge method, which should be implemented in child classes.

        Parameters:
        ----------
        *args : list
            The parameters of the stability condition. The number of parameters will depend on the catagory of the object.
            For P1, this will be a single complex number. For P2, this will be two real numbers. For K3, this will be
            two real numbers and one integer.

        Returns:
        -------
        bool
            True if the object is semistable with respect to the stability condition, False otherwise
        """
        pass

    

    def update_string_by_shift(self, my_str, n):
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
            


    


    

    
