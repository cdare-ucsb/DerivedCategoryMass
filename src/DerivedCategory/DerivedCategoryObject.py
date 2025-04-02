from abc import ABC, abstractmethod


from dotenv import load_dotenv
import os
import cmath
import math

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']


class DerivedCategoryObject(ABC):
    r"""!
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

    """
    def __init__(self, catagory : str, shift : int = 0):
        r"""!
        Initialize an instance of the DerivedCategoryObject with the specified catagory, string representation,
        and Chern Character. The Chern Character is optional and can be set later. 

        \param str catagory The catagory of the derived category object. This should be one of the implemented catagories ['P1', 'P2', 'K3'].

        \throws NotImplementedError If the catagory is not implemented
        \throws TypeError If the Chern Character is not an instance of Chern

        """
        
        if catagory not in IMPLEMENTED_CATAGORIES:    
            raise NotImplementedError(f"Catagory {catagory} is not implemented.")
        
        self.catagory = catagory ## The catagory of the derived category object. This should be one of the implemented catagories ['P1', 'P2', 'K3'].
        
        self.shift = shift ## The homological shift of the derived category object. This is an integer that represents the shift of the object in the derived category. The default value is 0, which means that the object is not shifted.
        

    def __str__(self):
        r"""!
        String representation of the derived category object

        \return str A string representation of the derived category object
        """
        return "0"
    

    @abstractmethod 
    def chernCharacter(self):
        r"""!
        Method to return the Chern Character of the derived category object

        \return ChernCharacter The Chern Character of the derived category object
        """

        pass
    

    @abstractmethod
    def shift(self, shift : int):
        r"""!
        Method to shift the derived category object by a given homological shift

        \param int shift The homological shift to apply to the derived category object

        \return DerivedCategoryObject The derived category object shifted by the homological shift
        """

        pass

    def central_charge(self, *args) -> complex:
        r"""!
        Compute the central charge of an object in the derived category of coherent sheaves. For all the current categories
        implemented, the only stability conditions considered are numerical stability conditions; in particular, they only
        depend on the Chern character of the object. Since DerivedCategoryObjects are the highest level objects which 
        have a chernCharacter method, this method will be implemented here (the implementation does not change for any of 
        the children classes).

        \param tuple args The arguments required to compute the central charge. The number of arguments and the type
                     of arguments will depend on the catagory of the sheaves in the complex. For P1, the central
                     charge requires a single complex number. For P2, the central charge requires two floating-point
                     numbers. For K3, the central charge requires two floating-point numbers and an integer.

        \return complex The central charge of the chain complex as a complex number

        \throws ValueError If the number of arguments is incorrect
        \throws TypeError If the type of the arguments is incorrect
        \throws NotImplementedError If the catagory of the sheaves in the complex is not implemented
        """

        
        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge for P1 requires exactly one argument.")
            if not isinstance(args[0], complex):
                raise TypeError("Central charge for P1 requires a complex number as an argument.")

            ch = self.chernCharacter()
            return complex(-1*ch[1] + args[0] * ch[0])


        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge for P2 requires exactly two arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Central charge for P2 requires two floating-point numbers as arguments.")
            
            ch = self.chernCharacter()
            return complex(-1*ch[2] + args[1] * ch[0],
                            ch[1] - args[0] * ch[0])
        
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
            ch = self.chernCharacter()
            
            return complex(2*d*alpha * ch[1] - ch[2] - ch[0] + (beta**2 - alpha**2)*d*ch[0], 
                           2*d*ch[1] - 2*d*alpha*beta*ch[0])

        else:
            raise NotImplementedError("Central charge not implemented for this variety.")
    

    def phase(self, *args) -> float:
        r"""!
        Computes the phase of the central charge of the coherent sheaf. The central charge
        is an element of the dual of the numerical Grothendieck group; in other words, a 
        funtction

        Z : K -> C

        where K is the numerical Grothendieck group, and C is the complex numbers. The phase
        of the central charge is the argument of this complex number.

        \param *args: float or int
            The parameters of the central charge. The number of parameters should be equal
            to the number of parameters required by the central charge for the given catagory.
            For example, a P1 object requires a single complex number parameter, while a P2
            object requires two real number parameters.

        \return float The phase of the central charge of the coherent sheaf, in units of pi
        """

        return cmath.phase(self.central_charge(*args)) / math.pi
    
        
    @abstractmethod
    def is_semistable(self, *args) -> bool:
        r"""!
        Method to determine if the derived category object is semistable with respect to a given stability condition. 
        This will simply act as a wrapper for the central charge method, which should be implemented in child classes.

        \param tuple args 
            The parameters of the stability condition. The number of parameters will depend on the catagory of the object.
            For P1, this will be a single complex number. For P2, this will be two real numbers. For K3, this will be
            two real numbers and one integer.

        \return bool True if the object is semistable with respect to the stability condition, False otherwise

        \throws ValueError
            If the DerivedCategoryObject is not initialized
            If the number of parameters is incorrect for the catagory
        \throws TypeError
            If the parameters are not of the correct type

        """

        pass

    @abstractmethod
    def mass(self, *_) -> float :
        r"""!
        Method to compute the mass of the derived category object with respect to a given stability condition. 
        This will simply act as a wrapper for the central charge method, which should be implemented in child classes.

        \param tuple args 
            The parameters of the stability condition. The number of parameters will depend on the catagory of the object.
            For P1, this will be a single complex number. For P2, this will be two real numbers. For K3, this will be
            two real numbers and one integer.

        \return float The mass of the derived category object with respect to the stability condition

        \throws ValueError
            If the DerivedCategoryObject is not initialized
            If the number of parameters is incorrect for the catagory
        \throws TypeError
            If the parameters are not of the correct type

        """

        pass



    def __hash__(self) -> int:
        r"""!
        Hash function for the derived category object. This is used to create a unique identifier for the object
        in dictionaries and sets. The hash is computed based on the catagory and string representation of the object.

        \return int The hash of the derived category object
        """
        return hash(self.chernCharacter())

    

    # def update_string_by_shift(self, my_str, n):
    #     r"""!
    #     Static helper function to update the possible string representations of the abstract
    #     DerivedCategoryObject by a homological shift. For example, if an object is called A[3],
    #     then a shift of 2 will yield A[5] (as a string). If the object is unshifted to start with,
    #     then shifting the object should return A[n] (as a string). The only object that should not
    #     be shifted is the zero object, which is represented by the string "0".

    #     \param str my_str
    #         The string representation of the derived category object
    #     \param int n
    #         The homological shift to apply to the derived category object

    #     \return str
    #         The string representation of the derived category object shifted by the homological shift

    #     \throws TypeError
    #         If my_str is not a string
    #         If n is not an integer

    #     """

    #     if not isinstance(my_str, str):
    #         raise TypeError("my_str must be a string.")
    #     if not isinstance(n, int):
    #         raise TypeError("n must be an integer.")

    #     # If the object is zero, it should remain zero
    #     if my_str == "0":
    #         return "0"
        
    #     # This regex checks for a pattern at the end of the string that looks like "[number]"
    #     pattern = r'\[(\d+)\]$'
    #     match = re.search(pattern, my_str)
        
    #     if match:
    #         # Extract the current number k, add n, and format the new number.
    #         k = int(match.group(1))
    #         new_k = k + n
    #         # Replace the [k] at the end with [new_k]
    #         updated_str = re.sub(pattern, f'[{new_k}]', my_str)
    #         return updated_str
    #     else:
    #         # If there's no [number] at the end, append [n]
    #         return my_str + f'[{n}]'
            


    


    

    
