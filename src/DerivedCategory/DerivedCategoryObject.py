from abc import ABC, abstractmethod

from src.DerivedCategory.GeometryContext import GeometryContext


from dotenv import load_dotenv
import os

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
    def __init__(self, geometry_context : GeometryContext , shift : int = 0):
        r"""!
        Initialize an instance of the DerivedCategoryObject with the specified catagory, string representation,
        and Chern Character. The Chern Character is optional and can be set later. 


        

        \throws TypeError If the Chern Character is not an instance of Chern

        """
        
        if not isinstance(geometry_context, GeometryContext):
            raise TypeError("Geometry context must be an instance of GeometryContext")
        
        if not isinstance(shift, int):
            raise TypeError("shift must be an integer")
        

        self.geometry_context = geometry_context ## The geometry context of the derived category object. This is an instance of the GeometryContext class, which contains information about the variety and its intersection form.  
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

    def catagory(self) -> str:
        return self.geometry_context.catagory



    def __hash__(self) -> int:
        r"""!
        Hash function for the derived category object. This is used to create a unique identifier for the object
        in dictionaries and sets. The hash is computed based on the catagory and string representation of the object.

        \return int The hash of the derived category object
        """
        return hash((self.chernCharacter(), self.geometry_context.catagory, self.shift))

    


    


    

    
