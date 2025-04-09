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
        Initialize an instance of a DerivedCategoryObject in the derived category that is specified by
        the GeometryContext parameter. Since not much else can be said except in more concrete subclasses,
        the only other information that an abstract object may possibly have is a homological shift, which is
        implemented by the shift parameter.

        \param GeometryContext geometry_context The geometry context of the derived category object. This is an instance of the GeometryContext class, which contains information about the variety and its intersection form.

        \param int shift The homological shift of the derived category object. This is an integer that represents the shift of the object in the derived category. The default value is 0, which means that the object is not shifted. 

        \throws TypeError if the geometry_context is not an instance of GeometryContext or if the shift is not an integer

        """
        
        if not isinstance(geometry_context, GeometryContext):
            raise TypeError("Geometry context must be an instance of GeometryContext")
        
        if not isinstance(shift, int):
            raise TypeError("shift must be an integer")
        

        self.geometry_context = geometry_context ## The geometry context of the derived category object. This is an instance of the GeometryContext class, which contains information about the variety and its intersection form.  
        
        self.shift = shift ## The homological shift of the derived category object. This is an integer that represents the shift of the object in the derived category. The default value is 0, which means that the object is not shifted.
        
    @abstractmethod
    def __str__(self):
        r"""!
        String representation of the derived category object.
        This method should be implemented in subclasses to provide a string representation of the derived category object.

        \return str A string representation of the derived category object
        """
        pass
    

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
        in dictionaries and sets. The hash is computed based off the Chern character, geometry context, and shift.

        NOTE: While the DerivedCatagory object is an ABC abstract class (meaning it should not be instantiated),
        python allows for subclasses to only partially implement all the abstract methods which is not checked until
        runtime. Hypothetically, a class could only implement shift() and not chernCharacter(). Thus, this method
        handles the potential NotImplementedError that could be raised by the chernCharacter() method. However, this
        may be overkill and is worth addressing in future versions.

        \return int The hash of the derived category object
        """
        try:
            ch = self.chernCharacter()
        except NotImplementedError:
            ch = None

        return hash((ch, self.geometry_context, self.shift))
    
    def __eq__(self, other):
        r"""!
        Equality operator for the derived category object. This is used to compare two derived category objects
        to see if they are equal.

        NOTE: While the DerivedCatagory object is an ABC abstract class (meaning it should not be instantiated),
        python allows for subclasses to only partially implement all the abstract methods which is not checked until
        runtime. Hypothetically, a class could only implement shift() and not chernCharacter(). Thus, this method
        handles the potential NotImplementedError that could be raised by the chernCharacter() method. However, this
        may be overkill and is worth addressing in future versions.

        \param DerivedCategoryObject other The other derived category object to compare to

        \return bool True if the derived category objects are equal, False otherwise
        """
        if not isinstance(other, DerivedCategoryObject):
            return False

        # Check if both raise NotImplementedError
        try:
            self_ch = self.chernCharacter()
        except NotImplementedError:
            self_ch = None

        try:
            other_ch = other.chernCharacter()
        except NotImplementedError:
            other_ch = None

        # If one is implemented but the other isn't → not equal
        if (self_ch is None) != (other_ch is None):
            return False

        # If neither implements chernCharacter, fallback to geometry + shift only
        if self_ch is None and other_ch is None:
            return self.geometry_context == other.geometry_context and self.shift == other.shift

        # Otherwise, compare everything
        return (self.geometry_context == other.geometry_context and
                self.shift == other.shift and
                self_ch == other_ch)



    


    

    
