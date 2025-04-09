
from .DerivedCategoryObject import DerivedCategoryObject
from src.DerivedCategory.ChernCharacter import ChernCharacter


class ZeroObject(DerivedCategoryObject):

    def __init__(self, geometry_context):
        r"""!
        Initialize an instance of a ZeroObject in the derived category that is specified by
        the GeometryContext parameter. Since not much else can be said except in more concrete subclasses,
        the only other information that an abstract object may possibly have is a homological shift, which is
        implemented by the shift parameter.

        \param GeometryContext geometry_context The geometry context of the derived category object. This is an instance of the GeometryContext class, which contains information about the variety and its intersection form.
        """
        
        super().__init__(geometry_context)

    def __str__(self):
        r"""!
        String representation of the zero object. Since the zero object is not a real object, it is
        represented by the string "0".
        
        \return str A string representation of the zero object
        """
        return "0"
    

    def shift(self, _: int):

        return self
    

    
    def chernCharacter(self):
        r"""!
        Method to return the Chern Character of the zero object. 
        
        \return ChernCharacter The Chern Character of the zero object
        \note The Chern character of the zero object is always zero, regardless of the geometry context.
        """
        return ChernCharacter(0,
                            dimension=self.geometry_context.divisor_data.variety_dimension,
                            basis=self.geometry_context.divisor_data.basis)
    
