
from . import DerivedCategoryObject
from src.DerivedCategory.GeometryContext import GeometryContext
from src.DerivedCategory.ChernCharacter import ChernCharacter


class NumericalObject(DerivedCategoryObject):

    def __init__(self, chern_char : ChernCharacter, geometry_context : GeometryContext):
        
        
        if not isinstance(geometry_context, GeometryContext):
            raise TypeError("Geometry context must be an instance of GeometryContext")
        
        if not isinstance(chern_char, ChernCharacter):
            raise TypeError("Chern character must be an instance of ChernCharacter")
        
        self.geometry_context = geometry_context ## The geometry context of the derived category object. This is an instance of the GeometryContext class, which contains information about the variety and its intersection form.

        self.chern_character = chern_char ## The Chern character of the derived category object. This is an instance of the ChernCharacter class, which contains information about the Chern character of the object.

    def __str__(self):
        r"""!
        
        """
        return "Object with Chern character: " + str(self.chern_character)
    

    def shift(self, n : int):

        if not isinstance(n, int):
            raise TypeError("Shift must be an integer")

        self.chern_character *= (-1)**n

        return self
    

    
    def chernCharacter(self):
        
        return self.chern_character
    