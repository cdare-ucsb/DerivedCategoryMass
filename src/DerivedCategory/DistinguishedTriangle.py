from DerivedCategory.DerivedCategoryObject import DerivedCategoryObject




class DistinguishedTriangle():
    """!
    A distinguished triangle is a sequence of three objects in the derived category of coherent sheaves
    (in our case) with morphisms of complexes between them such that the third complex is always 
    quasi-isomorphic to the mapping cone of the morphism between the first two complexes. In the derived
    category, distinguished triangles generalize the notion of exact sequences since they also give rise
    to long exact sequences on the level of cohomology.

    In our implementation, we will use the classes DerivedCategoryObject and ChainComplex to represent
    the objects in the distinguished triangle. 

    """

    def __init__(self, derived_object1 : DerivedCategoryObject, derived_object2 : DerivedCategoryObject, derived_object3 : DerivedCategoryObject):
        r"""!
        Initialize an instance of DistinguishedTriangle with the specified objects, and
        then update the Chern Characters of the objects if they are not already set.

        \param DerivedCategoryObject derived_object1 The first object in the distinguished triangle
        \param DerivedCategoryObject derived_object2: The second object in the distinguished triangle
        \param DerivedCategoryObject derived_object3: The third object in the distinguished triangle

        \throws TypeError: If any of the objects are not instances of DerivedCategoryObject
        \throws TypeError: If the objects are not from the same underlying category
        
        """

        if not all(isinstance(obj, DerivedCategoryObject) for obj in [derived_object1, derived_object2, derived_object3]):
            raise TypeError("All objects in the distinguished triangle must be instances of DerivedCategoryObject")

        self.geometric_context = derived_object1.geometry_context
        if not all( obj.geometry_context != self.geometric_context for obj in [derived_object2, derived_object3]):
            raise TypeError("All objects in the distinguished triangle must be instances of DerivedCategoryObject, from the same underlying geometric_context")

        self.object1 = derived_object1 ## The first object in the distinguished triangle

        self.object2 = derived_object2 ## The second object in the distinguished triangle
        
        self.object3 = derived_object3 ## The third object in the distinguished triangle



    def __str__(self):
        r"""!
        String representation of the distinguished triangle. Since the chain complexes are printed
        horizontally, the distinguished triangle will be printed vertically.

        \return str A string representation of the distinguished triangle
        """

        ret_str = str(self.object1) 
        ret_str += '\n|\n|\n|\nV\n'
        ret_str += str(self.object2)
        ret_str += '\n|\n|\n|\nV\n'
        ret_str += str(self.object3) 
        return ret_str

    
    def rotateRight(self):
        r"""!
        Method to shift the distinguished triangle to the right. This is done by wrapping the last
        object in the previous triangle around to the first object in the new triangle (homologically 
        shifted by -1) and then letting the previous first object be the new second object. The new
        third object is the previous second object.

        \return DistinguishedTriangle The distinguished triangle shifted to the right
        """
        return DistinguishedTriangle(self.object3.shift(-1), self.object1, self.object2)
    

    def rotateLeft(self):
        r"""!
        Method to shift the distinguished triangle to the left. This is done by wrapping the first
        object in the previous triangle around to the last object in the new triangle (homologically
        shifted by 1) and then letting the previous second object be the new first object. The new
        third object is the previous second object.

        \return DistinguishedTriangle The distinguished triangle shifted to the left
        """

        return DistinguishedTriangle(self.object2, self.object3, self.object1.shift(1))
    

    def __hash__(self):

        return hash((self.object1, self.object2, self.object3))

   

   

