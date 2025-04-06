from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject




class DistinguishedTriangle():
    r"""!
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
    

    def __getitem__(self, index):
        r"""!
        Get the object at the specified index in the distinguished triangle. The index is 0-based
        and can be 0, 1, or 2.

        \param int index: The index of the object to get

        \return DerivedCategoryObject The object at the specified index
        """

        if index == 0:
            return self.object1
        elif index == 1:
            return self.object2
        elif index == 2:
            return self.object3
        else:
            raise IndexError("Index must be 0, 1, or 2")
        
    def __setitem__(self, index, value):
        r"""!
        Set the object at the specified index in the distinguished triangle. The index is 0-based
        and can be 0, 1, or 2.

        \param int index: The index of the object to set
        \param DerivedCategoryObject value: The object to set at the specified index

        \throws TypeError: If the value is not an instance of DerivedCategoryObject
        \throws IndexError: If the index is not 0, 1, or 2
        """

        if not isinstance(value, DerivedCategoryObject):
            raise TypeError("Value must be an instance of DerivedCategoryObject")
        
        if index == 0:
            self.object1 = value
        elif index == 1:
            self.object2 = value
        elif index == 2:
            self.object3 = value
        else:
            raise IndexError("Index must be 0, 1, or 2")
        
    def __eq__(self, other):
        r"""!
        Check if two distinguished triangles are equal. Two distinguished triangles are equal if
        the objects in the triangles are equal (in any order).

        \param DistinguishedTriangle other: The other distinguished triangle to compare to

        \return bool True if the two distinguished triangles are equal, False otherwise
        """

        if not isinstance(other, DistinguishedTriangle):
            return False
        return self.object1 == other.object1 and self.object2 == other.object2 and self.object3 == other.object3

   

   

