from .DerivedCategoryObject import DerivedCategoryObject




class DistinguishedTriangle():
    """
    A distinguished triangle is a sequence of three objects in the derived category of coherent sheaves
    (in our case) with morphisms of complexes between them such that the third complex is always 
    quasi-isomorphic to the mapping cone of the morphism between the first two complexes. In the derived
    category, distinguished triangles generalize the notion of exact sequences since they also give rise
    to long exact sequences on the level of cohomology.

    In our implementation, we will use the classes DerivedCategoryObject and ChainComplex to represent
    the objects in the distinguished triangle. 

    Attributes:
    ----------
    object1 : DerivedCategoryObject
        The first object in the distinguished triangle
    object2 : DerivedCategoryObject
        The second object in the distinguished triangle
    object3 : DerivedCategoryObject
        The third object in the distinguished triangle

    """

    def __init__(self, derived_object1, derived_object2, derived_object3):
        """
        Initialize an instance of DistinguishedTriangle with the specified objects, and
        then update the Chern Characters of the objects if they are not already set.

        Parameters:
        ----------
        derived_object1 : DerivedCategoryObject
            The first object in the distinguished triangle
        derived_object2 : DerivedCategoryObject
            The second object in the distinguished triangle
        derived_object3 : DerivedCategoryObject
            The third object in the distinguished triangle

        Raises:
        -------
        TypeError
            If any of the objects are not instances of DerivedCategoryObject
        
        """

        if not all(isinstance(obj, DerivedCategoryObject) for obj in [derived_object1, derived_object2, derived_object3]):
            raise TypeError("All objects in the distinguished triangle must be instances of DerivedCategoryObject")

        first_catagory = derived_object1.catagory
        if not all( obj.catagory == first_catagory for obj in [derived_object2, derived_object3]):
            raise TypeError("All objects in the distinguished triangle must be instances of DerivedCategoryObject, from the same underlying category")

        self.object1 = derived_object1
        self.object2 = derived_object2
        self.object3 = derived_object3

        self.catagory = first_catagory

        # call on helper method to potentially update chern characters of any
        # DerivedCategoryObjects that have chern character (0,0,0)
        self._update_chern_characters()


    def __str__(self):
        """
        String representation of the distinguished triangle. Since the chain complexes are printed
        horizontally, the distinguished triangle will be printed vertically.

        Returns:
        -------
        str
            A string representation of the distinguished triangle in the form
                A
                |
                |
                |
                V
                B
                |
                |
                |
                V
                C
        """
        ret_str = str(self.object1) 
        ret_str += '\n|\n|\n|\nV\n'
        ret_str += str(self.object2)
        ret_str += '\n|\n|\n|\nV\n'
        ret_str += str(self.object3) 
        return ret_str

    
    def shiftRight(self):
        """
        Method to shift the distinguished triangle to the right. This is done by wrapping the last
        object in the previous triangle around to the first object in the new triangle (homologically 
        shifted by -1) and then letting the previous first object be the new second object. The new
        third object is the previous second object.

        Returns:
        -------
        DistinguishedTriangle
            The distinguished triangle shifted to the right
        """
        return DistinguishedTriangle(self.object3.shiftComplex(-1), self.object1, self.object2)
    def shiftLeft(self):
        """
        Method to shift the distinguished triangle to the left. This is done by wrapping the first
        object in the previous triangle around to the last object in the new triangle (homologically
        shifted by 1) and then letting the previous second object be the new first object. The new
        third object is the previous second object.

        Returns:
        -------
        DistinguishedTriangle
            The distinguished triangle shifted to the left
        """

        return DistinguishedTriangle(self.object2, self.object3, self.object1.shiftComplex(1))

   

    def _update_chern_characters(self):
        """
        Helper function to update the Chern Characters of the objects in the distinguished triangle
        if they are not already set. This is done by using the additivity of the Chern Character on
        exact sequences. In particular, if the Chern Character of the first object is (0,0,0), then
        the Chern Character of the first object is updated to be the difference of the Chern Characters
        of the second and third objects. Similarly, if the Chern Character of the second object is (0,0,0),
        then the Chern Character of the second object is updated to be the sum of the Chern Characters of
        the first and third objects. If the Chern Character of the third object is (0,0,0), then the Chern
        Character of the third object is updated to be the difference of the Chern Characters of the first
        and second objects.

        """

        if self.object1.chernCharacter() == None:
            # Update the Chern character of the first complex

            chern2 = self.object2.chernCharacter()
            chern3 = self.object3.chernCharacter()

            self.object1.chern_character = chern2 - chern3
        elif self.object2.chernCharacter() == None:
            # Update the Chern character of the second complex

            chern1 = self.object1.chernCharacter()
            chern3 = self.object3.chernCharacter()

            self.object2.chern_character = chern1 + chern3
        elif self.object3.chernCharacter() == None:
            # Update the Chern character of the third complex

            chern1 = self.object1.chernCharacter()
            chern2 = self.object2.chernCharacter()
            
            self.object3.chern_character = chern2 - chern1

