import numpy as np
import numbers
from collections.abc import Sequence


###############################################################################
#                                                                             #
#                           Numerical Characters                              #
# ----------------------------------------------------------------------------#
#  These objects are used to represent Chern classes of coherent sheaves /    #
#  vector bundles / line bundles, as well as more abstract constructions      #
#  in the derived category of coherent sheaves on some variety. As most       #
#  stability conditions will be numerical stability conditions, they will     #
#  only rely on the information of Chern characters as opposed to more        #
#  arbitrary geometric data.                                                  #
#                                                                             #
###############################################################################


class ChernCharacter():
    """!
    Generic parent class for the chern character. Chern characters should operate like
    lists of floating point numbers, since they are graded objects concentrated in degrees 
    going from 0 to the dimension of the variety. 

    Attributes:
    -------------

        graded_element (np.array): The graded element of the Chern Character object. This is a numpy array of floats or integers

    """

    def __init__(self, graded_element):
        r"""!
        Initialize a Chern Character object with the specified graded element. The graded element
        should be a list of floats or integers, representing the Chern character of a coherent sheaf
        or vector bundle. Since Chern characters are additive on exact sequences, this allows us to
        also consider Chern characters of complexes of coherent sheaves.

        In order to allow for quick computations, the graded element is stored as a numpy array.

        \param list graded_element A list of floats or integers representing the Chern character of a coherent sheaf or vector bundle. Theoretically this value should be rational, but this is not enforced in the class due to the additional complexity of enforcing this constraint.

        \throws TypeError If the graded element is not a list of floats or integers
        """

      
       

        if not isinstance(graded_element, (Sequence, np.ndarray)):
            raise TypeError("Chern Character must be initialized with a list of floats or integers")

        if not all(isinstance(x, numbers.Number) for x in graded_element):
            raise TypeError("Chern Character must be initialized with a list of floats or integers")

        self.graded_element = np.array(graded_element) ## Chern Characters stored as numpy array for quick computations

    def __add__(self, other):
        r"""!
        Add two Chern Character objects together. The Chern Character objects must have the same length
        in order to add them together. Since the addition of numpy vectors is already defined, this 
        method simply calls the addition of the numpy arrays.

        \param ChernCharacter other The Chern Character object to add to the current Chern Character object

        \return ChernCharacter The sum of the two Chern Character objects, which is simply the element-wise sum of the graded elements

        \throws TypeError If other is not a ChernCharacter object
        \throws ValueError If the Chern Character objects do not have the same length
        """

        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only add ChernCharacter objects together.")
        if len(self.graded_element) != len(other.graded_element):
            raise ValueError("Chern Characters must have the same length to add them together.")
        
        return ChernCharacter(self.graded_element + other.graded_element)
    
    def __sub__(self, other):
        r"""!
        Subtract two Chern Character objects together. The Chern Character objects must have the same length
        in order to subtract them. Since the subtraction of numpy vectors is already defined, this
        method simply calls the subtraction of the numpy arrays.

        \param ChernCharacter other The Chern Character object to subtract from the current Chern Character object

        \return ChernCharacter The difference of the two Chern Character objects, which is simply the element-wise difference of the graded elements

        \throws TypeError If other is not a ChernCharacter object

        \throws ValueError If the Chern Character objects do not have the same length
        """

        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only subtract ChernCharacter objects together.")
        if len(self.graded_element) != len(other.graded_element):
            raise ValueError("Chern Characters must have the same length to subtract them.")
        
        return ChernCharacter(self.graded_element - other.graded_element)
    
    def __mul__(self, scalar):
        r"""!
        Multiply the Chern Character object by a scalar. The scalar must be an integer, as the Chern Character
        is a graded object. Since the multiplication of numpy vectors is already defined, this method simply
        calls the multiplication of the numpy array by the scalar.

        \param int scalar The integer to multiply the Chern Character object by

        \return ChernCharacter The Chern Character object multiplied by the scalar. This is simply the element-wise multiplication of the graded element by the scalar

        \throws TypeError If scalar is not an integer
        """

        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")
        
        return ChernCharacter(self.graded_element * scalar)
    
    def __rmul__(self, scalar):
        r"""!
        This is effectively the same as the __mul__ method, though allows for the scalar to be on the left side of the multiplication.

        \param int scalar The integer to multiply the Chern Character object by

        \return ChernCharacter The Chern Character object multiplied by the scalar. This is simply the element-wise multiplication of the graded element by the scalar


        \throws TypeError If scalar is not an integer
        """

    
        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")
        
        return ChernCharacter(self.graded_element * scalar)
    
    def __eq__(self, other):
        r"""!
        Check if two Chern Character objects are equal. The Chern Character objects are equal if their graded elements are equal.

        \param ChernCharacter other The Chern Character object to compare to the current Chern Character object

        \return bool True if the Chern Character objects are equal, False otherwise

        """

        if not isinstance(other, ChernCharacter):
            return False
        
        if len(self.graded_element) != len(other.graded_element):
            return False

        return np.all(self.graded_element == other.graded_element)
    
    def __str__(self):
        r"""!
        String representation of the Chern Character object. This is simply the string representation of the numpy array,
        formatted as

        <a_0, a_1, ..., a_n>

        where a_i are the elements of the Chern Character object.

        \return str A string representation of the Chern Character object
        """

       
        formatted_str = ",".join(map(str, self.graded_element))
        return f'<{formatted_str}>'
    
    def __getitem__(self, key):
        r"""!
        Allow for indexing of the Chern Character object. This is simply a pass-through to the numpy array.

        \param int key The index of the Chern Character object to retrieve

        \return float The value of the Chern Character object at the specified index

        \throws TypeError If key is not an integer
        """

        if not isinstance(key, int):
            raise TypeError("Chern Character objects can only be indexed by integers.")

        return self.graded_element[key]
    
    def __len__(self):
        r"""!
        Return the length of the Chern Character object. This is simply the length of the numpy array.

        \return int The length of the Chern Character object
        """

        return len(self.graded_element)
    
    def __iter__(self):
        r"""!
        Allow for iteration over the Chern Character object. This is simply a pass-through to the numpy array.

        \return iter An iterator over the Chern Character object
        """
        return iter(self.graded_element)
    
    def __reversed__(self):
        r"""!
        Return the reversed Chern Character object. This is simply a pass-through to the numpy array.

        \return reversed The reversed Chern Character object
        """
        return reversed(self.graded_element)
    
    def __contains__(self, item):
        r"""!
        Check if an item is in the Chern Character object. This is simply a pass-through to the numpy array.

        \param item float The item to check if it is in the Chern Character object

        \return bool True if the item is in the Chern Character object, False otherwise
        """

 
        return item in self.graded_element
    
    def __setitem__(self, key, value):
        r"""!
        Allow for setting the value of an index in the Chern Character object. This is simply a pass-through to the numpy array.

        \param int key The index of the Chern Character object to set

        \param float value The value to set the Chern Character object at the specified index

        \throws TypeError If key is not an integer
        \throws TypeError If value is not a number
        """
        if not isinstance(key, int):
            raise TypeError("Chern Character objects can only be indexed by integers.")
        if not isinstance(value, numbers.Number):
            raise TypeError("Chern Character objects can only be set with numbers.")

        self.graded_element[key] = value
    





        

