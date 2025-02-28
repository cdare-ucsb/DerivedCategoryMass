import numpy as np

###############################################################################
#                                                                             #
#                           Numerical Characters                              #
# ----------------------------------------------------------------------------#
#  These objects are used to represent Chern classes of coherent sheaves,     #
#  vector bundles, and line bundles. They are used to determine if a given    #
#  object is a sum of line bundles or cotangent bundles.                      #
#                                                                             #
###############################################################################


class ChernCharacter():
    """
    Generic parent class for the chern character. Chern characters should operate like
    lists of floating point numbers, since they are graded objects concentrated in degrees 
    going from 0 to the dimension of the variety. 
    """

    def __init__(self, graded_element):

        if not all(isinstance(x, (int,float)) for x in graded_element):
            raise TypeError("Chern Character must be initialized with a list of floats or integers")

        self.graded_element = np.array(graded_element)

    def __add__(self, other):
        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only add ChernCharacter objects together.")
        if len(self.graded_element) != len(other.graded_element):
            raise ValueError("Chern Characters must have the same length to add them together.")
        
        return ChernCharacter(self.graded_element + other.graded_element)
    
    def __sub__(self, other):
        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only subtract ChernCharacter objects together.")
        if len(self.graded_element) != len(other.graded_element):
            raise ValueError("Chern Characters must have the same length to subtract them.")
        
        return ChernCharacter(self.graded_element - other.graded_element)
    
    def __mul__(self, scalar):
        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")
        
        return ChernCharacter(self.graded_element * scalar)
    
    def __rmul__(self, scalar):
        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")
        
        return ChernCharacter(self.graded_element * scalar)
    
    def __eq__(self, other):
        if not isinstance(other, ChernCharacter):
            return False
        
        if len(self.graded_element) != len(other.graded_element):
            return False

        return np.all(self.graded_element == other.graded_element)
    
    def __str__(self):
        formatted_str = ",".join(map(str, self.graded_element))
        return f'<{formatted_str}>'
    
    def __getitem__(self, key):
        return self.graded_element[key]
    
    def __len__(self):
        return len(self.graded_element)
    
    def __iter__(self):
        return iter(self.graded_element)
    
    def __reversed__(self):
        return reversed(self.graded_element)
    
    def __contains__(self, item):
        return item in self.graded_element
    
    def __setitem__(self, key, value):
        self.graded_element[key] = value
    





        

