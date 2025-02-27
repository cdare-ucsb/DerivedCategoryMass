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



############################
#  Generic Parent Class    #
############################


class ChernCharacter():
    """
    Generic parent class for the chern character. Chern characters should operate like
    lists of floating point numbers, since they are graded objects concentrated in degrees 
    going from 0 to the dimension of the variety. 
    """

    def __init__(self, chern_nums):

        if not all(isinstance(x, (int,float)) for x in chern_nums):
            raise TypeError("Chern Character must be initialized with a list of floats or integers")

        self.chern_nums = np.array(chern_nums)

    def __add__(self, other):
        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only add ChernCharacter objects together.")
        if len(self.chern_nums) != len(other.chern_nums):
            raise ValueError("Chern Characters must have the same length to add them together.")
        
        return ChernCharacter(self.chern_nums + other.chern_nums)
    
    def __sub__(self, other):
        if not isinstance(other, ChernCharacter):
            raise TypeError("Can only subtract ChernCharacter objects together.")
        if len(self.chern_nums) != len(other.chern_nums):
            raise ValueError("Chern Characters must have the same length to subtract them.")
        
        return ChernCharacter(self.chern_nums - other.chern_nums)
    
    def __mul__(self, scalar):
        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")
        
        return ChernCharacter(self.chern_nums * scalar)
    
    def __rmul__(self, scalar):
        if not isinstance(scalar, int):
            raise TypeError("Can only multiply ChernCharacter objects by integers.")
        
        return ChernCharacter(self.chern_nums * scalar)
    
    def __eq__(self, other):
        if not isinstance(other, ChernCharacter):
            return False
        
        if len(self.chern_nums) != len(other.chern_nums):
            return False

        return np.all(self.chern_nums == other.chern_nums)
    
    def __str__(self):
        formatted_str = ",".join(map(str, self.chern_nums))
        return f'<{formatted_str}>'
    




