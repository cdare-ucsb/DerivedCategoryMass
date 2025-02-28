from .ChernCharacter import ChernCharacter
import math
import cmath


###############################################################################
#                                                                             #
#                            Single-Degree Objects                            #
# ----------------------------------------------------------------------------#
#  These objects are used to represent coherent sheaves, vector bundles, and  #
#  line bundles. Line bundles (i.e. vector bundles of rank 1) will be our     #
#  building blocks for the sake of this project, since every coherent sheaf   #
#  in projective space  has a resolution by line bundles. Moreover, all       #         
#  spherical twists that we will consider can be represented as a composition #
#  of twists around (the pushforward of) line bundles.                        #    
#                                                                             #
###############################################################################


IMPLEMENTED_CATAGORIES = ['P1', 'P2']



class CoherentSheaf():
    
    def __init__(self, chern_character, catagory):
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {catagory} is not implemented.")
        
        if catagory == 'P1' and len(chern_character) != 2:
            raise ValueError("P1 objects should have a Chern Character of length 2")
        if catagory == 'P2' and len(chern_character) != 3:
            raise ValueError("P2 objects should have a Chern Character of length 3")
        
        if not isinstance(chern_character, ChernCharacter):
            raise TypeError("Chern Character must be a ChernCharacter object")


        self.catagory = catagory
        self.chern_character = chern_character

        

    def chernCharacter(self):
        return self.chern_character

    def phase(self, *args):
        return cmath.phase(self.central_charge(*args)) / math.pi

    def central_charge(self, *args):

        if self.catagory == 'P1':

            # check that args is a single complex number
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")

            return -1*self.chern_character.graded_element[1] + args[0]*self.chern_character.graded_element[0]
        
        elif self.catagory == 'P2':

            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")

            return complex(-1*self.chern_character.graded_element[2] +
                            args[1] * self.chern_character.graded_element[0],
                              self.chern_character.graded_element[1] - args[0] * self.chern_character.graded_element[0])
    
        else:
            # Catagory is not currently implemented
            raise NotImplementedError("Only P1 and P2 catagories are implemented")

    def __str__(self):
        return f'CoherentSheaf with Chern Character {self.chern_character}'


    def is_semistable(self, *args):
        pass




class LineBundle(CoherentSheaf):

    def __init__(self, degree, catagory):
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {catagory} is not implemented.")
        self.degree = degree
        self.catagory = catagory
        self.chern_character = ChernCharacter([1, self.degree, float(self.degree**2)/2])


    def is_semistable(self, *args):
        return True


    def __str__(self):
        return f'O({self.degree})'

