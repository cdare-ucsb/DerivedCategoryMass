from ChernCharacter import ChernCharacter
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


IMPLEMENTED_CATAGORIES = ['P1', 'P2', 'K3']



class CoherentSheaf():
    
    def __init__(self, chern_character, catagory):
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {catagory} is not implemented.")
        
        if catagory == 'P1' and len(chern_character) != 2:
            raise ValueError("P1 objects should have a Chern Character of length 2")
        if catagory == 'P2' and len(chern_character) != 3:
            raise ValueError("P2 objects should have a Chern Character of length 3")
        if catagory == 'K3' and len(chern_character) != 3:
            raise ValueError("K3 objects should have a Chern Character of length 3")
        
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
                raise TypeError("P1 central charges should have a single complex parameter")

            return -1*self.chern_character[1] + args[0]*self.chern_character[0]
        
        elif self.catagory == 'P2':

            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 central charges should have two real number parameters")

            return complex(-1*self.chern_character[2] +
                            args[1] * self.chern_character[0],
                              self.chern_character[1] - args[0] * self.chern_character[0])
    
        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")

            alpha = args[0]
            beta = args[1]
            d = args[2]
            
            return complex(2*d*alpha * self.chern_character[1] - self.chern_character[2] - self.chern_character[0] + (beta**2 - alpha**2)*d*self.chern_character[0], 
                           2*d*self.chern_character[1] - 2*d*alpha*beta*self.chern_character[0])
    
        else:
            # Catagory is not currently implemented
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")

    def __str__(self):
        return f'CoherentSheaf with Chern Character {self.chern_character}'


    def is_semistable(self, *args):
        pass

    def __hash__(self):
        return hash(self.chern_character)




class LineBundle(CoherentSheaf):

    def __init__(self, degree, catagory):
        if catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {catagory} is not implemented.")
        if not isinstance(degree, int):
            raise ValueError(f"degree must be an integer: currently passed {type(degree)}")
        
        self.degree = degree
        self.catagory = catagory
        if self.catagory == 'K3' or self.catagory == 'P2':
            self.chern_character = ChernCharacter([1, self.degree, float(self.degree**2)/2])
        elif self.catagory == 'P1':
            self.chern_character = ChernCharacter([1, self.degree])
        


    def is_semistable(self, *args):
        return True


    def __str__(self):
        return f'O({self.degree})'
    
    def __eq__(self, other):
        if not isinstance(other, LineBundle):
            return False
        return self.degree == other.degree
    
    def __hash__(self):
        return hash(self.degree)

