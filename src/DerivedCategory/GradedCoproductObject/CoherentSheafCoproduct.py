from .GradedCoproductObject import GradedCoproductObject
from src.DerivedCategory.CoherentSheaf import CoherentSheaf

import math

from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']


class CoherentSheafCoproduct(GradedCoproductObject):
    

    def __init__(self, sheaf_vector, shift_vector, dimension_vector = None):
        r"""!
        Initialize an instance of CoherentSheafCoproduct with the specified sheaf vector, shift vector,
        and potentially a dimension vector. If a dimension vector is not provided, it must 
        consist of non-negative integer values


        \param list sheaf_vector A list of coherent sheaves in the complex
        \param list shift_vector A list of homological shifts in the complex
        \param list dimension_vector A list of the number of direct sums of each coherent sheaf in the complex
        """

        # The main functionality is implemented in the parent class
        super().__init__(sheaf_vector, shift_vector, dimension_vector)

        # Check that the sheaf vector is valid and contains sheaves of only a single catagory
        if not all(isinstance(sheaf, CoherentSheaf) for sheaf in sheaf_vector):
            raise TypeError("All elements of sheaf_vector must be instances of CoherentSheaf.")

    
    


        

    


    

    def isShiftOfSheaf(self):
        r"""!
        Simple helper function which checks if the complex is a shift of a single sheaf

        \return bool True if the complex is a shift of a single sheaf, False otherwise
        """

        return len(self.sheaf_vector) == 1



        
    
    def get_smallest_phase(self, *args):
        r"""!
        Method to compute the smallest phase of the chain complex. This behaves as a sort of "smallest
        Harder-Narasimhan factor" for the complex, since Chain complexes will almost never be stable when
        they have objects in distinct shifts. The phase of an individual element of a chain complex generally
        requires that object to be stable, so that we typically use LineBundles for our current applications. 
        By definition of a slicing, the shift of each object in the complex should add to the respective phases;
        thus, this method computes the smallest sum of the phase of the sheaf and the shift of the sheaf in the
        complex.

        \param tuple args The arguments required to compute the phase. The number of arguments and the type of arguments will
                     depend on the catagory of the sheaves in the complex. For P1, the phase requires a single complex number.
                     For P2, the phase requires two floating-point numbers. For K3, the phase requires two floating-point
                     numbers and an integer.

        \return float The smallest phase of the chain complex

        \throws ValueError If the number of arguments is incorrect
        \throws TypeError If the type of the arguments is incorrect
        \throws NotImplementedError If the catagory of the sheaves in the complex is not implemented
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Phase of P1 requires exactly one argument.")
            if not isinstance(args[0], complex):
                raise TypeError("Phase of P1 requires a complex number as an argument.")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Phase of P2 requires exactly two arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Phase of P2 requires two floating-point numbers as arguments.")
        elif self.catagory == 'K3':
            if len(args) != 3:
                raise ValueError("Phase of K3 requires exactly three arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Phase of K3 requires three floating-point numbers as arguments.")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer.")
        else:
            raise NotImplementedError("Phase not implemented for this variety.")

         # Zip the three lists together and sort by descending shift
        min_shift = min(self.shift_vector)

        bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
        bundles_min_shift = filter(lambda x: x[2] == min_shift or x[2] == min_shift + 1 , bundles)

        min_phase = math.inf

        for sheaf, dim, shift in bundles_min_shift:
            if sheaf.phase(*args) + shift < min_phase:
                min_phase = sheaf.phase(*args) + shift

        return min_phase
    
    def get_largest_phase(self, *args):
        r"""!
        Method to compute the largest phase of the chain complex. This behaves as a sort of "largest
        Harder-Narasimhan factor" for the complex, since Chain complexes will almost never be stable when
        they have objects in distinct shifts. The phase of an individual element of a chain complex generally
        requires that object to be stable, so that we typically use LineBundles for our current applications.
        By definition of a slicing, the shift of each object in the complex should add to the respective phases;
        thus, this method computes the largest sum of the phase of the sheaf and the shift of the sheaf in the
        complex.

        \param tuple args The arguments required to compute the phase. The number of arguments and the type of arguments will
                     depend on the catagory of the sheaves in the complex. For P1, the phase requires a single complex number.
                     For P2, the phase requires two floating-point numbers. For K3, the phase requires two floating-point
                     numbers and an integer.

        \return float The largest phase of the chain complex

        \throws ValueError If the number of arguments is incorrect
        \throws TypeError If the type of the arguments is incorrect
        \throws NotImplementedError If the catagory of the sheaves in the complex is not implemented
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Phase of P1 requires exactly one argument.")
            if not isinstance(args[0], complex):
                raise TypeError("Phase of P1 requires a complex number as an argument.")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Phase of P2 requires exactly two arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Phase of P2 requires two floating-point numbers as arguments.")
        elif self.catagory == 'K3':
            if len(args) != 3:
                raise ValueError("Phase of K3 requires exactly three arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Phase of K3 requires three floating-point numbers as arguments.")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer.")
        else:
            raise NotImplementedError("Only local P1, local P2, and K3 catagories are implemented.")

        # Zip the three lists together and sort by descending shift
        max_shift = max(self.shift_vector)

        bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
        bundles_max_shift = filter(lambda x: x[2] == max_shift or x[2] == max_shift - 1, bundles)

        max_phase = -math.inf

        for sheaf, dim, shift in bundles_max_shift:
            if sheaf.phase(*args) + shift > max_phase:
                max_phase = sheaf.phase(*args) + shift

        return max_phase
    



    def is_semistable(self, *args):
        r"""!
        Method to compute whether the chain complex is semistable. This almost never occurs, since if
        the complex contains two or more stable objects of distinct phase, it will never be stable. For
        example, suppose E2 is a stable subobject of maximum phase and E1 is another stable object with
        strictly smaller phase. Then

                          E2 ----> Complex ------> Cone

        will destabilize the complex, and Cone will be nontrivial since it has a non-zero map to E1. 
        The easiest way to check that the complex is concentrated in only a single phase is to compare
        its largest and smallest phases from the previous methods.


        \param tuple args The arguments required to compute the phase. The number of arguments and the type of arguments will
                     depend on the catagory of the sheaves in the complex. For P1, the phase requires a single complex number.
                     For P2, the phase requires two floating-point numbers. For K3, the phase requires two floating-point
                     numbers and an integer.

        \return bool True if the chain complex is semistable, False otherwise

        \throws ValueError If the number of arguments is incorrect
        \throws TypeError If the type of the arguments is incorrect
        \throws NotImplementedError If the catagory of the sheaves in the complex is not implemented
        """

    
        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Phase of P1 requires exactly one argument.")
            if not isinstance(args[0], complex):
                raise TypeError("Phase of P1 requires a complex number as an argument.")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Phase of P2 requires exactly two arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Phase of P2 requires two floating-point numbers as arguments.")
        elif self.catagory == 'K3':
            if len(args) != 3:
                raise ValueError("Phase of K3 requires exactly three arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Phase of K3 requires three floating-point numbers as arguments.")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer.")
        else:
            raise NotImplementedError("Only local P1, local P2, and K3 catagories are implemented.")

        return self.get_largest_phase(*args) == self.get_smallest_phase(*args)



            