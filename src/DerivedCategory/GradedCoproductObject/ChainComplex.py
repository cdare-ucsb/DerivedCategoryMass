from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject
from src.DerivedCategory.CoherentSheaf.CoherentSheaf import CoherentSheaf
from src.DerivedCategory.ChernCharacter import ChernCharacter
import math

from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']


class ChainComplex(DerivedCategoryObject):
    """!
    For any abelian category, a chain complex is a sequence of objects and morphisms between them
    such that the composition of any two consecutive morphisms is the zero morphism. In the derived
    category of coherent sheaves on P^2, we can represent a chain complex as a sequence of coherent
    sheaves with a shift. For instance, a general complex will be of the form

                 i=-n       i=-n+1      i=-n+2    ...
        0 ------> E1 -------> E2 --------> E3 ---> ...

    (A priori, there is no reason the complexes cant also descend infinitely in the other direction). 
    For the purposes of this project, only finite complexes will be considered. Such a complex can be
    stored in a similar way to a DenseVector object --- namely, since the majority of entries in the
    complex will be zero, we can store the complex as a list of coherent sheaves and a shift vector.
    
    """

    def __init__(self, sheaf_vector, shift_vector, dimension_vector = None):
        r"""!
        Initialize an instance of ChainComplex with the specified sheaf vector, shift vector,
        and potentially a dimension vector. If a dimension vector is not provided, it must 
        consist of non-negative integer values


        \param list sheaf_vector A list of coherent sheaves in the complex
        \param list shift_vector A list of homological shifts in the complex
        \param list dimension_vector A list of the number of direct sums of each coherent sheaf in the complex


        \throws ValueError If the sheaf vector is empty
        \throws ValueError If the sheaf vector and shift vector have different lengths
        \throws TypeError If any element of the sheaf vector is not a CoherentSheaf object
        \throws TypeError If any element of the shift vector is not an integer
        \throws ValueError If the dimension vector is not the same length as the sheaf vector
        \throws TypeError If any element of the dimension vector is not an integer
        \throws ValueError If any element of the dimension vector is negative
        \throws ValueError If the catagory of the sheaves in the complex is not implemented
        \throws ValueError If the sheaf vector contains objects of different catagories
        \throws ValueError If the sheaf vector contains objects of different base spaces
        \throws ValueError If the sheaf vector contains objects of different projective spaces
        \throws TypeError If the sheaf vector contains objects of different projective spaces
        """


        ##### 
        # Check that sheaf vector is valid and contains object of only a single catagory
        #####

        if not sheaf_vector:
            raise ValueError("sheaf_vector cannot be empty.")

        if not all(isinstance(obj, CoherentSheaf) for obj in sheaf_vector) :
            raise TypeError("All elements of complex_vector must be instances of CoherentSheaf (of the same projective space).")

        sheaf_catagory = sheaf_vector[0].catagory
        if sheaf_catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {sheaf_catagory} is not implemented.")

        if not all(obj.catagory == sheaf_catagory for obj in sheaf_vector):
            raise ValueError("All elements of sheaf_vector must be sheaves over the same base space.")

        #####
        # Check that shift vector is valid
        #####
        if not all(isinstance(shift, int) for shift in shift_vector):
            raise TypeError("All elements of shift_vector must be integers.")

        if len(sheaf_vector) != len(shift_vector):
            raise ValueError("sheaf_vector and shift_vector must have the same length.") 

        #####
        # Check that dimension vector is valid
        #####
        if dimension_vector is None:
            dimension_vector = [1] * len(sheaf_vector)
        if not all(isinstance(dim, int) for dim in dimension_vector):
            raise TypeError("All elements of dimension_vector must be integers.")
        if len(sheaf_vector) != len(dimension_vector):
            raise ValueError("sheaf_vector and dimension_vector must have the same length.")
        if not all(dim >= 0 for dim in dimension_vector):
            # Dimension cannot be non-negative
            raise ValueError("All elements of dimension_vector must be non-negative integers.") 
        
        self.sheaf_vector = sheaf_vector ## List of coherent sheaves in the complex, so that the chain complex can operate similar to a DenseVector.

        self.dimension_vector = dimension_vector ## List of the number of direct sums of each coherent sheaf in the complex.

        self.shift_vector = shift_vector ## List of homological shifts in the complex.

        self.catagory = sheaf_catagory ## The catagory of the sheaves in the complex.

        # If an element of the complex has dimension 0, we can get rid of it using helper method
        self._remove_zeros_from_dimension_vector()
        # Combine repeated sheaves in the complex
        self._combine_repeats()

    
    

        

    def __str__(self):
        r"""!
        String representation of the chain complex. The complex is represented in cohomological order 
        (which technically would be descending order of the shifts, since IR[-2] means the complex with
        a copy of IR in index 2). The individual coherent sheaves in the complex are represented by their
        own respective print functinos --- this will generally be cumbersome for arbitrary vector bundles,
        but more clean for named instances like O(1) or Ω(3). 

        \return str A string representation of the chain complex
        """
         # Zip the three lists together and sort by descending shift
        bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
        bundles_sorted = sorted(bundles, key=lambda x: x[2], reverse=True)
        
        # Define the arrow string used to join bottom elements
        arrow = " --------> "
        # We'll use the same spacing (without arrow characters) to join top columns
        join_space = " " * len(arrow)
        
        top_columns = []
        bottom_columns = []
        for sheaf, dim, shift in bundles_sorted:
            # Create the bottom string: e.g., "O(5)[3]"
            bottom_elem = f"{sheaf}[{shift}]"
            col_width = len(bottom_elem)
            # Create the top string: e.g., "⊕7" and right-align it to match the bottom column width
            top_elem = f"⊕{dim}"
            top_columns.append(top_elem.rjust(col_width))
            bottom_columns.append(bottom_elem)
        
        top_line = join_space.join(top_columns)
        bottom_line = arrow.join(bottom_columns)
        return top_line + "\n" + bottom_line

    def __len__(self):
        r"""!
        The length of the chain complex is the number of sheaves in the complex

        \return int The number of sheaves in the complex
        """
        return len(self.sheaf_vector)
        
    

    def chernCharacter(self):
        r"""!
        Helper function to compute the Chern Character of the chain complex. The Chern Character of
        a chain complex is the alternating sum of the Chern Characters of the individual sheaves in
        the complex. Since the Chern character is additive, we may multiply the Chern Characters by
        the dimension of the sheaf to represent direct sums of sheaves.

        \return ChernCharacter The Chern Character of the chain complex
        """
        cherns = [sheaf.chernCharacter() for sheaf in self.sheaf_vector]

        rank_grothendieck_group = len(self.sheaf_vector[0].chernCharacter())
        chain_complex_chern = []

        for graded_piece in range(rank_grothendieck_group):

            chern_piece = 0
            for i in range(len(cherns)):
                # odd shifts get a negative sign, even shifts get a positive sign
                chern_piece += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i][graded_piece]
            chain_complex_chern.append(chern_piece)

        return ChernCharacter(chain_complex_chern)

    def central_charge(self, *args):
        r"""!
        Compute the central charge of the chain complex. The central charge of a chain complex is the
        alternating sum of the central charges of the individual sheaves in the complex. Since the
        central charge is additive, we may multiply the central charges by the dimension of the sheaf
        to represent direct sums of sheaves. However, most of this functionality is already defined
        in the chernCharacter() function, so we will simply call that function and then compute the
        central charge from the Chern Character.

        \param tuple args The arguments required to compute the central charge. The number of arguments and the type
                     of arguments will depend on the catagory of the sheaves in the complex. For P1, the central
                     charge requires a single complex number. For P2, the central charge requires two floating-point
                     numbers. For K3, the central charge requires two floating-point numbers and an integer.

        \return complex The central charge of the chain complex as a complex number

        \throws ValueError If the number of arguments is incorrect
        \throws TypeError If the type of the arguments is incorrect
        \throws NotImplementedError If the catagory of the sheaves in the complex is not implemented
        """

        
        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge for P1 requires exactly one argument.")
            if not isinstance(args[0], complex):
                raise TypeError("Central charge for P1 requires a complex number as an argument.")

            ch = self.chernCharacter()
            return complex(-1*ch[1] + args[0] * ch[0])


        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge for P2 requires exactly two arguments.")
            if not all(isinstance(arg, (float,int)) for arg in args):
                raise TypeError("Central charge for P2 requires two floating-point numbers as arguments.")
            
            ch = self.chernCharacter()
            return complex(-1*ch[2] + args[1] * ch[0],
                            ch[1] - args[0] * ch[0])
        
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
            ch = self.chernCharacter()
            
            return complex(2*d*alpha * ch[1] - ch[2] - ch[0] + (beta**2 - alpha**2)*d*ch[0], 
                           2*d*ch[1] - 2*d*alpha*beta*ch[0])

        else:
            raise NotImplementedError("Central charge not implemented for this variety.")


    
    def shift(self, shift):
        r"""!
        Method to shift the chain complex by a given homological shift

        \param int shift The homological shift to apply to the chain complex

        \return ChainComplex The chain complex shifted by the homological shift
        """

        new_shift_vector = [shift + s for s in self.shift_vector]

        return ChainComplex(self.sheaf_vector, new_shift_vector, self.dimension_vector)

    def isShiftOfSheaf(self):
        r"""!
        Simple helper function which checks if the complex is a shift of a single sheaf

        \return bool True if the complex is a shift of a single sheaf, False otherwise
        """

        return len(self.sheaf_vector) == 1


    def _remove_zeros_from_dimension_vector(self):
        r"""!
        Helper function which iterates through the dimension vector, and if a certain Coherent sheaf
        is only included 0 times, we may effectively erase it.
        """
        for i in range(len(self.dimension_vector)):
            if i < len(self.dimension_vector) and self.dimension_vector[i] == 0:
                del self.sheaf_vector[i]
                del self.shift_vector[i]
                del self.dimension_vector[i]

    def _combine_repeats(self):
        r"""!
        Helper function to combine repeated sheaves in the complex. This is useful for simplifying
        the complex, as we can combine repeated sheaves into a single sheaf with a larger dimension.
        This function specifically requires the __hash__ implementation for the CoherentSheaf and 
        LineBundle objects.
        """

        # Dictionary to hold combined dimensions for each (sheaf, shift) pair
        combined = {}

        # Iterate over the tuples
        for sheaf, dim, shift in zip(self.sheaf_vector, self.dimension_vector, self.shift_vector):
            key = (sheaf, shift)
            if key in combined:
                combined[key] += dim
            else:
                combined[key] = dim

        # Unpack the combined dictionary back into the class variables
        self.sheaf_vector = []
        self.dimension_vector = []
        self.shift_vector = []

        for (sheaf, shift), dim in combined.items():
            self.sheaf_vector.append(sheaf)
            self.dimension_vector.append(dim)
            self.shift_vector.append(shift)

        
    
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



            