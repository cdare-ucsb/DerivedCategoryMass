from .DerivedCategoryObject import DerivedCategoryObject, DerivedCategoryObjectP1, DerivedCategoryObjectP2
from .CoherentSheaf import CoherentSheafP1, CoherentSheafP2
from .ChernCharacter import ChernCharacterP2, ChernCharacterP1, ChernCharacter
import math


class ChainComplex(DerivedCategoryObject):
    """
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

    Attributes:
    ----------
    sheaf_vector : list
        A list of coherent sheaves in the complex
    shift_vector : list
        A list of homological shifts in the complex
    dimension_vector : list
        A list of the number of direct sums of each coherent sheaf in the complex
    
    """

    def __init__(self, sheaf_vector, shift_vector, dimension_vector = None):
        """
        Initialize an instance of ChainComplex with the specified sheaf vector, shift vector,
        and potentially a dimension vector. If a dimension vector is not provided, it must 
        consist of non-negative integer values

        Parameters:
        ----------
        sheaf_vector : list
            A list of coherent sheaves in the complex
        shift_vector : list
            A list of homological shifts in the complex
        dimension_vector : list
            A list of the number of direct sums of each coherent sheaf in the complex. By default,
            it will simply be set to be [1, 1, ..., 1] so that there is only one copy of each sheaf.

        Raises:
        -------
        ValueError
            If the sheaf vector is empty
            If the sheaf vector and shift vector have different lengths
        TypeError
            If any element of the sheaf vector is not a CoherentSheaf object
            If any element of the shift vector is not an integer
        """

        if not sheaf_vector:
            raise ValueError("sheaf_vector cannot be empty.")

        if not all(isinstance(obj, CoherentSheafP1) for obj in sheaf_vector) and not all(isinstance(obj, CoherentSheafP2) for obj in sheaf_vector):
            raise TypeError("All elements of complex_vector must be instances of CoherentSheaf (of the same projective space).")

        if not all(isinstance(shift, int) for shift in shift_vector):
            raise TypeError("All elements of shift_vector must be integers.")

        if len(sheaf_vector) != len(shift_vector):
            raise ValueError("sheaf_vector and shift_vector must have the same length.") 

        if dimension_vector is None:
            dimension_vector = [1] * len(sheaf_vector)
        elif not all(isinstance(dim, int) for dim in dimension_vector):
            raise TypeError("All elements of dimension_vector must be integers.")
        elif len(sheaf_vector) != len(dimension_vector):
            raise ValueError("sheaf_vector and dimension_vector must have the same length.")
        # Dimension cannot be non-negative
        elif not all(dim >= 0 for dim in dimension_vector):
            raise ValueError("All elements of dimension_vector must be non-negative integers.") 
        
        self.sheaf_vector = sheaf_vector
        self.dimension_vector = dimension_vector
        self.shift_vector = shift_vector

        # If an element of the complex has dimension 0, we can get rid of it using helper method
        self._remove_zeros_from_dimension_vector()

        

    def __str__(self):
        """
        String representation of the chain complex. The complex is represented in cohomological order 
        (which technically would be descending order of the shifts, since IR[-2] means the complex with
        a copy of IR in index 2). The individual coherent sheaves in the complex are represented by their
        own respective print functinos --- this will generally be cumbersome for arbitrary vector bundles,
        but more clean for named instances like O(1) or Ω(3). 

        Example
        -------
        For a complex with O(3)

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

        
    

    def chernCharacter(self):
        return ChernCharacter()


    
    def shiftComplex(self, shift):
        """
        Method to shift the chain complex by a given homological shift

        Parameters:
        ----------
        shift : int
            The homological shift to apply to the chain complex

        Returns:
        -------
        ChainComplex
            The chain complex shifted by the homological shift
        """
        new_shift_vector = [shift + s for s in self.shift_vector]

        return ChainComplex(self.sheaf_vector, new_shift_vector, self.dimension_vector)

    def isShiftOfSheaf(self):
        """
        Simple helper function which checks if the complex is a shift of a single sheaf

        Returns:
        -------
        bool
            True if the complex is a shift of a single sheaf, False otherwise
        """
        return len(self.complex) == 1


    def _remove_zeros_from_dimension_vector(self):
        """
        Helper function which iterates through the dimension vector, and if a certain Coherent sheaf
        is only included 0 times, we may effectively erase it.
        """
        for i in range(len(self.dimension_vector)):
            if self.dimension_vector[i] == 0:
                del self.sheaf_vector[i]
                del self.shift_vector[i]
                del self.dimension_vector[i]

    

class ChainComplexP1(ChainComplex, DerivedCategoryObjectP1):

    def __init__(self, sheaf_vector, shift_vector, dimension_vector = None):    
        if not all(isinstance(obj, CoherentSheafP1) for obj in sheaf_vector):
            raise TypeError("All elements of complex_vector must be instances of CoherentSheafP1.")

        super().__init__(sheaf_vector, shift_vector, dimension_vector)

    def __str__(self):
        return super().__str__()
    
    def chernCharacter(self):
        """
        Helper function to compute the Chern Character of the chain complex. The Chern Character of
        a chain complex is the alternating sum of the Chern Characters of the individual sheaves in
        the complex. Since the Chern character is additive, we may multiply the Chern Characters by
        the dimension of the sheaf to represent direct sums of sheaves.
        """
        cherns = [sheaf.chernCharacter() for sheaf in self.sheaf_vector]

        ch0 = 0
        ch1 = 0

        for i in range(len(cherns)):
            ch0 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch0
            ch1 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch1

        return ChernCharacterP1(ch0, ch1)
    
    def central_charge(self, w):
        """
        Method to compute the central charge of the chain complex. The central charge of a chain complex
        is the alternating sum of the central charges of the individual sheaves in the complex. Since the
        central charge is additive, we may multiply the central charges by the dimension of the sheaf to
        represent direct sums of sheaves.

        Parameters:
        ----------
        w : complex
            The complex parameter for a geometric stability condition in P1
        

        Returns:
        -------
        complex
            The central charge of the chain complex as a complex number

        Raises:
        -------
        TypeError
            If w is not a complex number

        """
        
        if not isinstance(w, complex):
            raise TypeError("w must be a complex number.")
    
        chern_char = self.chernCharacter()
        return complex(-chern_char.ch1 + w * chern_char.ch0)
    
    def get_smallest_phase(self, w):

         # Zip the three lists together and sort by descending shift
        min_shift = min(self.shift_vector)

        bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
        bundles_min_shift = filter(lambda x: x[2] == min_shift or x[2] == min_shift + 1 , bundles)

        min_phase = math.inf

        for sheaf, dim, shift in bundles_min_shift:
            if sheaf.phase(w) + shift < min_phase:
                min_phase = sheaf.phase(w) + shift

        return min_phase
    
    def get_largest_phase(self, w):

            # Zip the three lists together and sort by descending shift
            max_shift = max(self.shift_vector)
    
            bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
            bundles_max_shift = filter(lambda x: x[2] == max_shift or x[2] == max_shift - 1, bundles)
    
            max_phase = -math.inf
    
            for sheaf, dim, shift in bundles_max_shift:
                if sheaf.phase(w) + shift > max_phase:
                    max_phase = sheaf.phase(w) + shift
    
            return max_phase
    
    def is_semistable(self, w):

        return self.get_largest_phase(w) == self.get_smallest_phase(w)



            

class ChainComplexP2(ChainComplex, DerivedCategoryObjectP2):

    def __init__(self, sheaf_vector, shift_vector, dimension_vector = None):  
        if not all(isinstance(obj, CoherentSheafP2) for obj in sheaf_vector):
            raise TypeError("All elements of complex_vector must be instances of CoherentSheafP2.")  
        super().__init__(sheaf_vector, shift_vector, dimension_vector)

    def __str__(self):
        return super().__str__()
    
    def chernCharacter(self):
        """
        Helper function to compute the Chern Character of the chain complex. The Chern Character of
        a chain complex is the alternating sum of the Chern Characters of the individual sheaves in
        the complex. Since the Chern character is additive, we may multiply the Chern Characters by
        the dimension of the sheaf to represent direct sums of sheaves.
        """
        cherns = [sheaf.chernCharacter() for sheaf in self.sheaf_vector]

        ch0 = 0
        ch1 = 0
        ch2 = 0

        for i in range(len(cherns)):
            ch0 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch0
            ch1 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch1
            ch2 += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i].ch2

        return ChernCharacterP2(ch0, ch1, ch2)
    


    
    def central_charge(self, s, q):
        """
        Method to compute the central charge of the chain complex. The central charge of a chain complex
        is the alternating sum of the central charges of the individual sheaves in the complex. Since the
        central charge is additive, we may multiply the central charges by the dimension of the sheaf to
        represent direct sums of sheaves.

        Parameters:
        ----------
        s : float
            The parameter controlling the imaginary part of the central charge
        q : float
            The parameter controlling the real part of the central charge

        Returns:
        -------
        complex
            The central charge of the chain complex as a complex number

        Raises:
        -------
        TypeError
            If s or q are not floating-point decimals

        """
        
        if not isinstance(s, float) or not isinstance(q, float):
            raise TypeError("s and q must be floating-point decimals.")
    
        chern_char = self.chernCharacter()
        return complex(-chern_char.ch2 + q * chern_char.ch0, chern_char.ch1 - s * chern_char.ch0)
                
    def get_smallest_phase(self, s, q):

         # Zip the three lists together and sort by descending shift
        min_shift = min(self.shift_vector)

        bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
        bundles_min_shift = filter(lambda x: x[2] == min_shift or x[2] == min_shift + 1 , bundles)

        min_phase = math.inf

        for sheaf, dim, shift in bundles_min_shift:
            if sheaf.phase(s, q) + shift < min_phase:
                min_phase = sheaf.phase(s, q) + shift

        return min_phase
    
    def get_largest_phase(self, s, q):

            # Zip the three lists together and sort by descending shift
            max_shift = max(self.shift_vector)
    
            bundles = list(zip(self.sheaf_vector, self.dimension_vector, self.shift_vector))
            bundles_max_shift = filter(lambda x: x[2] == max_shift or x[2] == max_shift - 1, bundles)
    
            max_phase = -math.inf
    
            for sheaf, dim, shift in bundles_max_shift:
                if sheaf.phase(s, q) + shift > max_phase:
                    max_phase = sheaf.phase(s, q) + shift
    
            return max_phase
    
    def is_semistable(self, s, q):

        return self.get_largest_phase(s,q) == self.get_smallest_phase(s,q)



