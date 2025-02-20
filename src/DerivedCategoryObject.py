from .CoherentSheaf import CoherentSheaf
from .ChernCharacter import ChernCharacter
import re
import math



###############################################################################
#                                                                             #
#                             Multi-Degree Objects                            #
# ----------------------------------------------------------------------------#
#  These objects are used to represent chain complexes over a the category of #              
#  coherent sheaves; these objects are stored primarily as an array of degree #
#  , ranks, and homological shifts. Despite the fact that ChainComplexes are  #
#  partially implemented, LineBundleComplex objects are much more useful for  #
#  several methods in this project (note in particular that every coherent    #
#  sheaf on P^2 has a resolution by line bundles). Consequently, it is        #
#  important to make sure every function in ChainComplex is overridden in     #
#  LineBundleComplex.                                                         #
#                                                                             #
#  The Distinguished triangle further generalizes the notion of ordered       #                                                     #
#  exact sequences in the Derived category. As our main object of interest,   #
#  SphericalTwistComposition, exists in the Derived category, it will be      #
#  crucial to understand the passage from complexes of coherent sheaves to    #
#  distinguished triangles.                                                   #
#                                                                             #
###############################################################################



class DerivedCategoryObject():
    """
    General parent class for both ChainComplex and LineBundleComplex objects. This class
    is used to represent objects in the derived category of coherent sheaves on P^2. Since
    objects in the derived category are technically equivalence classes of complexes of
    coherent sheaves (up to quasi-isomorphism), we will use an abstract parent class to 
    capture only string information and numerical information.

    Attributes:
    ----------
    string : str
        A string representation of the derived category
    chern_character : ChernCharacter
        The Chern Character of the derived category object
    """


    def __init__(self, string = "0", chern_character = ChernCharacter(0,0,0)):
        """
        Initialize an instance of DerivedCategoryObject with the specified string and Chern Character.

        Parameters:
        ----------
        string : str
            The string representation of the derived category object
        chern_character : ChernCharacter
            The Chern Character of the derived category object
        """
        self.string = string
        self.chern_character = chern_character
    
    def __str__(self):
        """
        String representation of the derived category object

        Returns:
        -------
        str
            A string representation of the derived category
        """
        if self.string is not None:
            return self.string
        else:
            return "0"


    def chernCharacter(self):
        """
        Method to return the Chern Character of the derived category object

        Returns:
        -------
        ChernCharacter
            The Chern Character of the derived category object
        """
        return self.chern_character
    

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
        
        if self.chern_character == ChernCharacter(0,0,0):
            raise ValueError("DerivedCategoryObject not initialized, cannot compute central charge.")
    
        chern_char = self.chernCharacter()
        return complex(-chern_char.ch2 + q * chern_char.ch0, chern_char.ch1 - s * chern_char.ch0)


    def shiftComplex(self, shift):
        """
        Method to shift the derived category object by a given homological shift

        Parameters:
        ----------
        shift : int
            The homological shift to apply to the derived category object

        Returns:
        -------
        DerivedCategoryObject
            The derived category object shifted by the homological shift
        """
        new_string =  self._update_string_by_shift(self.string, shift)
        new_chern = int((-1)**shift) * self.chern_character 
        return DerivedCategoryObject(new_string, new_chern)


    def _update_string_by_shift(self, my_str, n):
        """
        Static helper function to update the possible string representations of the abstract
        DerivedCategoryObject by a homological shift. For example, if an object is called A[3],
        then a shift of 2 will yield A[5] (as a string). If the object is unshifted to start with,
        then shifting the object should return A[n] (as a string). The only object that should not
        be shifted is the zero object, which is represented by the string "0".

        Parameters:
        ----------
        my_str : str
            The string representation of the derived category object
        n : int
            The homological shift to apply to the derived category object

        Returns:
        -------
        str
            The string representation of the derived category object shifted by the homological shift

        Raises:
        -------
        TypeError
            If my_str is not a string
            If n is not an integer

        """

        if not isinstance(my_str, str):
            raise TypeError("my_str must be a string.")
        if not isinstance(n, int):
            raise TypeError("n must be an integer.")

        # If the object is zero, it should remain zero
        if my_str == "0":
            return "0"
        
        # This regex checks for a pattern at the end of the string that looks like "[number]"
        pattern = r'\[(\d+)\]$'
        match = re.search(pattern, my_str)
        
        if match:
            # Extract the current number k, add n, and format the new number.
            k = int(match.group(1))
            new_k = k + n
            # Replace the [k] at the end with [new_k]
            updated_str = re.sub(pattern, f'[{new_k}]', my_str)
            return updated_str
        else:
            # If there's no [number] at the end, append [n]
            return my_str + f'[{n}]'







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

        if not all(isinstance(obj, CoherentSheaf) for obj in sheaf_vector):
            raise TypeError("All elements of complex_vector must be instances of CoherentSheaf.")

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

        return ChernCharacter(ch0, ch1, ch2)


    
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

        if not isinstance(derived_object1, DerivedCategoryObject):
            raise TypeError("derived_object1 must be an instance of DerivedCategoryObject.")
        if not isinstance(derived_object2, DerivedCategoryObject):
            raise TypeError("derived_object2 must be an instance of DerivedCategoryObject.")
        if not isinstance(derived_object3, DerivedCategoryObject):
            raise TypeError("derived_object3 must be an instance of DerivedCategoryObject.")

        self.object1 = derived_object1
        self.object2 = derived_object2
        self.object3 = derived_object3

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

        if self.object1.chernCharacter() == ChernCharacter(0,0,0) and not isinstance(self.object1, ChainComplex):
            # Update the Chern character of the first complex

            chern2 = self.object2.chernCharacter()
            chern3 = self.object3.chernCharacter()

            self.object1.chern_character = chern2 - chern3
        elif self.object2.chernCharacter() == ChernCharacter(0,0,0) and not isinstance(self.object2, ChainComplex):
            # Update the Chern character of the second complex

            chern1 = self.object1.chernCharacter()
            chern3 = self.object3.chernCharacter()

            self.object2.chern_character = chern1 + chern3
        elif self.object3.chernCharacter() == ChernCharacter(0,0,0) and not isinstance(self.object3, ChainComplex):
            # Update the Chern character of the third complex

            chern1 = self.object1.chernCharacter()
            chern2 = self.object2.chernCharacter()
            
            self.object3.chern_character = chern2 - chern1