from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject
from typing import List

from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']



class GradedCoproductObject(DerivedCategoryObject):

    def __init__(self, object_vector : List[DerivedCategoryObject], shift_vector : List[int], dimension_vector : List[int] = None):
        r"""!
        Initialize an instance of GradedCoproductObject with the specified object_vector, shift vector,
        and potentially a dimension vector. If a dimension vector is provided, it must consist of non-negative;
        if it is not provided, it is assumed all dimensions are 1.


        \param list object_vector A list of DerivedCategoryObjects 
        \param list shift_vector A list of homological shifts in the complex
        \param list dimension_vector A list of the number of direct sums of each object in the sum


        \throws ValueError If the object vector is empty
        \throws ValueError If the object vector and shift vector have different lengths
        \throws TypeError If any element of the object vector is not a DerivedCategoryObject
        \throws TypeError If any element of the shift vector is not an integer
        \throws ValueError If the dimension vector is not the same length as the object vector
        \throws TypeError If any element of the dimension vector is not an integer

        \throws ValueError If any element of the dimension vector is negative
        \throws ValueError If the catagory of the objects in the complex is not implemented
        \throws ValueError If the object vector contains objects of different catagories
        """


        ##### 
        # Check that object vector is valid and contains object of only a single catagory
        #####

        if not object_vector:
            raise ValueError("object_vector cannot be empty.")

        if not all(isinstance(obj, DerivedCategoryObject) for obj in object_vector) :
            raise TypeError("All elements of complex_vector must be instances of DerivedCategoryObject (of the same projective space).")

        object_catagory = object_vector[0].catagory
        if object_catagory not in IMPLEMENTED_CATAGORIES:
            raise ValueError(f"Catagory {object_catagory} is not implemented.")

        if not all(obj.catagory == object_catagory for obj in object_vector):
            raise ValueError("All elements of object_vector must be sheaves over the same base space.")

        #####
        # Check that shift vector is valid
        #####
        if not all(isinstance(shift, int) for shift in shift_vector):
            raise TypeError("All elements of shift_vector must be integers.")

        if len(object_vector) != len(shift_vector):
            raise ValueError("object_vector and shift_vector must have the same length.") 

        #####
        # Check that dimension vector is valid
        #####
        if dimension_vector is None:
            dimension_vector = [1] * len(object_vector)
        if not all(isinstance(dim, int) for dim in dimension_vector):
            raise TypeError("All elements of dimension_vector must be integers.")
        if len(object_vector) != len(dimension_vector):
            raise ValueError("object_vector and dimension_vector must have the same length.")
        if not all(dim >= 0 for dim in dimension_vector):
            # Dimension cannot be non-negative
            raise ValueError("All elements of dimension_vector must be non-negative integers.") 
        
        self.object_vector = object_vector ## List of coherent sheaves in the complex, so that the chain complex can operate similar to a DenseVector.

        self.dimension_vector = dimension_vector ## List of the number of direct sums of each object in the complex.

        self.shift_vector = shift_vector ## List of homological shifts in the complex.

        self.catagory = object_catagory ## The catagory of the sheaves in the complex.

        # If an element of the complex has dimension 0, we can get rid of it using helper method
        self._remove_zeros_from_dimension_vector()
        # Combine repeated sheaves in the complex
        self._combine_repeats()



    def __str__(self):
        r"""!
        String representation of the coproduct object. The complex is represented in cohomological order 
        (which technically would be descending order of the shifts, since IR[-2] means the complex with
        a copy of IR in index 2). The individual coherent sheaves in the complex are represented by their
        own respective print functinos --- this will generally be cumbersome for arbitrary vector bundles,
        but more clean for named instances like O(1) or Ω(3). 

        \return str A string representation of the chain complex
        """
         # Zip the three lists together and sort by descending shift
        bundles = list(zip(self.object_vector, self.dimension_vector, self.shift_vector))
        bundles_sorted = sorted(bundles, key=lambda x: x[2], reverse=True)
        
        # Define the arrow string used to join bottom elements
        arrow = " --------- "
        # We'll use the same spacing (without arrow characters) to join top columns
        join_space = " " * len(arrow)
        
        top_columns = []
        bottom_columns = []
        for object, dim, shift in bundles_sorted:
            # Create the bottom string: e.g., "O(5)[3]"
            bottom_elem = f"{object}[{shift}]"
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
        return len(self.object_vector)
    
    def chernCharacter(self):
        r"""!
        Helper function to compute the Chern Character of the chain complex. The Chern Character of
        a chain complex is the alternating sum of the Chern Characters of the individual sheaves in
        the complex. Since the Chern character is additive, we may multiply the Chern Characters by
        the dimension of the object to represent direct sums of sheaves.

        \return ChernCharacter The Chern Character of the chain complex
        """
        cherns = [object.chernCharacter() for object in self.object_vector]

        # We store the first Chern Character instead of making a "Zero" Chern Character and iterating
        # normally since there is currently no functionality to make a "Zero" chern character (the 
        # dimension attribute makes it difficult to do so).
        chern_to_return = (-1)**(self.shift_vector[0]) * self.dimension_vector[0] * cherns[0]

        for i in range(1, len(cherns)):
            # odd shifts get a negative sign, even shifts get a positive sign
            chern_piece += (-1)**(self.shift_vector[i]) * self.dimension_vector[i] * cherns[i]

        return chern_to_return
    
    
    def shift(self, shift : int):
        r"""!
        Method to shift the chain complex by a given homological shift

        \param int shift The homological shift to apply to the chain complex

        \return GradedCoproductObject shifted by the homological shift
        """

        new_shift_vector = [shift + s for s in self.shift_vector]

        return GradedCoproductObject(self.object_vector, new_shift_vector, self.dimension_vector)
    
    def _remove_zeros_from_dimension_vector(self):
        r"""!
        Helper function which iterates through the dimension vector, and if a certain object
        is only included 0 times, we may effectively erase it.
        """
        for i in range(len(self.dimension_vector)):
            if i < len(self.dimension_vector) and self.dimension_vector[i] == 0:
                del self.object_vector[i]
                del self.shift_vector[i]
                del self.dimension_vector[i]

    def _combine_repeats(self):
        r"""!
        Helper function to combine repeated sheaves in the complex. This is useful for simplifying
        the complex, as we can combine repeated sheaves into a single object with a larger dimension.
        This function specifically requires the __hash__ implementation for the DerivedCategoryObject
        """

        # Dictionary to hold combined dimensions for each (object, shift) pair
        combined = {}

        # Iterate over the tuples
        for object, dim, shift in zip(self.object_vector, self.dimension_vector, self.shift_vector):
            key = (object, shift)
            if key in combined:
                combined[key] += dim
            else:
                combined[key] = dim

        # Unpack the combined dictionary back into the class variables
        self.object_vector = []
        self.dimension_vector = []
        self.shift_vector = []

        for (object, shift), dim in combined.items():
            self.object_vector.append(object)
            self.dimension_vector.append(dim)
            self.shift_vector.append(shift)
    
