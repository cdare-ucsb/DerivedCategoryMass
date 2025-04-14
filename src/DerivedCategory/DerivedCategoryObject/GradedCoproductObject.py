from src.DerivedCategory.DerivedCategoryObject.DerivedCategoryObject import DerivedCategoryObject

from typing import List
from functools import reduce
from operator import add

from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']



class GradedCoproductObject(DerivedCategoryObject):

    _instances = {} ## Memoization: Dictionary to hold instances of GradedCoproductObject


    def __new__(cls, object_vector : List[DerivedCategoryObject], shift_vector : List[int] = None, dimension_vector : List[int] = None):
        
        shift_key = tuple(shift_vector) if shift_vector is not None else None
        dim_key = tuple(dimension_vector) if dimension_vector is not None else None
            

        key = (
            tuple(object_vector),
            shift_key,
            dim_key
        )

        if key not in cls._instances:
            instance = super().__new__(cls)
            cls._instances[key] = instance

        return cls._instances[key]
    


    def __init__(self, object_vector : List[DerivedCategoryObject], shift_vector : List[int] = None, dimension_vector : List[int] = None ):
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

        if hasattr(self, '_initialized') and self._initialized:
            return


        ##### 
        # Check that object vector is valid and contains object of only a single catagory
        #####

        if not object_vector:
            raise ValueError("object_vector cannot be empty.")

        if not all(isinstance(obj, DerivedCategoryObject) for obj in object_vector) :
            raise TypeError("All elements of complex_vector must be instances of DerivedCategoryObject (of the same projective space).")

        first_geometry_context = object_vector[0].geometry_context
        if not all(obj.geometry_context == first_geometry_context for obj in object_vector):
            raise ValueError("All elements of object_vector must have the same geometry context.")
        
        self.geometry_context = first_geometry_context ## Geometry context of the sheaves in the complex, including the catagory and the intersection form

        #####
        # Check that shift vector is valid
        #####

        if shift_vector and not all(isinstance(shift, int) for shift in shift_vector):
            raise TypeError("All elements of shift_vector must be integers.")

        if shift_vector and len(object_vector) != len(shift_vector):
            raise ValueError("object_vector and shift_vector must have the same length.") 
        
        if shift_vector is None:
            shift_vector = [0] * len(object_vector)
        

        #####
        # Check that dimension vector is valid
        #####

        if dimension_vector and not all(isinstance(dim, int) for dim in dimension_vector):
            raise TypeError("All elements of dimension_vector must be integers.")
        if dimension_vector and len(object_vector) != len(dimension_vector):
            raise ValueError("object_vector and dimension_vector must have the same length.")
        if dimension_vector and not all(dim >= 0 for dim in dimension_vector):
            # Dimension cannot be non-negative
            raise ValueError("All elements of dimension_vector must be non-negative integers.") 
        
        if dimension_vector is None:
            dimension_vector = [1] * len(object_vector)
        

        ########
        # Set member variables
        ########
        
        self.object_vector = object_vector ## List of coherent sheaves in the complex, so that the chain complex can operate similar to a DenseVector.

        self.dimension_vector = dimension_vector ## List of the number of direct sums of each object in the complex.

        self.shift_vector = shift_vector ## List of homological shifts in the complex.

        # If an element of the complex has dimension 0, we can get rid of it using helper method
        self._remove_zeros_from_dimension_vector()
        # Combine repeated sheaves with same shift number in the complex, e.g. O^4[0] ⊕ O^4[0] = O^8[0]
        self._combine_repeats()

        self._initialized = True ## Flag to indicate that the object has been initialized



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
    
    def __contains__(self, other):

        # Case (c): Not a DerivedCategoryObject
        if not isinstance(other, DerivedCategoryObject):
            return False

        # Case (b): A non-graded DerivedCategoryObject
        if not isinstance(other, GradedCoproductObject):
            return other in self.object_vector

        # Case (a): A GradedCoproductObject
        # Step 1: check all objects in `other` appear in `self`
        if not all(obj in self.object_vector for obj in other.object_vector):
            return False

        # Step 2: for each object in `other`, match by object identity,
        # then check that shift and dimension are compatible
        for o_obj, o_shift, o_dim in zip(other.object_vector, other.shift_vector, other.dimension_vector):
            try:
                # Find index in self where the object matches
                idx = self.object_vector.index(o_obj)
            except ValueError:
                return False  # shouldn't happen due to all(...) check, but safe

            # Check shift and dimension compatibility
            if self.shift_vector[idx] != o_shift:
                return False
            if self.dimension_vector[idx] < o_dim:
                return False

        return True
    
    def __sub__(self, other) -> 'GradedCoproductObject':

        if isinstance(other, GradedCoproductObject):

            if not other in self:
                raise ValueError("The GradedCoproductObject we are subtracting has different underlying objects or too many copies of the same object.")
            
            new_object_vector = self.object_vector.copy()
            new_dimension_vector = self.dimension_vector.copy()

            for o_obj, o_dim in zip(other.object_vector,  other.dimension_vector):

                idx = new_object_vector.index(o_obj)
                new_dimension_vector[idx] -= o_dim

            
            return GradedCoproductObject(object_vector=new_object_vector,
                                        shift_vector=self.shift_vector,
                                        dimension_vector=new_dimension_vector)
        elif isinstance(other, DerivedCategoryObject):

            if not other in self:
                raise ValueError("The DerivedCategoryObject we are subtracting is not in the GradedCoproductObject.")
            
            new_dimension_vector = self.dimension_vector.copy()
            idx = self.object_vector.index(other)
            new_dimension_vector[idx] -= 1

            return GradedCoproductObject(object_vector=self.object_vector,
                                        shift_vector=self.shift_vector,
                                        dimension_vector=new_dimension_vector)
        else:
            raise TypeError("The object we are subtracting must be a DerivedCategoryObject or a GradedCoproductObject.")
    
    def get(self, shift : int, other_val = None):
        
        try:
            idx = self.shift_vector.index(shift)
            return self.object_vector[idx], self.dimension_vector[idx]
        except ValueError:
            # If the shift is not found, return None or other_val
            if other_val is not None:
                return other_val
            else:
                return None
            
            
    def __iter__(self):
        """
        Allows iteration over (object, shift, dimension) triples in the graded coproduct.
        """
        return iter(zip(self.object_vector, self.shift_vector, self.dimension_vector))
    
    def chernCharacter(self):
        r"""!
        Helper function to compute the Chern Character of the chain complex. The Chern Character of
        a chain complex is the alternating sum of the Chern Characters of the individual sheaves in
        the complex. Since the Chern character is additive, we may multiply the Chern Characters by
        the dimension of the object to represent direct sums of sheaves.

        \return ChernCharacter The Chern Character of the chain complex
        """

        ## use the alternating sum formula: Σ (-1)^shift * mult * ch(obj)
        terms = [
            (-1)**shift * mult * obj.chernCharacter()
            for obj, shift, mult in zip(self.object_vector, self.shift_vector, self.dimension_vector)
        ]
        if not terms:
            raise ValueError("Cannot compute Chern character of an empty complex.")

        return reduce(add, terms)

    
    
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

        filtered = [
            (obj, shift, dim)
            for obj, shift, dim in zip(self.object_vector, self.shift_vector, self.dimension_vector)
            if dim != 0
        ]
        self.object_vector, self.shift_vector, self.dimension_vector = map(list, zip(*filtered)) if filtered else ([], [], [])

    def _combine_repeats(self):
        r"""!
        Helper function to combine repeated sheaves in the complex. This is useful for simplifying
        the complex, as we can combine repeated sheaves into a single object with a larger dimension.
        This function specifically requires the __hash__ implementation for the DerivedCategoryObject
        """

        # Dictionary to hold combined dimensions for each (object, shift) pair
        combined = {}

        # Iterate over the tuples
        for obj, dim, shift in zip(self.object_vector, self.dimension_vector, self.shift_vector):
            key = (obj, shift)
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
    






