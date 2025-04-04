from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject

from typing import List

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
    






# class LineBundleCoproduct(GradedCoproductObject):
    

#     def __init__(self, line_bundle_vector : List[LineBundle], shift_vector : List[int], dimension_vector : List[int] = None, degree_K3 : int = 1):
#         r"""!
#         Initialize an instance of LineBundleCoproduct with the specified sheaf vector, shift vector,
#         and potentially a dimension vector. If a dimension vector is not provided, it must 
#         consist of non-negative integer values


#         \param list line_bundle_vector A list of coherent sheaves in the complex
#         \param list shift_vector A list of homological shifts in the complex
#         \param list dimension_vector A list of the number of direct sums of each coherent sheaf in the complex
#         """

#         # The main functionality is implemented in the parent class
#         super().__init__(object_vector=line_bundle_vector, shift_vector=shift_vector, dimension_vector=dimension_vector, degree_K3=degree_K3)

#         # Check that the sheaf vector is valid and contains sheaves of only a single catagory
#         if not all(isinstance(sheaf, LineBundle) for sheaf in line_bundle_vector):
#             raise TypeError("All elements of line_bundle_vector must be instances of LineBundle.")

    


#     def isShiftOfSheaf(self):
#         r"""!
#         Simple helper function which checks if the complex is a shift of a single sheaf

#         \return bool True if the complex is a shift of a single sheaf, False otherwise
#         """

#         return len(self.object_vector) == 1


        
    
    



    # def is_semistable(self, *args):
    #     r"""!
    #     Method to compute whether the chain complex is semistable. This almost never occurs, since if
    #     the complex contains two or more stable objects of distinct phase, it will never be stable. For
    #     example, suppose E2 is a stable subobject of maximum phase and E1 is another stable object with
    #     strictly smaller phase. Then

    #                       E2 ----> Complex ------> Cone

    #     will destabilize the complex, and Cone will be nontrivial since it has a non-zero map to E1. 
    #     The easiest way to check that the complex is concentrated in only a single phase is to compare
    #     its largest and smallest phases from the previous methods.


    #     \param tuple args The arguments required to compute the phase. The number of arguments and the type of arguments will
    #                  depend on the catagory of the sheaves in the complex. For P1, the phase requires a single complex number.
    #                  For P2, the phase requires two floating-point numbers. For K3, the phase requires two floating-point
    #                  numbers and an integer.

    #     \return bool True if the chain complex is semistable, False otherwise

    #     \throws ValueError If the number of arguments is incorrect
    #     \throws TypeError If the type of the arguments is incorrect
    #     \throws NotImplementedError If the catagory of the sheaves in the complex is not implemented
    #     """

    
    #     if self.catagory == 'P1':
    #         if len(args) != 1:
    #             raise ValueError("Phase of P1 requires exactly one argument.")
    #         if not isinstance(args[0], complex):
    #             raise TypeError("Phase of P1 requires a complex number as an argument.")
    #     elif self.catagory == 'P2':
    #         if len(args) != 2:
    #             raise ValueError("Phase of P2 requires exactly two arguments.")
    #         if not all(isinstance(arg, (float,int)) for arg in args):
    #             raise TypeError("Phase of P2 requires two floating-point numbers as arguments.")
    #     elif self.catagory == 'K3':
    #         if len(args) != 3:
    #             raise ValueError("Phase of K3 requires exactly three arguments.")
    #         if not all(isinstance(arg, (float,int)) for arg in args):
    #             raise TypeError("Phase of K3 requires three floating-point numbers as arguments.")
    #         if not isinstance(args[2], int):
    #             raise TypeError("The degree of the K3 surface must be an integer.")
    #     else:
    #         raise NotImplementedError("Only local P1, local P2, and K3 catagories are implemented.")

    #     return self.get_largest_phase(*args) == self.get_smallest_phase(*args)
    

    # def get_HN_factors(self, *args) ->  List[Tuple[DerivedCategoryObject, float]]:
        
    #     return_list = []
    #     for sheaf, shift in zip(self.object_vector, self.shift_vector):
    #         if not sheaf.is_stable(*args):
    #             raise ValueError("Sheaf is not stable.")
            
    #         phase = sheaf.phase(*args) + shift
    #         return_list.append((sheaf, phase))
    #     return_list.sort(key=lambda x: x[1], reverse=True)
    #     return return_list


# class SphericalTwistCoproduct(GradedCoproductObject):
#     r"""!
#     This class acts similar to the ChainComplex class for CoherentSheaf; specifically, when considering
#     double (or any n>1) spherical twists, there are always triangles that the twist will fit into that are 
#     not necessarily the defining triangle and can be obtained by applying the twist functor to the defining
#     triangle of the individual twists --- such triangles will often involve the sum of shifts of spherical twists,
#     so we need an implementation for keeping track of such objects in a single argument of the DistinguishedTriangle.

#     For successive spherical twist applications, this entire procedure will likely need to be generalized since
#     it is inefficient to encode DoubleSphericalTwistSum, TripeSphericalTwistSum, etc. as separate classes.

#     """



#     def __init__(self, sph_twists_vector : List[SphericalTwistComposition], shift_vector : List[int], dimension_vector : List[int] = None, degree_K3 : int = 1):
#         r"""!
#         Initialize an instance of SphericalTwistSum with the specified line bundles. This class is used
#         to represent objects of the form

#         Tw_{a_1} O(b_1)^{n_1}[s_1] ⊕ Tw_{a_2} O(b_2)^{n_2}[s_2]}

#         where the spherical twists are defined as the cone of the evaluation morphism

#                 Hom(O(a), O(b)) ⊗ O(a) ---->  O(b) ----> Tw_a O(b)

#         This class is primarily used for DoubleSphericalTwist, since the second canonical triangle
#         will often be a sum of spherical twists.

#         \param list line_bundle_pairs_vector: A list of tuples where each tuple is a pair of line bundles (lb1, lb2) that the spherical twist is
#                                          defined
#         \param list dimension_vector A list of non-negative integers representing the number of times each spherical twist is applied
#         \param list shift_vector A list of integers representing the shift of each spherical twist
#         \param int degree An optional argument for the degree of the variety, which is relevant to the dimension of the
#                        derived RHom space for K3 surfaces of picard rank 1. This does not affect the P1 or P2 implementations.
#                        Default is 1.
#         """


#         super().__init__(object_vector=sph_twists_vector, dimension_vector=dimension_vector, shift_vector=shift_vector, degree=degree_K3)


#         if not all(isinstance(x, SphericalTwistComposition) for x in sph_twists_vector):
#             raise TypeError("All objects in the spherical twist sum must be instances of SphericalTwistComposition")


#     # def __str__(self):
#     #     r"""!
#     #     Returns a string representation of the spherical twist by printing the defining triangle. 
#     #     The string representation is similar to that of the chain complex, where 2 lines are printed.
#     #     The first line contains the number of times the spherical twist is applied, and the second line
#     #     contains the actual twist. For example, the string representation of a spherical twist sum given
#     #     by the data [(O(1), O(1)), [3], [-2]] would be

#     #             ⊕3
#     #     Tw_1 O(1)[-2]

#     #     \return str A string representation of the spherical twist sum
#     #     """

#     #     # Initialize lists to hold each formatted component
#     #     first_line_components = []
#     #     second_line_components = []
        
#     #     # Iterate over the zipped input lists
#     #     for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
#     #         # Format the second line component
#     #         second_line_component = ""
#     #         if s != 0:
#     #             second_line_component = f"Tw_{lb1.degree} O({lb2.degree})[{s}]"
#     #         else:
#     #             second_line_component = f"Tw_{lb1.degree} O({lb2.degree})"

#     #         second_line_components.append(second_line_component)
            
#     #         # Calculate the position to place the '⊕n' above 'O'
#     #         o_position = second_line_component.index(')')

#     #         # Create a string with spaces up to the 'O' position, then add '⊕n'
#     #         first_line_component = ""
#     #         if n != 1:
#     #             first_line_component = ' ' * o_position + f'⊕{n}'
#     #         else:
#     #             first_line_component = ' ' * o_position + '  '
#     #         first_line_components.append(first_line_component)
            
        
#     #     # Join all components with ' ⊕ ' separator
#     #     first_line = ' '.join(first_line_components)
#     #     second_line = ' ⊕ '.join(second_line_components)
        
#     #     # Combine the two lines with a newline character
#     #     if first_line.isspace():
#     #         return second_line
#     #     result = f"{first_line}\n{second_line}"
        
#     #     return result
    
    

    
    
#     def get_lowest_shift_component(self):
#         r"""!
#         Similar to the ChainComplex class, we wish to recover the lowest shift component of the spherical twist sum.
#         This is useful for computing the Harder-Narasimhan factors of the DoubleSphericalTwistSum class.

#         \return SphericalTwistSum The spherical twist sum with the lowest shift
#         """

#         # Bundle each component into a tuple, then use max with custom key
#         triple = min(
#             zip(self.object_vector, self.dimension_vector, self.shift_vector),
#             key=lambda t: t[2]  # t = (lb_pair, dim, shift)
#         )

#         max_lb_pair, max_dimension, max_shift = triple
#         return SphericalTwistCoproduct([max_lb_pair], [max_dimension], [max_shift], self.degree)

            
#     def get_highest_shift_component(self):
#         r"""!
#         Similar to the ChainComplex class, we wish to recover the highest shift component of the spherical twist sum.
#         This is useful for computing the Harder-Narasimhan factors of the DoubleSphericalTwistSum class.

#         \return SphericalTwistSum The spherical twist sum with the highest shift
#         """
        
#         # Bundle each component into a tuple, then use max with custom key
#         triple = max(
#             zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector),
#             key=lambda t: t[2]  # t = (lb_pair, dim, shift)
#         )

#         max_lb_pair, max_dimension, max_shift = triple
#         return SphericalTwistCoproduct([max_lb_pair], [max_dimension], [max_shift], self.degree)

    
#     def is_semistable(self, *args):
#         r"""!
#         Method to check if the spherical twist sum is stable. The sum of shifts of distinct objects will generally
#         never be stable unless the objects and shifts are all the same; for example, one can consider the trivial case
#         where E is a stable object, but

#         E[1] ---> E[1] ⊕ E ---> E

#         destabilizes the direct sum. In general, the spherical twist sum is stable if the Harder-Narasimhan
#         filtration is trivial, i.e. just the object itself. 

#         \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
#                         For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
#                         and an integer representing the degree of the K3 surface.

#         \return bool True if the spherical twist sum is stable, False otherwise

#         \throws TypeError If the args are not of the correct type
#         \throws ValueError If the number of args is incorrect
#         \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
#         """

#         if self.catagory == 'P1':
#             if len(args) != 1:
#                 raise ValueError("Central charge of P1 requires single complex number parameter")
#             if not isinstance(args[0], complex):
#                 raise TypeError("P1 objects should have a single complex parameter")
#         elif self.catagory == 'P2':
#             if len(args) != 2:
#                 raise ValueError("Central charge of P2 requires two real number parameters")
#             if not all(isinstance(x, (float, int)) for x in args):
#                 raise TypeError("P2 objects should have two real number parameters")
#         elif self.catagory == 'K3':

#             if len(args) != 3:
#                 raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
#             if not all(isinstance(x, (float, int)) for x in args):
#                 raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
#             if not isinstance(args[2], int):
#                 raise TypeError("The degree of the K3 surface must be an integer")
#         else:
#             raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
            
            
#         return len(self.get_HN_factors_ordered(*args)) == 1 
            
#     def get_HN_factors_ordered(self, *args):
#         r"""!
#         This is a slightly modified version of the get_HN_factors method, where the Harder-Narasimhan factors
#         are ordered by phase. In particular, this method is a slight misnomer in the sense that it is not actually
#         claiming that the list returned is the HN filtration of the object, but rather a concatenated list of all 
#         the individual semistable factors. This is still nonetheless useful for computing the HN factors of the
#         DoubleSphericalTwistSum class.

#         \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
#                         For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
#                         and an integer representing the degree of the K3 surface.

#         \return A list of tuples where the first element is a DerivedCategoryObject and the second element is a float
#                     representing the phase of the object. The list is always returned in such a way that the largest phase
#                     HN factor is first and smallest is last.

#         \throws TypeError If the args are not of the correct type
#         \throws ValueError If the number of args is incorrect
#         \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
#         """

#         if self.catagory == 'P1':
#             if len(args) != 1:
#                 raise ValueError("Central charge of P1 requires single complex number parameter")
#             if not isinstance(args[0], complex):
#                 raise TypeError("P1 objects should have a single complex parameter")
#         elif self.catagory == 'P2':
#             if len(args) != 2:
#                 raise ValueError("Central charge of P2 requires two real number parameters")
#             if not all(isinstance(x, (float, int)) for x in args):
#                 raise TypeError("P2 objects should have two real number parameters")
#         elif self.catagory == 'K3':

#             if len(args) != 3:
#                 raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
#             if not all(isinstance(x, (float, int)) for x in args):
#                 raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
#             if not isinstance(args[2], int):
#                 raise TypeError("The degree of the K3 surface must be an integer")
#         else:
#             raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        

#         HN_factors = []
#         for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
#             sph_twist = SphericalTwistComposition(lb1, lb2, self.degree)

#             individual_twist_HN_factors = sph_twist.get_HN_factors(*args)
#             for (obj, phase) in individual_twist_HN_factors:
#                 if isinstance(obj, LineBundleCoproduct):
#                     HN_factors.append((LineBundleCoproduct(sheaf_vector=obj.sheaf_vector,
#                                                     shift_vector=[shift + s for shift in obj.shift_vector],
#                                                     dimension_vector=[dim * n for dim in obj.dimension_vector]),
#                                                     phase + s))
#                 else:
#                     new_chern = obj.chernCharacter() * int(n * (-1)**s)
#                     HN_factors.append((DerivedCategoryObject(string=obj.update_string_by_shift(str(obj), s),
#                                                             catagory=self.catagory,
#                                                             chern_character=new_chern),
#                                                                 phase + s))
                    
#         # Return the list in such a way that the highest phase comes first and the lowest phase comes last
#         return sorted(HN_factors, key=lambda x: x[1], reverse=True)

            