from src.DerivedCategory import DerivedCategoryObject
from src.DerivedCategory.SphericalTwist import SphericalTwist
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.GradedCoproductObject import GradedCoproductObject


class SphericalTwistSum(GradedCoproductObject):
    r"""!
    This class acts similar to the ChainComplex class for CoherentSheaf; specifically, when considering
    double (or any n>1) spherical twists, there are always triangles that the twist will fit into that are 
    not necessarily the defining triangle and can be obtained by applying the twist functor to the defining
    triangle of the individual twists --- such triangles will often involve the sum of shifts of spherical twists,
    so we need an implementation for keeping track of such objects in a single argument of the DistinguishedTriangle.

    For successive spherical twist applications, this entire procedure will likely need to be generalized since
    it is inefficient to encode DoubleSphericalTwistSum, TripeSphericalTwistSum, etc. as separate classes.

    """


    def __init__(self, line_bundle_pairs_vector, dimension_vector, shift_vector, degree=1):
        r"""!
        Initialize an instance of SphericalTwistSum with the specified line bundles. This class is used
        to represent objects of the form

        Tw_{a_1} O(b_1)^{n_1}[s_1] ⊕ Tw_{a_2} O(b_2)^{n_2}[s_2]}

        where the spherical twists are defined as the cone of the evaluation morphism

                Hom(O(a), O(b)) ⊗ O(a) ---->  O(b) ----> Tw_a O(b)

        This class is primarily used for DoubleSphericalTwist, since the second canonical triangle
        will often be a sum of spherical twists.

        \param list line_bundle_pairs_vector: A list of tuples where each tuple is a pair of line bundles (lb1, lb2) that the spherical twist is
                                         defined
        \param list dimension_vector A list of non-negative integers representing the number of times each spherical twist is applied
        \param list shift_vector A list of integers representing the shift of each spherical twist
        \param int degree An optional argument for the degree of the variety, which is relevant to the dimension of the
                       derived RHom space for K3 surfaces of picard rank 1. This does not affect the P1 or P2 implementations.
                       Default is 1.

        \throws TypeError If line_bundle_pairs_vector is not a list of tuples
        \throws TypeError If dimension_vector is not a list of integers
        \throws TypeError If shift_vector is not a list of integers
        \throws ValueError If line_bundle_pairs_vector is not a list of tuples with exactly two elements
        \throws ValueError If line_bundle_pairs_vector is not a list of tuples where both objects are line bundles
        \throws ValueError If dimension_vector is not a list of non-negative integers
        \throws ValueError If line_bundle_pairs_vector, dimension_vector, and shift_vector do not have the same length
        \throws ValueError If line_bundle_pairs_vector is empty
        """

        if not all(isinstance(x, tuple) for x in line_bundle_pairs_vector):
            raise TypeError("line_bundle_pairs_vector must be a list of tuples")
        if not all(len(x) == 2 for x in line_bundle_pairs_vector):
            raise ValueError("line_bundle_pairs_vector must be a list of tuples with exactly two elements")
        if not all(isinstance(x[0], LineBundle) and isinstance(x[1], LineBundle) for x in line_bundle_pairs_vector):
            raise TypeError("line_bundle_pairs_vector must be a list of tuples where both objects are line bundles")
        if not all(x[0].catagory == x[1].catagory for x in line_bundle_pairs_vector):
            raise ValueError("Line bundles must be defined on the same variety; not all catagories currently match.")
        
        self.catagory = line_bundle_pairs_vector[0][0].catagory

        if not all(x[0].catagory == self.catagory for x in line_bundle_pairs_vector):
            raise ValueError("All line bundles pairs must be defined on the same catagory")

        if not all(isinstance(x, int) for x in dimension_vector):
            raise TypeError("dimension_vector must be a list of integers")
        if not all( x >= 0 for x in dimension_vector):
            raise ValueError("dimension_vector must be a list of non-negative integers")
        if not all(isinstance(x, int) for x in shift_vector):
            raise TypeError("shift_vector must be a list of integers")
        
        
        if len(line_bundle_pairs_vector) != len(dimension_vector) or len(line_bundle_pairs_vector) != len(shift_vector):
            raise ValueError("line_bundle_pairs_vector, dimension_vector, and shift_vector must have the same length")
        
        if len(line_bundle_pairs_vector) == 0:
            raise ValueError("line_bundle_pairs_vector must have at least one element")
        
        self.line_bundle_pairs_vector = line_bundle_pairs_vector ## A list of tuples where each tuple is a pair of line bundles (lb1, lb2) that the spherical twist is defined around

        self.dimension_vector = dimension_vector ## A list of non-negative integers representing the number of times each spherical twist is applied

        self.shift_vector = shift_vector ## A list of integers representing the shift of each spherical twist

        self.degree = degree ## An optional argument for the degree of the K3 surface, if applicable


        self._remove_zeros_from_dimension_vector()
        self._combine_repeats()



    # def __str__(self):
    #     r"""!
    #     Returns a string representation of the spherical twist by printing the defining triangle. 
    #     The string representation is similar to that of the chain complex, where 2 lines are printed.
    #     The first line contains the number of times the spherical twist is applied, and the second line
    #     contains the actual twist. For example, the string representation of a spherical twist sum given
    #     by the data [(O(1), O(1)), [3], [-2]] would be

    #             ⊕3
    #     Tw_1 O(1)[-2]

    #     \return str A string representation of the spherical twist sum
    #     """

    #     # Initialize lists to hold each formatted component
    #     first_line_components = []
    #     second_line_components = []
        
    #     # Iterate over the zipped input lists
    #     for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
    #         # Format the second line component
    #         second_line_component = ""
    #         if s != 0:
    #             second_line_component = f"Tw_{lb1.degree} O({lb2.degree})[{s}]"
    #         else:
    #             second_line_component = f"Tw_{lb1.degree} O({lb2.degree})"

    #         second_line_components.append(second_line_component)
            
    #         # Calculate the position to place the '⊕n' above 'O'
    #         o_position = second_line_component.index(')')

    #         # Create a string with spaces up to the 'O' position, then add '⊕n'
    #         first_line_component = ""
    #         if n != 1:
    #             first_line_component = ' ' * o_position + f'⊕{n}'
    #         else:
    #             first_line_component = ' ' * o_position + '  '
    #         first_line_components.append(first_line_component)
            
        
    #     # Join all components with ' ⊕ ' separator
    #     first_line = ' '.join(first_line_components)
    #     second_line = ' ⊕ '.join(second_line_components)
        
    #     # Combine the two lines with a newline character
    #     if first_line.isspace():
    #         return second_line
    #     result = f"{first_line}\n{second_line}"
        
    #     return result
    
    

    
    def shift(self, n):
        r"""!
        Cohomological shift of the complex by some fixed amount. This is one of the main methods one wishes
        to override for the DerivedCategoryObject class, since the shift of a spherical twist sum is simply
        the sum of the shifts of the individual spherical twists.

        \param int n The amount to shift the complex by  

        \return SphericalTwistSum The shifted spherical twist sum
        """


        new_shift_vector = [x + n for x in self.shift_vector]

        return SphericalTwistSum(self.line_bundle_pairs_vector, self.dimension_vector, new_shift_vector, self.degree)
    
    def central_charge(self, *args):
        r"""!
        Method to compute the central charge of the spherical twist sum. Since all stability conditions are numerical
        for our current implementations, it suffices to compute the central charge of the Chern character

        \param args: The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return complex The central charge of the spherical twist sum as a complex number

        \throws TypeError If the args are not of the correct type
        \throws ValueError If the number of args is incorrect
        \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")

        central_charge = 0

        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            sph_twist = SphericalTwist(lb1, lb2, self.degree)
            central_charge += n * (-1)**s * sph_twist.central_charge(*args)
        
        return central_charge
    
    def get_lowest_shift_component(self):
        r"""!
        Similar to the ChainComplex class, we wish to recover the lowest shift component of the spherical twist sum.
        This is useful for computing the Harder-Narasimhan factors of the DoubleSphericalTwistSum class.

        \return SphericalTwistSum The spherical twist sum with the lowest shift
        """

        min_shift = float('inf')
        min_lb_pair = None
        min_dimension = None

        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            if s < min_shift:
                min_shift = s
                min_lb_pair = (lb1, lb2)
                min_dimension = n
        
        return SphericalTwistSum([min_lb_pair], [min_dimension], [min_shift], self.degree)
            
    def get_highest_shift_component(self):
        r"""!
        Similar to the ChainComplex class, we wish to recover the highest shift component of the spherical twist sum.
        This is useful for computing the Harder-Narasimhan factors of the DoubleSphericalTwistSum class.

        \return SphericalTwistSum The spherical twist sum with the highest shift
        """
        
        max_shift = float('-inf')
        max_lb_pair = None
        max_dimension = None

        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            if s > max_shift:
                max_shift = s
                max_lb_pair = (lb1, lb2)
                max_dimension = n
        
        return SphericalTwistSum([max_lb_pair], [max_dimension], [max_shift], self.degree)
    
    def is_semistable(self, *args):
        r"""!
        Method to check if the spherical twist sum is stable. The sum of shifts of distinct objects will generally
        never be stable unless the objects and shifts are all the same; for example, one can consider the trivial case
        where E is a stable object, but

        E[1] ---> E[1] ⊕ E ---> E

        destabilizes the direct sum. In general, the spherical twist sum is stable if the Harder-Narasimhan
        filtration is trivial, i.e. just the object itself. 

        \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return bool True if the spherical twist sum is stable, False otherwise

        \throws TypeError If the args are not of the correct type
        \throws ValueError If the number of args is incorrect
        \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
            
            
        return len(self.get_HN_factors_ordered(*args)) == 1 
            
    def get_HN_factors_ordered(self, *args):
        r"""!
        This is a slightly modified version of the get_HN_factors method, where the Harder-Narasimhan factors
        are ordered by phase. In particular, this method is a slight misnomer in the sense that it is not actually
        claiming that the list returned is the HN filtration of the object, but rather a concatenated list of all 
        the individual semistable factors. This is still nonetheless useful for computing the HN factors of the
        DoubleSphericalTwistSum class.

        \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return A list of tuples where the first element is a DerivedCategoryObject and the second element is a float
                    representing the phase of the object. The list is always returned in such a way that the largest phase
                    HN factor is first and smallest is last.

        \throws TypeError If the args are not of the correct type
        \throws ValueError If the number of args is incorrect
        \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        

        HN_factors = []
        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            sph_twist = SphericalTwist(lb1, lb2, self.degree)

            individual_twist_HN_factors = sph_twist.get_HN_factors(*args)
            for (obj, phase) in individual_twist_HN_factors:
                if isinstance(obj, ChainComplex):
                    HN_factors.append((ChainComplex(sheaf_vector=obj.sheaf_vector,
                                                    shift_vector=[shift + s for shift in obj.shift_vector],
                                                    dimension_vector=[dim * n for dim in obj.dimension_vector]),
                                                    phase + s))
                else:
                    new_chern = obj.chernCharacter() * int(n * (-1)**s)
                    HN_factors.append((DerivedCategoryObject(string=obj.update_string_by_shift(str(obj), s),
                                                            catagory=self.catagory,
                                                            chern_character=new_chern),
                                                                phase + s))
                    
        # Return the list in such a way that the highest phase comes first and the lowest phase comes last
        return sorted(HN_factors, key=lambda x: x[1], reverse=True)

            
    def _remove_zeros_from_dimension_vector(self):
        """!
        Helper function which iterates through the dimension vector, and if a certain Coherent sheaf
        is only included 0 times, we may effectively erase it.
        """
        for i in range(len(self.dimension_vector)):
            if i < len(self.dimension_vector) and self.dimension_vector[i] == 0:
                del self.line_bundle_pairs_vector[i]
                del self.shift_vector[i]
                del self.dimension_vector[i]

    def _combine_repeats(self):
        """!
        Helper function to combine repeated sheaves in the complex. This is useful for simplifying
        the complex, as we can combine repeated sheaves into a single sheaf with a larger dimension.
        """

        # Dictionary to hold combined dimensions for each (sheaf, shift) pair
        combined = {}

        # Iterate over the tuples
        for (lb1,lb2), dim, shift in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            key = (lb1, lb2, shift)
            if key in combined:
                combined[key] += dim
            else:
                combined[key] = dim

        # Unpack the combined dictionary back into the class variables
        self.line_bundle_pairs_vector = []
        self.dimension_vector = []
        self.shift_vector = []

        for (lb1, lb2, shift), dim in combined.items():
            self.line_bundle_pairs_vector.append((lb1, lb2))
            self.dimension_vector.append(dim)
            self.shift_vector.append(shift)


    