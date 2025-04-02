from src.DerivedCategory import DerivedCategoryObject
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition
from src.DerivedCategory.GradedCoproductObject import GradedCoproductObject, LineBundleCoproduct

from typing import List


class SphericalTwistCoproduct(GradedCoproductObject):
    r"""!
    This class acts similar to the ChainComplex class for CoherentSheaf; specifically, when considering
    double (or any n>1) spherical twists, there are always triangles that the twist will fit into that are 
    not necessarily the defining triangle and can be obtained by applying the twist functor to the defining
    triangle of the individual twists --- such triangles will often involve the sum of shifts of spherical twists,
    so we need an implementation for keeping track of such objects in a single argument of the DistinguishedTriangle.

    For successive spherical twist applications, this entire procedure will likely need to be generalized since
    it is inefficient to encode DoubleSphericalTwistSum, TripeSphericalTwistSum, etc. as separate classes.

    """


    def __init__(self, sph_twists_vector : List[SphericalTwistComposition], shift_vector : List[int], dimension_vector : List[int] = None, degree_K3 : int = 1):
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
        """


        super().__init__(object_vector=sph_twists_vector, dimension_vector=dimension_vector, shift_vector=shift_vector, degree=degree_K3)

        if not all(isinstance(x, SphericalTwistComposition) for x in sph_twists_vector):
            raise TypeError("All objects in the spherical twist sum must be instances of SphericalTwistComposition")


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
    
    

    
    
    def get_lowest_shift_component(self):
        r"""!
        Similar to the ChainComplex class, we wish to recover the lowest shift component of the spherical twist sum.
        This is useful for computing the Harder-Narasimhan factors of the DoubleSphericalTwistSum class.

        \return SphericalTwistSum The spherical twist sum with the lowest shift
        """

        # Bundle each component into a tuple, then use max with custom key
        triple = min(
            zip(self.object_vector, self.dimension_vector, self.shift_vector),
            key=lambda t: t[2]  # t = (lb_pair, dim, shift)
        )

        max_lb_pair, max_dimension, max_shift = triple
        return SphericalTwistCoproduct([max_lb_pair], [max_dimension], [max_shift], self.degree)

            
    def get_highest_shift_component(self):
        r"""!
        Similar to the ChainComplex class, we wish to recover the highest shift component of the spherical twist sum.
        This is useful for computing the Harder-Narasimhan factors of the DoubleSphericalTwistSum class.

        \return SphericalTwistSum The spherical twist sum with the highest shift
        """
        
        # Bundle each component into a tuple, then use max with custom key
        triple = max(
            zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector),
            key=lambda t: t[2]  # t = (lb_pair, dim, shift)
        )

        max_lb_pair, max_dimension, max_shift = triple
        return SphericalTwistCoproduct([max_lb_pair], [max_dimension], [max_shift], self.degree)

    
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
            sph_twist = SphericalTwistComposition(lb1, lb2, self.degree)

            individual_twist_HN_factors = sph_twist.get_HN_factors(*args)
            for (obj, phase) in individual_twist_HN_factors:
                if isinstance(obj, LineBundleCoproduct):
                    HN_factors.append((LineBundleCoproduct(sheaf_vector=obj.sheaf_vector,
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

    