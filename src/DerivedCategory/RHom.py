from src.DerivedCategory.SphericalTwist import SphericalTwistComposition
from src.DerivedCategory.CoherentSheaf import LineBundle, CoherentSheaf
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject, ZeroObject

from typing import Dict
import math




class LongExactSequenceException(Exception):
    r"""!
    Exception raised when the long exact sequence cannot be computed
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args)

        self.message = kwargs.get('message') ## The error message

        self.sequence_str = kwargs.get('sequence_str') ## The string representation of the sequence




def RHom(object1 : DerivedCategoryObject, object2 : DerivedCategoryObject) -> Dict[int,int]:
    r"""!
    Primary function for computing the RHom space between two derived category objects, as a graded C-vector space. The actual behavior of this
    wrapper method is handled by a series of helper methods, which implement the logic for computing the RHom space between
    various subclasses of derived category objects. The most computationally expensive computation is finding the RHom space between
    two SphericalTwistCompositions; any time a SphericalTwistComposition is done, the function makes roughly 2**N recursive subcalls where
    N is the number of line bundles twisted around in a single parameter. 

    \param object1 A DerivedCategoryObject representing the first object in the RHom space
    \param object2 A DerivedCategoryObject representing the second object in the RHom space

    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space

    \throws NotImplementedError If the DerivedCategoryObject is not a LineBundle, SphericalTwistComposition, GradedCoproductObject, or ZeroObject
    """


    if isinstance(object1, LineBundle) and isinstance(object2, LineBundle):
        ## We are computing the RHom space between two line bundles; this
        ## should just be handled by the standard Ext group of the respective
        ## abelian category
        return Ext(object1, object2)
    
    elif isinstance(object1, LineBundle) and isinstance(object2, SphericalTwistComposition):
        return _compute_rhom_line_bundle_to_sph_helper(object1, object2)
    elif isinstance(object1, SphericalTwistComposition) and isinstance(object2, LineBundle):
        return _compute_rhom_sph_to_line_bundle_helper(object1, object2)
    elif isinstance(object1, GradedCoproductObject):
        
        if isinstance(object2, LineBundle) or isinstance(object2, SphericalTwistComposition):
            return _compute_rhom_graded_coproduct_to_derived_ob_helper(object1, object2)
        elif isinstance(object2, GradedCoproductObject):
            return _compute_rhom_graded_coproduct_to_graded_coproduct_helper(object1, object2)
        elif isinstance(object2, ZeroObject):
            return {0 : 0}
        else:
            raise NotImplementedError(f"Cannot compute RHom between {type(object1)} and {type(object2)}")
    elif isinstance(object2, GradedCoproductObject):
        if isinstance(object1, LineBundle) or isinstance(object1, SphericalTwistComposition):
            return _compute_rhom_derived_ob_to_graded_coproduct_helper(object1, object2)
        elif isinstance(object1, ZeroObject):
            return {0 : 0}
        else:
            raise NotImplementedError(f"Cannot compute RHom between {type(object1)} and {type(object2)}")
    elif isinstance(object1, ZeroObject) or isinstance(object2, ZeroObject):
        ## If one of the objects is a zero object, then the RHom space is just the zero object
        return {0 : 0}

    else:
        raise NotImplementedError(f"Cannot compute RHom between {type(object1)} and {type(object2)}")
    





def Ext(object1 : CoherentSheaf, object2 : CoherentSheaf) -> Dict[int,int]:
    r"""!
    Primary function for computing the Ext space between two coherent sheaves. While the Ext space is technically the 
    derived Hom space for objects which lie in the heart of the standard t-structure, it is easier to treat one as a 
    special case of the other (especially since we consider CoherentSheaf and DerivedCategoryObject two different modules).

    One important characteristic is that the dictionary returned should only be concentrated in degrees 0 to (-1 times) the dimension of the
    underlying variety. \
    
    \param object1 A CoherentSheaf representing the first object in the Ext space
    \param object2 A CoherentSheaf representing the second object in the Ext space

    \return dict A dictionary representing the dimensions of the Ext space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
    \throws NotImplementedError If the CoherentSheaf is not a LineBundle. Currently, we only have the logic of Hirzebruch-Riemann-Roch to implement the Ext space
    between two line bundles. We hope to implement the general case in the future.
    """

    if isinstance(object1, LineBundle) and isinstance(object2, LineBundle):
        return _ext_line_bundles(object1, object2)
    else:
        # TODO: Figure out how to implement for more general coherent sheaves
        return NotImplementedError("Currently, Ext(A, B) is only implemented for line bundles. \n\t\t\t We hope to implement general coherent sheaves in the future.")






def _ext_line_bundles(line_bundle_1 : LineBundle, line_bundle_2 : LineBundle) -> Dict[int,int]:
    r"""!

    Helper method which implements the logic for computing the Ext space between two line bundles. The logic is based on Hirzebruch-Riemann-Roch, and 
    checks the different GeometryContexts possible for the line bundles. Since we utilize SymPy for the combinatorial logic, the method is not as 
    efficient as if purely numerical invariants were used.
    
    """


    if not isinstance(line_bundle_1, LineBundle) or not isinstance(line_bundle_2, LineBundle):
        raise TypeError("Both objects must be LineBundles")
    if line_bundle_1.geometry_context != line_bundle_2.geometry_context:
        raise ValueError("LineBundles must be defined on the same underlying variety")
    
    if line_bundle_1.geometry_context.catagory == 'LocalP1':
        # The dimensions of the pushforwards are computed using the triangle
        #   
        #       i^* i_* E -> E -> E x O(2)[2]
        #
        # and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        # the Hom spaces between line bundles on P1, the computation reduces to a case-by-case
        # combinatorial problem. Since the homological index of the hom-space on P1 is bounded between
        # 0 and 1, the hom-space for local P1 is concentrated between degrees 0 and 2.


        lb1_coeff = line_bundle_1.getCoefficient(line_bundle_1.geometry_context.polarization)
        lb2_coeff = line_bundle_2.getCoefficient(line_bundle_2.geometry_context.polarization)

        degree_dif = lb2_coeff - lb1_coeff

        if degree_dif == 0:
            return {0 : 1, -2: 1}
        if degree_dif >= 2:
            return {0 : degree_dif + 1, -1 : degree_dif - 1}
        elif degree_dif == 1:
            return {0 : 2}
        elif degree_dif == -1:
            return {-2 : 2}
        else:
            return {-1 : -1*degree_dif - 1, -2 : -1*degree_dif + 1}
    elif line_bundle_1.geometry_context.catagory == 'LocalP2':
        ##  The dimensions of the pushforwards are computed using the triangle
        ##
        ##      i^* i_* E -> E -> E x O(3)[2]
        ##
        ##  and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        ##  the Hom spaces between line bundles on P^2, the computation reduces to a case-by-case
        ##  combinatorial problem. Since the homological index of the hom-space on P^2 is bounded between
        ##  0 and 2, the hom-space for local P2 is concentrated between degrees 0 and 3.

        lb1_coeff = line_bundle_1.getCoefficient(line_bundle_1.geometry_context.polarization)
        lb2_coeff = line_bundle_2.getCoefficient(line_bundle_2.geometry_context.polarization)
        degree_dif = line_bundle_2 - line_bundle_1

        if degree_dif == 0:
            return {0 : 1,  -3 : 1}
        elif degree_dif > -3 and degree_dif < 0:
            rank3 = math.comb(line_bundle_1 - line_bundle_2 + 2, 2)
            return {-3 : rank3}
        elif degree_dif > 0 and degree_dif < 3:
            rank0 = math.comb(degree_dif + 2, 2)
            return {0 : rank0}
        elif degree_dif >= 3:
            rank0 = math.comb(degree_dif + 2, 2)
            rank1 = math.comb(degree_dif - 1, 2)
            return {0 : rank0, -1 : rank1}
        else:
            rank2 = math.comb(line_bundle_1 - line_bundle_2 -1, 2)
            rank3 = math.comb(line_bundle_1 - line_bundle_2 + 2, 2)   
            return {-2 : rank2, -3 : rank3} 
        
    elif line_bundle_1.geometry_context.catagory == 'K3':
        ##  Serre duality implies that the complex must
        ##  be concentrated in a single degree, which is either 0 or 2 and corresponds to the cases that b > a and
        ##   a < b, respectively. 
        ##  
        ##  We return a tuple of the form (n0,n1,n2) indicating the dimension of the graded RHom space. Aside from the
        ##  instance that a=b, this will only include one nonzero element. Hirzebruch-Riemann-Roch shows that
        ##  
        ##  dim RHom(O(a), O(b)) = 2 + d(b-a)^2 

        D1 = line_bundle_1.divisor
        D2 = line_bundle_2.divisor
        D = D2 - D1
        
        D_squared = line_bundle_2.geometry_context.divisor_data.evaluate(D, D)
        chi = int(2 + D_squared / 2)

        if D == 0:
            return {0: 1 , -2: 1}
        elif line_bundle_2.geometry_context.divisor_data.is_effective(D):
            return {0 : chi}
        elif line_bundle_2.geometry_context.divisor_data.is_effective(-D):
            return {-2 : chi}
        else:
            raise ValueError(f"Divisor arising from {line_bundle_2} - {line_bundle_1} is neither effective nor anti-effective. Currrently there is no way to compute this for K3 surfaces.")



def _compute_rhom_graded_coproduct_to_graded_coproduct_helper(gc1 : GradedCoproductObject, gc2 : GradedCoproductObject) -> Dict[int, int]:


    direct_sum_dict = {}
    for obj1, shift1, dim1 in zip(gc1.object_vector, gc1.shift_vector, gc1.dimension_vector):

        hom_dict = {}

        if isinstance(obj1, LineBundle) or isinstance(obj1, SphericalTwistComposition):
            hom_dict = _compute_rhom_derived_ob_to_graded_coproduct_helper(obj1, gc2)
        else:
            raise NotImplementedError(f"RHom is not yet implemented for {type(obj1)} --> {type(gc2)}")
        
        temp_dict = _shift_dict(_multiply_dict(hom_dict, dim1), -1*shift1)
        direct_sum_dict = _add_dicts(direct_sum_dict, temp_dict)

    return direct_sum_dict






def _compute_rhom_graded_coproduct_to_derived_ob_helper(gc : GradedCoproductObject, d_ob : DerivedCategoryObject) -> Dict[int, int]:

    r"""!
    Helper method which implements the general homological algebra for computing the right-derived Hom space of a trivial direct sum of 
    DerivedCategoryObject objects with another DerivedCategoryObject. In paricular, this method simply encodes the fact that the derived Hom splits
    with respect to direct sums and commutes with shifts — this we may extract the sum and shift data from the first parameter and pass it 
    to some result obtained by applying RHom to the individual objects. 

    \param GradedCoproductObject gc The graded coproduct object which is the first argument of RHom(-, -)
    \param DerivedCategoryObject d_ob The derived category object which is the second argument of RHom(-, -)
    
    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
    """


    direct_sum_dict = {}
    for obj, shift, dim in zip(gc.object_vector, gc.shift_vector, gc.dimension_vector):

        hom_dict = RHom(obj, d_ob)

        temp_dict = _shift_dict(_multiply_dict(hom_dict, dim), -1*shift)
        direct_sum_dict = _add_dicts(direct_sum_dict, temp_dict)

    return direct_sum_dict








def _compute_rhom_derived_ob_to_graded_coproduct_helper(d_ob : DerivedCategoryObject, gc : GradedCoproductObject) -> Dict[int, int]:
    r"""!
    Helper method which implements the general homological algebra for computing the right-derived Hom space of a general object with
    a trivial direct sum of DerivedCategoryObject objects. In paricular, this method simply encodes the fact that the derived Hom splits
    with respect to direct sums and commutes with shifts — this we may extract the sum and shift data from the second parameter and pass it 
    to some result obtained by applying RHom to the individual objects. 

    \param DerivedCategoryObject d_ob The derived category object which is the first argument of RHom(-, -)
    \param GradedCoproductObject gc The graded coproduct object which is the second argument of RHom(-, -)

    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
    """


    direct_sum_dict = {}
    for obj, shift, dim in zip(gc.object_vector, gc.shift_vector, gc.dimension_vector):

        hom_dict = RHom(d_ob, obj)

        temp_dict = _shift_dict(_multiply_dict(hom_dict, dim), shift)
        direct_sum_dict = _add_dicts(direct_sum_dict, temp_dict)

    return direct_sum_dict
        




def _compute_rhom_line_bundle_to_sph_helper(lb : LineBundle, sph : SphericalTwistComposition) -> Dict[int, int]:
    r"""!
    Recursive helper method which implements the general homological-algebraic logic for computing the right-derived Hom space
    of a line bundle with a spherical twist composition. The base case logic is handled by the general RHom wrapper. This method
    primarily employs long-exact sequence logic to try to compute the dimensions of the RHom space when possible. Ultimately, the 
    method is much less robust than one would hope for — even when short exact sequences exist in the long exact cohomology sequence,
    they sometimes cause LongExactSequenceExceptions due to the fact that we mathematically cant assume 0 -> ? -> B -> C -> ? -> 0 
    gives any information about the outside two terms.

    \param LineBundle lb The first argument of RHom(-, -), which is assumed to be a LineBundle for this helper method
    \param SphericalTwistComposition sph The second argument of RHom(-, -), which is assumed to be a SphericalTwistComposition for this helper method

    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
    
    \throws LongExactSequenceException If we encounter a sequence of non-zero terms in our Long exact sequence that do not yield a definitive answer for the dimension of our desired RHom space
    """


    ################################################################################################################
    #                                                                                                              #       
    #                      Recursive Calls from applying RHom(O(a_n) , ----) to                                    #
    #                                                                                                              #
    #      O(a_{n-1}) ---> Tw_{a_{n-2}}... Tw_{a_1} O(a_0) ----> Tw_{a_{n-1}} Tw_{a_{n-2}}... Tw_{a_1} O(a_0)      #
    #                                                                                                              #
    ################################################################################################################
    

    ## The result of the previous step RHom( O(a_{n-1}) ,  Tw_{a_{n-2}} ... Tw_{a_1} O(a_0) ) will directly
    ## be applied to RHom( O(a_n) , O(a_{n-1}) )
    first_term_dict = RHom(lb, sph.defining_triangle[0])
    

    ## This is a second call to recursion that does not complete in constant time. In particular, this means our function
    ## must make roughly 2^n calls to itself, where n is the number of elements in the sequence. The runtime can be improved
    ## via memoization, but this is not implemented here.
    middle_term_dict = RHom(lb, sph.defining_triangle[1])


    ###########################################################################
    #                                                                         #
    #               Attempt to resolve long-exact sequence                    #
    #                                                                         #
    ###########################################################################



    # Combine first_term_dict and middle_term_dict into return_dict
    return_dict = {}

    raw_keys = set(first_term_dict) | set(middle_term_dict)

    ######
    # Extend the cohomological degrees we look over to be 1 less and 1 more than the keys from other terms
    ######

    if raw_keys:  # guard against empty dicts
        extended_keys = raw_keys | {min(raw_keys) - 1, max(raw_keys) + 1}
    else:
        extended_keys = set()

    keys = sorted(extended_keys, reverse=True)

    for k in keys:
        ## At each iteration, we must check that we are not overwriting a previous value
        ## in returned dictionary. Since there are only three terms in a triangle contributing
        ## to the long exact sequence, the only possible way that a value could be overwritten is
        ## if there is a long-exact sequence with at least 4 consecutive non-zero terms.
        ##
        ## (Resolving such an error would require more specialized care outside of pure diagram chasing;
        ## in particular, one would need to know the explicit maps in the long exact sequence.)
        ##
        ## Algorithmically, this is simply handled by saving the index we wish to overwrite to and 
        ## the value stored, and then checking to make sure that index is not already in the dictionary.


        long_ex_str = f"""
                [{k+1}] \t:\t {first_term_dict.get(k+1,0)} ---> {middle_term_dict.get(k+1,0)} ---> {return_dict.get(k+1,0)}
                [{k}] \t:\t {first_term_dict.get(k,0)}   ---> {middle_term_dict.get(k,0)} ---> ?
                [{k-1}] \t:\t {first_term_dict.get(k-1, 0)} ---> {middle_term_dict.get(k-1,0)} ---> ?
                """


        if middle_term_dict.get(k, 0) != 0:
            ##                                      ((prev set))
            ##           A[k+1]  ----> B[k+1] ------> C[k+1]
            ##           A[k]    ----> B[k]   ------> (sum of surrounding terms)
            ##           A[k-1]  ----> B[k-1] ------>    ?                   
            ##

            if first_term_dict.get(k, 0) != 0:

                ## There are exactly two kinds of long-exact sequences we can resolve; otherwise raise 
                ## an exception.

                if return_dict.get(k+1, 0) == 0 and middle_term_dict.get(k-1, 0) == 0:
                    ##                                      ((prev set))
                    ##           A[k+1]  ----> B[k+1] ------> 0
                    ##           A[k]    ----> B[k]   ------> (sum of surrounding terms)
                    ##           A[k-1]  ----> 0 ------>    ?                   
                    ##
                    ##  This four-term exact sequence can be resolved since three of the dimensions are known

                    return_dict[k] = middle_term_dict.get(k, 0) - first_term_dict.get(k, 0) + first_term_dict.get(k-1, 0)
                    continue

                elif middle_term_dict.get(k+1, 0) == 0 and first_term_dict.get(k-1, 0) == 0:

                    ##                                      ((prev set))
                    ##           A[k+1]  ----> 0 ------> C[k+1]
                    ##           A[k]    ----> B[k]   ------> (sum of surrounding terms)
                    ##           0       ----> B[k-1] ------>    ?                   
                    ##
                    ##  This four-term exact sequence can be resolved since three of the dimensions are known

                    return_dict[k] = return_dict.get(k+1, 0) + middle_term_dict.get(k, 0) - first_term_dict.get(k, 0)
                    continue
                else:
                    raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            else:

                ## first_term_dict[k] (i.e. A[k]) is zero

                ##                                      ((prev set))
                ##           A[k+1]  ----> B[k+1] ------> C[k+1]
                ##           0       ----> B[k]   ------> (sum of surrounding terms)
                ##           A[k-1]  ----> B[k-1] ------>    ?                   
                ##

                if middle_term_dict.get(k-1, 0) == 0:
                    ## we simply have a short-exact sequence
                    return_dict[k] = middle_term_dict.get(k, 0) + first_term_dict.get(k-1, 0)
                    continue

                elif first_term_dict.get(k-1, 0) == 0:
                    return_dict[k] = middle_term_dict.get(k, 0) 
                    continue

                else: 
                    ## A[k-1] and B[k-1] are both non-zero, so we cannot resolve the sequence
                    raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
        else:

            ## middle_term_dict[k] (i.e. B[k]) is zero

            ##                                      ((prev set))
            ##           A[k+1]  ----> B[k+1] ------> C[k+1]
            ##           A[k]    ----> 0   ------> (sum of surrounding terms)
            ##           A[k-1]  ----> B[k-1] ------>    ?                   
            ##


            #####
            # TODO: This logic actually doesnt work for an exact sequence of >3 terms
            #####
            if first_term_dict.get(k-1, 0) >= middle_term_dict.get(k-1, 0):
                return_dict[k] = first_term_dict.get(k-1, 0) - middle_term_dict.get(k-1, 0)
                continue
            else:
                return_dict[k] = 0
                    
            

    # Remove zero entries before returning
    return {k: v for k, v in return_dict.items() if v != 0}





def _compute_rhom_sph_to_line_bundle_helper(sph : SphericalTwistComposition, lb : LineBundle) -> Dict[int, int]:
    r"""!
    Recursive helper method which implements the general homological-algebraic logic for computing the right-derived Hom space
    of a spherical twist composition with a line bundle. The base case logic is handled by the general RHom wrapper. This method
    primarily employs long-exact sequence logic to try to compute the dimensions of the RHom space when possible. Ultimately, the 
    method is much less robust than one would hope for — even when short exact sequences exist in the long exact cohomology sequence,
    they sometimes cause LongExactSequenceExceptions due to the fact that we mathematically cant assume 0 -> ? -> B -> C -> ? -> 0 
    gives any information about the outside two terms.

    \param SphericalTwistComposition sph The first argument of RHom(-, -), which is assumed to be a SphericalTwistComposition for this helper method
    \param LineBundle lb The second argument of RHom(-, -), which is assumed to be a LineBundle for this helper method

    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
    
    \throws LongExactSequenceException If we encounter a sequence of non-zero terms in our Long exact sequence that do not yield a definitive answer for the dimension of our desired RHom space
    """


    ################################################################################################################
    #                                                                                                              #       
    #                      Recursive Calls from applying RHom(-----, O(a_n)) to                                    #
    #                                                                                                              #
    #      O(a_{n-1}) ---> Tw_{a_{n-2}}... Tw_{a_1} O(a_0) ----> Tw_{a_{n-1}} Tw_{a_{n-2}}... Tw_{a_1} O(a_0)      #
    #                                                                                                              #
    ################################################################################################################

    
    middle_term_dict = RHom(sph.defining_triangle[1], lb)

    third_term_dict = RHom(sph.defining_triangle[0], lb)


    ###########################################################################
    #                                                                         #
    #               Attempt to resolve long-exact sequence                    #
    #                                                                         #
    ###########################################################################



    # Combine first_term_dict and middle_term_dict into return_dict
    return_dict = {}

    raw_keys = set(third_term_dict) | set(middle_term_dict)

    if raw_keys:  # guard against empty dicts
        extended_keys = raw_keys | {min(raw_keys) - 1, max(raw_keys) + 1}
    else:
        extended_keys = set()

    keys = sorted(extended_keys, reverse=True)

    for k in keys:
        ## At each iteration, we must check that we are not overwriting a previous value
        ## in returned dictionary. Since there are only three terms in a triangle contributing
        ## to the long exact sequence, the only possible way that a value could be overwritten is
        ## if there is a long-exact sequence with at least 4 consecutive non-zero terms.
        ##
        ## (Resolving such an error would require more specialized care outside of pure diagram chasing;
        ## in particular, one would need to know the explicit maps in the long exact sequence.)
        ##
        ## Algorithmically, this is simply handled by saving the index we wish to overwrite to and 
        ## the value stored, and then checking to make sure that index is not already in the dictionary.

        long_ex_str = f"""
                [{k+1}] \t:\t {return_dict.get(k+1,0)} ---> {middle_term_dict.get(k+1,0)} ---> {third_term_dict.get(k+1,0)}
                [{k}] \t:\t ? ---> {middle_term_dict.get(k,0)} ---> {third_term_dict.get(k,0)}
                [{k-1}] \t:\t ? ---> {middle_term_dict.get(k-1,0)} ---> {third_term_dict.get(k-1,0)}
                """

        if middle_term_dict.get(k, 0) != 0:

            ##        ((prev set))
            ##           A[k+1]            ----> B[k+1] ------> C[k+1]
            ##  (sum of surrounding terms) ----> B[k]   ------> C[k]
            ##             ?               ----> B[k-1] ------> C[k-1]
            ##

            if third_term_dict.get(k, 0) != 0:
                ## ? ---> B[k+1] ---> C[k+1] ---> ? cannot be resolved
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            else:
                ##
                ##                 ?                ----> B[k+1] ------> C[k+1]
                ##    (sum of surrounding terms) ---->   B[k]    ------> 0
                ##
                if return_dict.get(k+1, 0) != 0:
                    raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
                else:
                    ##
                    ##                 0                ----> B[k+1] ------> C[k+1]
                    ##    (sum of surrounding terms) ---->   B[k]    ------> 0
                    ##
                    ##    Can be resolved since exactness implies dim B[k+1] - dim B[k] = dim C[k+1] - dim A[k]

                    return_dict[k] = third_term_dict.get(k+1, 0) - middle_term_dict.get(k+1, 0) + middle_term_dict.get(k, 0)
                    continue

        else:
            ##        ((prev set))
            ##           A[k+1]            ----> B[k+1] ------> C[k+1]
            ##  (sum of surrounding terms) ----> 0   ------> C[k]
            ##             ?               ----> B[k-1] ------> C[k-1]
            ##
            ##   From the previous case, we know that one of B[k+1] and C[k+1] must be 0. If C[k+1] is 0, then
            ##   A[k] is surrounded by 0s and is thus 0. If B[k+1] is 0, then A[k] = C[k+1]. In either case,
            ##   we have A[k] = C[k+1] and thus the long-exact sequence is resolved.

            return_dict[k] = third_term_dict.get(k+1, 0) 


    # Remove zero entries before returning
    return {k: v for k, v in return_dict.items() if v != 0}


def _shift_dict(d: Dict[int, int], shift: int) -> Dict[int, int]:
    r"""!
    Helper function to shift the keys of a dictionary by a given amount. This is used to shift the
    degrees of the RHom space by some cohomological degree. 

    \param dict d The dictionary to shift
    \param int shift The amount to shift the keys by

    \return dict A new dictionary with the keys shifted by the given amount
    """

    return {k + shift: v for k, v in d.items()}

def _add_dicts(d1: Dict[int, int], d2: Dict[int, int]) -> Dict[int, int]:
    r"""!
    Helper function to combine two dictionaries by adding their values together. This is used to combine
    the dimensions of the RHom spaces in the long exact sequence.

    \param dict d1 The first dictionary to combine
    \param dict d2 The second dictionary to combine
    \return dict A new dictionary with the keys and values combined
    """
    result = d1.copy()
    for k, v in d2.items():
        result[k] = result.get(k, 0) + v
    return result

def _multiply_dict(d: Dict[int, int], scalar: int) -> Dict[int, int]:
    r"""!
    Helper function to multiply the values of a dictionary by a given scalar. This is used to scale
    the dimensions of the RHom spaces in the long exact sequence.

    \param dict d The dictionary to scale
    \param int scalar The scalar to multiply the values by

    \return dict A new dictionary with the values scaled by the given scalar
    """
    return {k: v * scalar for k, v in d.items()}