from src.DerivedCategory.SphericalTwist import SphericalTwistComposition
from src.DerivedCategory.CoherentSheaf import LineBundle, CoherentSheaf
from src.DerivedCategory.DistinguishedTriangle import DistinguishedTriangle
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject, ZeroObject

from typing import List, Dict
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
    Primary function for computing the RHom space between two derived category objects, as a graded C-vector space. The function primarily relies on the helper method _compute_rhom_helper to recursively
    compute the RHom space between shorter, less complicated series of twists. At each step roughly
    three recursive calls are made, only one of which is constant time since it defaults to the base case; thus, the
    expected complexity of this function is O(2^n) where n is the number of line bundles in the list.
    The function is designed to be called with a list of line bundles, which are assumed to be defined on the same
    variety; however, for some of the local Calabi-Yau varieties, the obstruction of the RHom just between two line bundles
    makes it impossible to compute the RHom for a higher number of twists using basic diagram chasing.

    \param object1 A DerivedCategoryObject representing the first object in the RHom space
    \param object2 A DerivedCategoryObject representing the second object in the RHom space

    
    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space

    \throws TypeError If object1 or object2 is not a DerivedCategoryObject
    \throws ValueError If object1 or object2 is not a LineBundle
    \throws ValueError If there are not at least two line bundles in the list
    \throws ValueError If the line bundles are not defined on the same variety
    \throws TypeError If degree_K3 is not an integer
    \throws ValueError If degree_K3 is not a positive integer
    """


    if isinstance(object1, LineBundle) and isinstance(object2, LineBundle):
        ## We are computing the RHom space between two line bundles; this
        ## should just be handled by the standard Ext group of the respective
        ## abelian category
        return Ext(object1, object2)
    
    elif isinstance(object1, LineBundle) and isinstance(object2, SphericalTwistComposition):
        return _compute_rhom_line_bundle_to_sph_helper(object1, object2)
    elif isinstance(object1, LineBundle) and isinstance(object2, GradedCoproductObject):
        return _compute_rhom_line_bundle_to_graded_coproduct_helper(object1, object2)
    elif isinstance(object1, GradedCoproductObject) and isinstance(object2, LineBundle):
        return _compute_rhom_graded_coproduct_to_line_bundle_helper(object1, object2)
    elif isinstance(object1, GradedCoproductObject) and isinstance(object2, GradedCoproductObject):
        return _compute_rhom_graded_coproduct_to_graded_coproduct_helper(object1, object2)
    





def Ext(object1 : CoherentSheaf, object2 : CoherentSheaf) -> Dict[int,int]:

    if isinstance(object1, LineBundle) and isinstance(object2, LineBundle):
        return _ext_line_bundles(object1, object2)
    else:
        # TODO: Figure out how to implement for more general coherent sheaves
        return NotImplementedError("Currently, Ext(A, B) is only implemented for line bundles. \n\t\t\t We hope to implement general coherent sheaves in the future.")






def _ext_line_bundles(line_bundle_1 : LineBundle, line_bundle_2 : LineBundle) -> Dict[int,int]:


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

        if isinstance(obj1, LineBundle):
            hom_dict = _compute_rhom_line_bundle_to_graded_coproduct_helper(obj1, gc2)
        else:
            raise NotImplementedError(f"RHom is not yet implemented for {type(obj1)} --> {type(gc2)}")






def _compute_rhom_graded_coproduct_to_line_bundle_helper(gc : GradedCoproductObject, lb : LineBundle) -> Dict[int, int]:

    r"""!
    Recursive helper method which implements the general homological-algebraic logic for 
    computing the right-derived Hom space of a graded coproduct with a line bundle.
    The base case is handled by the _dimHom_LineBundlesP1, _dimHom_LineBundlesP2, and
    _dimHom_LineBundlesK3 methods, which simply utilize Hirzebruch-Riemann-Roch and some easy
    diagram chasing.

    \param GradedCoproductObject gc The graded coproduct object to compute the RHom space for
    \param LineBundle lb The line bundle to compute the RHom space with

    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space

    """



    direct_sum_dict = {}
    for obj, shift, dim in zip(gc.object_vector, gc.shift_vector, gc.dimension_vector):

        hom_dict = RHom(obj, lb)

        temp_dict = _shift_dict(_multiply_dict(hom_dict, dim), -1*shift)
        direct_sum_dict = _add_dicts(direct_sum_dict, temp_dict)

    return direct_sum_dict








def _compute_rhom_line_bundle_to_graded_coproduct_helper(lb : LineBundle, gc : GradedCoproductObject) -> Dict[int, int]:


    direct_sum_dict = {}
    for obj, shift, dim in zip(gc.object_vector, gc.shift_vector, gc.dimension_vector):

        hom_dict = RHom(lb, obj)

        temp_dict = _shift_dict(_multiply_dict(hom_dict, dim), shift)
        direct_sum_dict = _add_dicts(direct_sum_dict, temp_dict)

    return direct_sum_dict
        




def _compute_rhom_line_bundle_to_sph_helper(lb : LineBundle, sph : SphericalTwistComposition) -> Dict[int, int]:
    r"""!
    Recursive helper method which implements the general homological-algebraic logic for 
    computing the right-derived Hom space

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

    keys = sorted(set(first_term_dict) | set(middle_term_dict), reverse=True)
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

        a = first_term_dict.get(k)
        b = middle_term_dict.get(k)

        if a is not None and b is None:
            ##
            ##         ----> B[k+1] ------> (sum of surrounding terms)
            ##    A[k] ---->   0    ------>
            ##


            if first_term_dict.get(k+1, 0) != 0 and middle_term_dict.get(k+1, 0) != 0:
                ## A[k+1] ---> B[k+1] ---> C[k+1] ---> A[k] cannot be resolved
                long_ex_str = f"\n\t[{k+1}]\t:\t {first_term_dict[k+1]} ---> {middle_term_dict[k+1]} ---> ?\n\t[{k}]\t:\t {first_term_dict[k]} ---> 0"
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            
            target_k = k + 1 
            value = a + middle_term_dict.get(k + 1, 0)

        elif a is None and b is not None:
            ##
            ##     0    ----> B[k] ------> (sum of surrounding terms)
            ##  A[k-1]  ----> 
            ##

            if first_term_dict.get(k-1, 0) != 0 and middle_term_dict.get(k-1, 0) != 0:
                ## B[k] ---> C[k] ---> A[k-1] ---> B[k-1] cannot be resolved
                long_ex_str = f"\n\t[{k}]\t:\t 0 ---> {middle_term_dict[k]} ---> ?\n\t[{k-1}]\t:\t {first_term_dict[k-1]} ---> {middle_term_dict[k-1]}"
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)

            target_k = k
            value = b + first_term_dict.get(k - 1, 0)
        elif a is not None and b is not None:
            ## Theoretically it is allowed for 
            ##     -----> 0 ------> C[k+1] ----->
            ## A[k]-----> B[k] -----> 0
            ##
            ## But we must be careful to check which 
            if return_dict.get(k + 1, 0) != 0 and middle_term_dict.get(k+1, 0) != 0:
                long_ex_str = f"\n\t[{k+1}]\t:\t        ---->{middle_term_dict [k + 1]} ---> {return_dict[k + 1]}\n\t[{k}]\t:\t {first_term_dict[k]} ---> {middle_term_dict[k]} ---> ?"
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            elif return_dict.get(k+1, 0) != 0 and a != b + return_dict[k+1]:
                long_ex_str = f"\n\t[{k+1}]\t:\t                      ---> {return_dict[k + 1]}\n\t[{k}]\t:\t {first_term_dict[k]} ---> {middle_term_dict[k]} ---> ?"
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            elif return_dict.get(k+1, 0) != 0:
                return_dict[k] = 0
                continue
            elif return_dict.get(k+1, 0) == 0 and first_term_dict.get(k-1, 0) != 0:
                long_ex_str = f"\n\t[{k+1}]\t:\t                      ---> 0\n\t[{k}]\t:\t {first_term_dict[k]} ---> {middle_term_dict[k]} ---> ?\n\t[{k-1}]\t:\t {first_term_dict[k-1]} ---> "
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            else:
                if b < a:
                    raise LongExactSequenceException("Encountered sequence 0 -> A[k] -> B[k] -> C[k] -> 0  with dim B < dim A")

                target_k = k
                value = b - a
                
            
        else:
            ## Both entries on this line are none; while the result may still be non-zero for this degree, it will be resolved
            ## on another line
            target_k = k

            if first_term_dict.get(k-1, 0) >= middle_term_dict.get(k-1, 0):
                value = first_term_dict.get(k-1, 0) - middle_term_dict.get(k-1, 0)
            else:
                # next map should be an injection
                value = 0

        if target_k in return_dict:
            long_ex_str = f"\n\t[{target_k-1}]\t:\t {first_term_dict.get(target_k-1,0)} ---> {middle_term_dict.get(target_k-1,0)} ---> {return_dict.get(target_k-1,0)}\n\t[{target_k}]\t:\t {first_term_dict.get(target_k,0)} ---> {middle_term_dict.get(target_k,0)} ---> {return_dict.get(target_k,0)}\n\t[{target_k+1}]\t:\t {first_term_dict.get(target_k+1,0)} ---> {middle_term_dict.get(target_k+1,0)} ---> {return_dict.get(target_k+1,0)}"
            raise LongExactSequenceException(f"Overwriting return_dict[{target_k}] --- a long exact sequence was not caught", sequence_str=long_ex_str)
        ## We are guaranteed that this is the first time we are writing to this index
        return_dict[target_k] = value

    return return_dict





def _compute_rhom_sph_to_lb_helper(sph : SphericalTwistComposition, lb : LineBundle) -> Dict[int, int]:
    r"""!
    Recursive helper method which implements the general homological-algebraic logic for 
    computing the right-derived Hom space of a
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

    keys = sorted(set(third_term_dict) | set(middle_term_dict), reverse=True)
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

        b = middle_term_dict.get(k)
        c = third_term_dict.get(k)

        if b is not None and c is None:
            ##
            ##                                 ----> B[k+1] ------> C[k+1]
            ##    (sum of surrounding terms) ---->   B[k]    ------> 0
            ##


            if middle_term_dict.get(k+1, 0) != 0 and third_term_dict.get(k+1, 0) != 0:
                ## A[k+1] ---> B[k+1] ---> C[k+1] ---> A[k] cannot be resolved
                long_ex_str = f"\n\t[{k+1}]\t:\t ? ---> {middle_term_dict[k+1]} ---> {third_term_dict.get(k+1,0)}\n\t[{k}]\t:\t ? ---> {middle_term_dict[k]} ---> 0"
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            
            target_k = k
            value = b + third_term_dict.get(k + 1, 0)

        elif b is None and c is not None:
            ##
            ##                             ---->   0    ------> C[k]
            ##  (sum of surrounding terms) ----> B[k-1] ------> C[k-1]
            ##

            if third_term_dict.get(k-1, 0) != 0 and middle_term_dict.get(k-1, 0) != 0:
                ## B[k] ---> C[k] ---> A[k-1] ---> B[k-1] cannot be resolved
                long_ex_str = f"\n\t[{k}]\t:\t     ---> 0 ---> {c}\n\t[{k-1}]\t:\t? ---> {middle_term_dict[k-1]} ---> {third_term_dict[k-1]}"   
                raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)

            target_k = k-1
            value = c + middle_term_dict.get(k - 1, 0)
        elif c is not None and b is not None:
            ## Both terms on this line are non-zero; we must be careful to check which dimension is larger so
            ## that we are obeying the dimensionality constraints of exactness
            if return_dict.get(k + 1, 0) != 0:
                    long_ex_str = f"\n\t[{k+1}]\t:\t              ---> {return_dict[k + 1]}\n\t[{k}]\t:\t {first_term_dict[k]} ---> {middle_term_dict[k]} ---> ?"
                    raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
            
            if b >= a:
                if first_term_dict.get(k-1, 0) != 0:
                    long_ex_str = f"\n\t[{k}]\t:\t {first_term_dict[k]} ---> {middle_term_dict[k]} ---> ?\n\t[{k-1}]\t:\t {first_term_dict[k-1]} ---> "
                    raise LongExactSequenceException("Cannot resolve long-exact sequence", sequence_str=long_ex_str)
                target_k = k
                value = b - a
            else:
                ## this is an edge case when we are computing spherical twists with consecutive degrees, e.g. Tw_1 Tw_2 O(3)
                ## or Tw_3 Tw_2 O(1). It generally does not occur
                target_k = k - 1
                value = a - b
        else:
            ## Both entries on this line are none; while the result may still be non-zero for this degree, it will be resolved
            ## on another line
            continue

        if target_k in return_dict:
            long_ex_str = f"\n\t[{target_k-1}]\t:\t {first_term_dict.get(target_k-1,0)} ---> {middle_term_dict.get(target_k-1,0)} ---> {return_dict.get(target_k-1,0)}\n\t[{target_k}]\t:\t {first_term_dict.get(target_k,0)} ---> {middle_term_dict.get(target_k,0)} ---> {return_dict.get(target_k,0)}\n\t[{target_k+1}]\t:\t {first_term_dict.get(target_k+1,0)} ---> {middle_term_dict.get(target_k+1,0)} ---> {return_dict.get(target_k+1,0)}"
            raise LongExactSequenceException(f"Overwriting return_dict[{target_k}] --- a long exact sequence was not caught", sequence_str=long_ex_str)
        ## We are guaranteed that this is the first time we are writing to this index
        return_dict[target_k] = value

    return return_dict


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