from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.DistinguishedTriangle import DistinguishedTriangle
from src.DerivedCategory.GeometryContext import GeometryContext


import math
import cmath
import json
from typing import List, Dict
from functools import cached_property


from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']








#                           Spherical Twist                                   #
# ----------------------------------------------------------------------------#
#  This object is used to represent the composition of spherical twists in    #
#  in the derived category of local projective space  consisting of sheaves   #
#  supported on the zero divisor (the original projective space)              #
#  In particular, the only compositions of twists we consier are twists       #
#  around (pushforwards) of line bundles. Theoretically, since the derived    #
#  category of coherent sheaves on projective space is constructible, the     #
#  braid relations between spherical twists allow any object obtained as a    #
#  series of spherical twists to be represented as a sequence of spherical    #
#  twists specifically around line bundles.                                   #
#                                                                             #
#  In order to determine possible Harder-Narasimhan filtrations of the        #
#  spherical twists, we need to iteratively keep track of previous Harder-    #
#  Narasimhan filtrations. At the moment, this is only computed up to two     #
#  successive spherical twists.                                               #






class LongExactSequenceException(Exception):
    r"""!
    Exception raised when the long exact sequence cannot be computed
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args)

        self.message = kwargs.get('message') ## The error message

        self.sequence_str = kwargs.get('sequence_str') ## The string representation of the sequence







class SphericalTwistComposition(DerivedCategoryObject):
    r"""!
    A spherical twist is a non-standard autoequivalence of the derived category of coherent sheaves, in the 
    sense that it does not arise from any composition of (1) standard autoequivalences on the variety (2) 
    tensoring by line bundles and (3) (co)homological shifts. Such autoequivalences typically only arise in the
    case of Calabi-Yau categories and toric varieties (where (-2)-curves can exist in the Fano setting as well, 
    e.g. P^2 blown up at 2 points), and often control the structure of the stability manifold. They are also 
    relevant to homological mirror symmetry since they are the derived equivalent of Dehn twists in the Fukaya
    category (i.e. the symplectic setting). They are explicitly defined as the cone of the derived evaluation morphism

            Hom(A, B) ⊗ A ---->  B ----> Tw_A B

    where A is a spherical object in the sense that RHom(A,A) is isomorphic as a graded-vector space to the 
    singular cohomology of an n-sphere.


    Currently, we only consider spherical twists around line bundles in the derived category of coherent sheaves
    since they always yield examples of spherical objects for Local P^n and K3 surfaces. On K3 surfaces, it is not
    true that the spherical twists account for all spherical objects, but they are still provide a rich source of
    examples to help predict mass growth.
    """

    _instances = {} ## List of all instances of the SphericalTwist class

    _canonical_triangle_cache = {} ## Cache of the canonical triangle for the spherical twist


    def __new__(cls, line_bundle_vector : List[LineBundle], geometry_context : GeometryContext):

        key = (line_bundle_vector, geometry_context)
        if key not in cls._instances:
            instance = super().__new__(cls)
            cls._instances[key] = instance

        return cls._instances[key]

    
    def __init__(self, line_bundle_vector : List[LineBundle], geometry_context : GeometryContext):
        r"""!
        Initialize an instance of SphericalTwist with the specified line bundles. The spherical twist
        is defined as the cone of the evaluation morphism 

                Hom(O(a), O(b)) ⊗ O(a) ---->  O(b) ----> Tw_a O(b)

        where O(a) and O(b) denote either line bundles or pushforwards of line bundles along the inclusion of the zero section. 
        The spherical twist is represented as a distinguished triangle in the derived category of coherent
        sheaves. 

        Several helper methods are used to compute the dimensions of the Hom spaces between the pushforwards
        of the line bundles, and then to construct the distinguished triangle.

        """

        if hasattr(self, '_initialized'):
            return

        #######
        # Parameter checks
        #######

        if not line_bundle_vector:
            raise ValueError("line_bundle_vector cannot be empty.")
        elif len(line_bundle_vector) < 2:
            raise ValueError("line_bundle_vector must contain at least two line bundles.")
        
        if not isinstance(geometry_context, GeometryContext):
            raise TypeError("geometry_context must be an instance of GeometryContext.")
        
        if not all(isinstance(obj, LineBundle) for obj in line_bundle_vector):
            raise TypeError("All elements of line_bundle_vector must be instances of LineBundle.")
        
        if not all(lb.catagory() ==  geometry_context.catagory for lb in line_bundle_vector):
            raise TypeError("All elements of line_bundle_vector must arise from the same underlying category")

        #######
        # Set remaining member variables
        #######
        
        self._initialized = True ## Flag to indicate that the object has been initialized

        self.line_bundle_vector = line_bundle_vector ## The vector of line bundles in the Hom space

        self.geometry_context = geometry_context ## The geometric context of the spherical twist


    def __str__(self):
        r"""!
        Returns a string representation of the spherical twist via its line bundle degrees

        \return str A string representation of the spherical twist
        """
        ret_str = ""

        for line_bundle in self.line_bundle_vector[:0:-1]:
            ret_str += f"Tw_{line_bundle} "
        
        ret_str += str(self.line_bundle_vector[0])
        return ret_str
    
    @cached_property
    def defining_triangle(self) -> DistinguishedTriangle:

        r"""!
        Method to compute the defining triangle for the spherical twist. The defining triangle is given by

         RHom(O(a_n), Tw_O(a_n-1)...Tw_O(a_1) O(a_0)) ⊗ O(a_n) ----> Tw_O(a_n-1)...Tw_O(a_1) O(a_0) ----> Tw_O(a_n) Tw_O(a_n-1)...Tw_O(a_1) O(a_0)

        where O(a_0), ..., O(a_n) are the line bundles in the Hom space.

        \return DistinguishedTriangle The distinguished triangle representing the spherical twist
        """

        RHom_dict_first_line_bundle = RHom(self.line_bundle_vector, geometric_context=self.geometry_context)
        # Convert the RHom_dict into a CoherentSheafCoproduct object
        shift_vector, dimension_vector = zip(*RHom_dict_first_line_bundle.items())
        shift_vector = list(shift_vector)
        dimension_vector = list(dimension_vector)
        bundle_vector = [LineBundle(self.line_bundle_vector[0].divisor, self.geometry_context) ] * len(shift_vector) 

        first_triangle_object = LineBundleCoproduct(sheaf_vector=bundle_vector,
                                        shift_vector=shift_vector,
                                        dimension_vector=dimension_vector)
        
        second_triangle_object = None


        if len(self.line_bundle_vector) == 2:
            second_triangle_object = LineBundle(self.line_bundle_vector[1].divisor, self.geometry_context)
        else:
            second_triangle_object = SphericalTwistComposition(line_bundle_vector=self.line_bundle_vector[1:], geometry_context=self.geometry_context)

        return DistinguishedTriangle(derived_object1=first_triangle_object,
                                    derived_object2=second_triangle_object,
                                    derived_object3=self)
    

    # def defining_triangle_to_json(self):
    #     r"""!
    #     Method to convert the spherical twist to a JSON string. This is used to pass the chain complex
    #     data of the spherical twist in the browser. Since a spherical twist is defined by its defining
    #     triangle, 
        
    #             RHom(O(a), O(b)) ⊗ O(a) ----> O(b) ----> Tw_O(a) O(b)

    #     and Tw_O(a) O(b) is merely a symbol, the only data we really need to pass is what the first RHom
    #     object looks like (dimension, shifts, etc) and the degrees [a,b].

    #     \return str A JSON string representation of the chain complex. The first object contains the information
    #                 of the first RHom object, and the degrees of the line bundles are stored in the degrees key.
    #     """

    #     object1 = {
    #         "shift_vector" : self.defining_triangle.object1.shift_vector,
    #         "dimension_vector" : self.defining_triangle.object1.dimension_vector
    #     }

    #     chain_complex_data = {
    #         "object1" : object1,
    #         "degrees" : [self.line_bundle_1.degree,
    #                     self.line_bundle_2.degree]
    #     }

    #     return json.dumps(chain_complex_data)
    

    

    @cached_property
    def chernCharacter(self):
        r"""!
        Method to compute the Chern Character of the spherical twist. The Chern Character of the
        spherical twist is the Chern Character of the third object in the distinguished triangle.

        \return ChernCharacter The Chern Character of the spherical twist
        """

        defining_triangle = self.defining_triangle()

        return defining_triangle.object2.chernCharacter() - defining_triangle.object1.chernCharacter()
    
    def shift(self, n : int):
        r"""!
        Method to shift the spherical twist by n units. As a spherical twist is initially only
        specified as a string until its defining triangle is computed, the shift method simply
        relies on the implementation in the parent class DerivedCategoryObject.

        \param int n: The number of units to shift the object by

        \return DerivedCategoryObject The shifted object
        """

        if not isinstance(n, int):
            raise TypeError("Shift must be an integer")

        return GradedCoproductObject(sph_twists_vector=[self], dimension_vector=[1], shift_vector=[n], degree=self.degree)
    

    

                
                
                
                












###################################################################
#                  RHom Computation                               #
###################################################################



def ApplySphericalTwist(target, line_bundle : LineBundle, geometry_context : GeometryContext) -> DerivedCategoryObject:
    r"""!
    Helper function to apply a spherical twist to an object in the derived category. 


    \return DerivedCategoryObject The target after applying the spherical twist
    """

    if not isinstance(target, DerivedCategoryObject) and not isinstance(target, DistinguishedTriangle):
        raise TypeError("object must be an instance of DerivedCategoryObject or DistinguishedTriangle")

    if not isinstance(line_bundle, LineBundle):
        raise TypeError("twist must be an instance of LineBundle")
    if target.catagory != line_bundle.catagory:
        raise ValueError("target and twist must be defined on the same variety")


    if isinstance(target, SphericalTwistComposition):
        
        new_list = target.line_bundle_vector.copy()
        new_list.append(line_bundle)

        return SphericalTwistComposition(line_bundle_vector=new_list, geometry_context=geometry_context)
    
    elif isinstance(target, LineBundle):
        return SphericalTwistComposition(line_bundle_vector=[target, line_bundle], geometry_context=geometry_context)
    
    elif isinstance(target, GradedCoproductObject):
        
        sph_twist_vector = []

        for graded_obj in target.object_vector:
            sph_twist_vector.append( ApplySphericalTwist(target=graded_obj, line_bundle=line_bundle, geometry_context=geometry_context) )
            

        return GradedCoproductObject(sph_twists_vector=sph_twist_vector,
                                dimension_vector=target.dimension_vector,
                                shift_vector=target.shift_vector, 
                                degree=degree_K3)
    
    
    elif isinstance(target, DistinguishedTriangle):
        return DistinguishedTriangle(derived_object1=ApplySphericalTwist(target=target.object1, line_bundle=line_bundle, degree_K3=degree_K3),
                                    derived_object2=ApplySphericalTwist(target=target.object2, line_bundle=line_bundle, degree_K3=degree_K3),
                                    derived_object3=ApplySphericalTwist(target=target.object3, line_bundle=line_bundle, degree_K3=degree_K3))
    
    else:
        raise NotImplementedError(f"ApplySphericalTwist not implemented for type {type(target)};\n\t\t\t(only implemented for DistinguishedTriangles, SphericalTwistComposition, LineBundle, and GradedCoproductObject where the underlying objects are SphericalTwistComposition or LineBundle)")

        




def RHom(line_bundles : List[LineBundle], geometry_context : GeometryContext) -> Dict[int,int]:
    r"""!
    Primary function for computing the RHom space between a line bundle and a composition of spherical twists, as a 
    graded C-vector space. The function primarily relies on the helper method _computer_rhom_helper to recursively
    compute the RHom space between shorter, less complicated series of twists. At each step roughly
    three recursive calls are made, only one of which is constant time since it defaults to the base case; thus, the
    expected complexity of this function is O(2^n) where n is the number of line bundles in the list.
    The function is designed to be called with a list of line bundles, which are assumed to be defined on the same
    variety; however, for some of the local Calabi-Yau varieties, the obstruction of the RHom just between two line bundles
    makes it impossible to compute the RHom for a higher number of twists using basic diagram chasing.

    \param list line_bundles A list of LineBundles, which are assumed to be defined on the same variety

    
    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space

    \throws TypeError If line_bundles is not a list of LineBundles
    \throws ValueError If there are not at least two line bundles in the list
    \throws ValueError If the line bundles are not defined on the same variety
    \throws TypeError If degree_K3 is not an integer
    \throws ValueError If degree_K3 is not a positive integer
    """


    # Input validation
    if not isinstance(line_bundles, list):
        raise TypeError("line_bundles must be a list of LineBundles.")
    if not all(isinstance(x, LineBundle) for x in line_bundles):
        raise TypeError("line_bundles must be a list of LineBundles.")
    
    if len(line_bundles) < 2:
        raise ValueError("It does not make sense to compute RHom with less than 1 object; the sequence must have at least 2 elements")
    

    if not isinstance(geometry_context, GeometryContext):
        raise TypeError("geometry_context must be an instance of GeometryContext")

    if not all(lb.catagory() == geometry_context.catagory for lb in line_bundles):
        raise ValueError("Line bundles must be defined on the same underlying catagory")
    

    

    ## Convert the list of line bundles to a list of integers for our helper method
    line_bundle_degrees = [x.degree for x in line_bundles]

    return _compute_rhom_helper(line_bundle_degrees, catagory, degree_K3)





def _get_canonical_triangles_helper(line_bundles : List[LineBundle], degree_K3 : int = 1) -> List[DistinguishedTriangle]:

    if len(line_bundles) < 2:
        raise ValueError("It does not make sense to compute canonical triangles with less than 1 object; the sequence must have at least 2 elements")

    elif len(line_bundles) == 2:
        return [ SphericalTwistComposition(line_bundle_vector=line_bundles, degree_K3=degree_K3).defining_triangle() ]
    
    else:

        key = (line_bundles, degree_K3, line_bundles[0].catagory)

        if key not in SphericalTwistComposition._canonical_triangle_cache:

            previous_step_list = _get_canonical_triangles_helper(line_bundles[:-1], degree_K3=degree_K3)
            new_list = [SphericalTwistComposition(line_bundle_vector=line_bundles, degree_K3=degree_K3).defining_triangle()]

            for triangle in previous_step_list:
                new_list.append( ApplySphericalTwist(line_bundle_vector=line_bundles[:-1], degree_K3=degree_K3) )

            SphericalTwistComposition._canonical_triangle_cache[key] = new_list

        return SphericalTwistComposition._canonical_triangle_cache[key]






def _dimHom_LineBundlesP1(line_bundle_1 : int, line_bundle_2 : int) -> Dict[int,int]:
        r"""!
        Helper method which computes the dimension of the hom spaces between the pushforwards of the
        line bundles O(a) and O(b). The dimensions of the pushforwards are computed using the triangle

        i^* i_* E -> E -> E x O(2)[2]

        and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        the Hom spaces between line bundles on P1, the computation reduces to a case-by-case
        combinatorial problem. Since the homological index of the hom-space on P1 is bounded between
        0 and 1, the hom-space for local P1 is concentrated between degrees 0 and 2. Thus, we return
        a tuple of the form (a,b,c)

        \param int line_bundle_1 The degree of the first line bundle in the Hom space
        \param int line_bundle_2 The degree of the second line bundle in the Hom space

        \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
        """

        degree_dif = line_bundle_2 - line_bundle_1

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



def _dimHom_LineBundlesP2(line_bundle_1 : int, line_bundle_2 : int) -> Dict[int,int]:
        r"""!
        Helper method which computes the dimension of the hom spaces between the pushforwards of the
        line bundles O(a) and O(b). The dimensions of the pushforwards are computed using the triangle

        i^* i_* E -> E -> E x O(3)[2]

        and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        the Hom spaces between line bundles on P^2, the computation reduces to a case-by-case
        combinatorial problem. Since the homological index of the hom-space on P^2 is bounded between
        0 and 2, the hom-space for local P2 is concentrated between degrees 0 and 3. Thus, we return
        a tuple of the form (a,b,c,d)

        \param int line_bundle_1 The first line bundle in the Hom space
        \param int line_bundle_2 The second line bundle in the Hom space

        \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
        """

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
        

        
def _dimHom_LineBundlesK3(line_bundle_1 : int, line_bundle_2 : int, degree_K3 : int = 1) -> Dict[int,int]:
    r"""!
    Helper method which computes the dimension of the hom spaces between line bundles on a 
    degree d K3 surface. Since the only K3s we consider are Picard rank 1, if O(a) and O(b) are
    distinct then either O(b-a) or O(a-b) must be ample; by Serre duality and a result of Huybrechts,
    one may always argue that Ext1(O(a), O(b)) = 0. In fact, Serre duality implies that the complex must
    be concentrated in a single degree, which is either 0 or 2 and corresponds to the cases that b > a and
    a < b, respectively. 

    We return a tuple of the form (n0,n1,n2) indicating the dimension of the graded RHom space. Aside from the
    instance that a=b, this will only include one nonzero element. Hirzebruch-Riemann-Roch shows that

    dim RHom(O(a), O(b)) = 2 + d(b-a)^2 

    \param int line_bundle_1 The degree of the first line bundle in the Hom space
    \param int line_bundle_2 The degree of the second line bundle in the Hom space
    \param int degree_K3 The degree of the K3 surface; default is 1

    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space
    """

    degree_dif = line_bundle_2 - line_bundle_1

    if degree_dif == 0:
        return {0 : 1, -2 : 1}
    elif degree_dif > 0:
        return {0 : degree_K3 * degree_dif**2 + 2}
    else:
        return {-2: degree_K3 * degree_dif**2 + 2}
        




def _compute_rhom_helper(seq: List[int], geometry_context : GeometryContext) -> Dict[int, int]:
    r"""!
    Recursive helper method which implements the general homological-algebraic logic for 
    computing the right-derived Hom space of a line bundle with a successive number of twists.
    The base case is handled by the _dimHom_LineBundlesP1, _dimHom_LineBundlesP2, and
    _dimHom_LineBundlesK3 methods, which simply utilize Hirzebruch-Riemann-Roch and some easy
    diagram chasing.

    The order of the line bundles is the opposite of what would be traditionally read in English; if the
    function is called with [1, 2, 3, 4], then the output will represent the dimensions of the RHom space
    for RHom(O(4),  Tw_3 Tw_2 O(1) ). In particular, it should be noted that the last element of the sequence is
    not included in the actual spherical twist, but in fact determines the line bundle that we are mapping from.

    \param list seq The sequence of line bundles to compute the RHom space for

    

    \return dict A dictionary representing the dimensions of the RHom space, where the keys are the cohomological degrees and the values are the dimensions of the corresponding vector space

    \throws NotImplementedError If the catagory is not P1, P2, or K3
    \throws ValueError If computing the RHom of what are expected to be two line bundles is not concentrated in 2 degrees
    \throws LongExactSequenceError If the long-exact sequence cannot be resolved
    """

    


    ###########################################################################
    #                                                                         #
    #                            Base Case for Recursion                      #
    #                                                                         #
    ###########################################################################

    ## This is the base case we wish to start with. It is simply computing RHom(O(a), O(b)) for two line
    ## bundles. In the notation of this function, we have O(a) = O(seq[1]) and O(b) = O(seq[0])

    if len(seq) == 2:
        if geometry_context.catagory == 'P1':
            return _dimHom_LineBundlesP1(seq[1], seq[0])
        elif geometry_context.catagory == 'P2':
            return _dimHom_LineBundlesP2(seq[1], seq[0])
        elif geometry_context.catagory == 'K3':

            degree_K3 = 1

            if geometry_context.polarization is None:
                raise ValueError("K3 surfaces require a polarization to compute the RHom space")
            try:
                H = geometry_context.polarization
                degree_K3 = geometry_context.divisor_data.evaluate(H**2)
            except:
                raise ValueError("Could not evaluate the polarization ** 2 on the K3 surface")


            return _dimHom_LineBundlesK3(seq[1], seq[0], degree_K3=degree_K3)
        else:
            raise NotImplementedError("Only P1, P2 and K3 catagories are implemented")
        



    else:
        ## This is the recursive step. We start by computing the necessary sequences to compute the
        ## RHoms
    


        prev_step_seq = seq[:-1] ## This should be the sequence without the last element
        last_two_seq = [seq[-2], seq[-1]] ## This should be the last two elements of the sequence
        remove_2nd_last = seq[:-2] + [seq[-1]] ## This should be the sequence with the second to last element removed


        ################################################################################################################
        #                                                                                                              #       
        #                      Recursive Calls from applying RHom(O(a_n) , ----) to                                    #
        #                                                                                                              #
        #      O(a_{n-1}) ---> Tw_{a_{n-2}}... Tw_{a_1} O(a_0) ----> Tw_{a_{n-1}} Tw_{a_{n-2}}... Tw_{a_1} O(a_0)      #
        #                                                                                                              #
        ################################################################################################################


        ## The result of the previous step RHom( O(a_{n-1}) ,  Tw_{a_{n-1}} ... Tw_{a_1} O(a_0) ) will directly
        ## be applied to RHom( O(a_n) , O(a_{n-1}) )
        first_term_dict = _compute_rhom_helper(prev_step_seq, catagory, deg) 

        ## This is simply RHom( O(a_n) , O(a_{n-1}) ), which should be handled by the base case (e.g. this call is constant-time)
        case_dict = _compute_rhom_helper(last_two_seq, catagory, deg)

        ## This is a second call to recursion that does not complete in constant time. In particular, this means our function
        ## must make roughly 2^n calls to itself, where n is the number of elements in the sequence. The runtime can be improved
        ## via memoization, but this is not implemented here.
        middle_term_dict = _compute_rhom_helper(remove_2nd_last, catagory, deg)





        ################################################################################################################
        #                                                                                                              #       
        #                 Simplify RHom(O(a_n) , O(a_{n-1})) terms in first object of triangle using                   #
        #                                                                                                              #
        #   RHom( O(a_n),  O(a_{n-1})[b] ⊕ O(a_{n-1})[c] ) =  RHom( O(a_n),  O(a_{n-1}) )[b]  ⊕                        #
        #                                                                   RHom( O(a_n),  O(a_{n-1}) )[c]             #
        #                                                                                                              #
        ################################################################################################################    


        # Process first_term_dict depending on case_dict since RHom( - , - ) splits across direct sums and commutes with shifts
        direct_sum_dict = {}
        for shift_deg, dimension in case_dict.items():
            temp_dict = _shift_dict(_multiply_dict(first_term_dict, dimension), shift_deg)
            direct_sum_dict = _add_dicts(direct_sum_dict, temp_dict)

        first_term_dict = direct_sum_dict




        ###########################################################################
        #                                                                         #
        #               Attempt to resolve long-exact sequence                    #
        #                                                                         #
        ###########################################################################



        # Combine first_term_dict and middle_term_dict into return_dict
        return_dict = {}

        keys = set(first_term_dict) | set(middle_term_dict)
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


        






class DoubleSphericalTwist(DerivedCategoryObject):
    """!
    A class to represent the composition of successive spherical twists applied to a line bundle. The double
    spherical twist is given as a distinguished triangle similar to the case of the single spherical twist; however,
    for higher numbers of spherical twists, there are often multiple triangles that the object fits into. The added
    functionality that this class provides is the ability to account for both triangles when computing the Harder-
    Narasimhan filtration of the object. Specifically, one must account for the defining triangle

        RHom(O(a), Tw_b O(c)) ⊗ O(a) ---->  Tw_b O(c) ----> Tw_a Tw_b O(c)

    as well as what we refer to as the 'secondary canonical triangle' given by

        Tw_a (RHom(O(b), O(c)) ⊗ O(b)) ----> Tw_a O(c) ----> Tw_a Tw_b O(c)

    The Harder-Narasimhan factors of the double spherical twist are computed by first computing the Harder-Narasimhan
    factors of the defining triangle, and then the secondary canonical triangle. Unlike the single SphericalTwist class,
    we do not actually provide the Harder-Narasimhan filtration in all cases; there are edge cases where nothing can 
    currently be said and we must return an empty list leading to a mass of 0. 
    """
    

    def __init__(self, line_bundle_1, line_bundle_2, line_bundle_3, degree=1):
        r"""!
        Initialize an instance of DoubleSphericalTwist with the specified line bundles. The spherical twist
        is defined as the cone of the evaluation morphism 

                RHom(O(a), Tw_b O(c)) ⊗ O(a) ---->  Tw_b O(c) ----> Tw_a Tw_b O(c)

        where O(a), O(b), and O(c) denote line bundles. 
        The double spherical twist is represented as a distinguished triangle in the derived category of coherent
        sheaves. 

        Several helper methods are used to compute the dimensions of the RHom spaces between the pushforwards
        of the line bundles, and then to construct the distinguished triangle.

        \param LineBundleline_bundle_1 The last line bundle twisted around; i.e. O(a) where we are computing Tw_a Tw_b O(c)
        \param LineBundle line_bundle_2 The first line bundle twisted around; i.e. O(b) where we are computing Tw_a Tw_b O(c)
        \param LineBundle line_bundle_3 The line bundle we are applying the spherical twist to; i.e. O(c) where we are computing Tw_a Tw_b O(c)
        \param int degree An optional argument for the degree of the variety, which is relevant to the dimension of the
                       derived RHom space for K3 surfaces of picard rank 1. This does not affect the P1 or P2 implementations.
                       Default is 1.

        \throws TypeError If line_bundle_1, line_bundle_2, or line_bundle_3 are not instances of LineBundle
        \throws ValueError If the line bundles are not defined on the same catagory
        """

        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundleP1.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundleP1.")
        if not isinstance(line_bundle_3, LineBundle):
            raise TypeError("line_bundle_3 must be an instance of LineBundleP1.")
        
        if line_bundle_1.catagory != line_bundle_2.catagory or line_bundle_1.catagory != line_bundle_3.catagory:
            raise ValueError("Line bundles must be defined on the same variety")

        self.line_bundle_1 = line_bundle_1 ## The last line bundle twisted around; i.e. O(a) where we are computing Tw_a Tw_b O(c)

        self.line_bundle_2 = line_bundle_2 ## The first line bundle twisted around; i.e. O(b) where we are computing Tw_a Tw_b O(c)

        self.line_bundle_3 = line_bundle_3 ## The line bundle we are applying the spherical twist to; i.e. O(c) where we are computing Tw_a Tw_b O(c)

        self.catagory = line_bundle_1.catagory ## The catagory of the line bundles, e.g. P1, P2, K3
        
        self.degree = degree ## An optional argument for the degree of the K3 surface, if applicable

        self.defining_triangle = self._sph_twist_DoubleLineBundles(line_bundle_1, line_bundle_2, line_bundle_3) ## The distinguished triangle of the double spherical twist


        # UPDATE AS MORE CATAGORIES ARE IMPLEMENTED
        if self.catagory not in __CURRENT_DOUBLE_TWIST_IMPLEMENTED__:
            raise ValueError(f"Double spherical twists are currently only implemented for {__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ }")
        

    def _sph_twist_DoubleLineBundles(self, line_bundle_1, line_bundle_2, line_bundle_3):
        r"""!
        Helper method to compute the distinguished triangle of the double spherical twist. The distinguished triangle
        is given by the cone of the evaluation morphism 

            RHom(O(a), Tw_b O(c)) ⊗ O(a) ---->  Tw_b O(c) ----> Tw_a Tw_b O(c)

        \param LineBundle line_bundle_1 The last line bundle twisted around; i.e. O(a) where we are computing Tw_a Tw_b O(c)
        \param LineBundle line_bundle_2 The first line bundle twisted around; i.e. O(b) where we are computing Tw_a Tw_b O(c)
        \param LineBundle line_bundle_3 The line bundle we are applying the spherical twist to; i.e. O(c) where we are computing Tw_a Tw_b O(c)

        \return DistinguishedTriangle The distinguished triangle of the double spherical twist

        \throws TypeError If line_bundle_1, line_bundle_2, or line_bundle_3 are not instances of LineBundle
        """

        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundle.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundle.")
        if not isinstance(line_bundle_3, LineBundle):
            raise TypeError("line_bundle_3 must be an instance of LineBundle.")
            

        homDims = _dimHom_Line_and_SingleTwist(line_bundle_1, line_bundle_2, line_bundle_3, self.degree)  

        bundle_vector = []
        dimension_vector = [] 
        shift_vector = []

        # create the necessary lists for the ChainComplex constructor 
        for i in range(len(homDims)):
            if homDims[i] == 0:
                continue
            dimension_vector.append(homDims[i])
            # The homDims tuple represents the degrees (-1, 0, 1, 2, 3). This is represented by cohomological
            # shifts [1], [0], [-1], [-2], and [-3] respectively.
            shift_vector.append(1 - i)
            bundle_vector.append(LineBundle(line_bundle_1.degree, self.catagory))

        if not bundle_vector:
            print(f"line_bundle_1 = {line_bundle_1}\nline_bundle_2={line_bundle_2}\nline_bundle_3={line_bundle_3}\n_dimHom_Line_and_SingleTwist={homDims}")

        object1 = ChainComplex(sheaf_vector=bundle_vector, shift_vector=shift_vector, dimension_vector=dimension_vector)

        object2 = SphericalTwistCoproduct([(line_bundle_2, line_bundle_3)], dimension_vector=[1], shift_vector=[0], degree=self.degree)
        object3 = DerivedCategoryObject(string=f"Tw_{line_bundle_1.degree} Tw_{line_bundle_2.degree} O({line_bundle_3.degree})", catagory=self.catagory)

        return DistinguishedTriangle(object1, object2, object3)
    
    def chernCharacter(self):
        r"""!
        Method to compute the Chern character of the double spherical twist. The Chern character of the double
        spherical twist is the Chern character of the third object in the distinguished triangle.

        \return ChernCharacter The Chern Character of the double spherical twist
        """


        return self.defining_triangle.object3.chernCharacter()
    
    def central_charge(self, *args):
        r"""!
        Method to compute the central charge of the spherical twist. The central charge of the spherical
        twist is the central charge of the third object in the distinguished triangle.

        \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return complex The central charge of the spherical twist as a complex number

        \throws TypeError If the args are not of the correct type
        \throws ValueError If the number of args is incorrect
        \throws NotImplementedError If the catagory of the object is not P1, P2, or K3

        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
            
            ch = self.chernCharacter()
            
            return -1*ch[1] + args[0]*ch[0]

        # elif self.catagory == 'P2':
        #     if len(args) != 2:
        #         raise ValueError("Central charge of P2 requires two real number parameters")
        #     if not all(isinstance(x, (float, int)) for x in args):
        #         raise TypeError("P2 objects should have two real number parameters")
            
        #     ch = self.chernCharacter()
            
        #     return complex(-1*ch[2] +
        #                     args[1] * ch[0],
        #                       ch[1] - args[0] * ch[0])
            

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
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")

        
    
    def __str__(self):
        r"""!
        Returns a string representation of the spherical twist by printing the defining triangle

        \return str A string representation of the spherical twist
        """

        return str(self.defining_triangle.object3)
    

    def defining_triangle_to_json(self):
        r"""!
        Helper function to convert the data of the canonical triangle for a double spherical twist to a JSON string.
        The data includes the dimensions, shifts, and line bundles from the first object of the triangle (i.e. the SphericalTwistSum
        object), as well as a triple of the line bundle integers.

        \return str A JSON string representation of the spherical twist triangle data
        """

        object1 = {
            "shift_vector" : self.defining_triangle.object1.shift_vector,
            "dimension_vector" : self.defining_triangle.object1.dimension_vector
        }

        chain_complex_data = {
            "object1" : object1,
            "degrees" : [self.line_bundle_1.degree,
                        self.line_bundle_2.degree,
                        self.line_bundle_3.degree]
        }

        return json.dumps(chain_complex_data)  



    

    


    def secondary_triangle_to_json(self):
        r"""!
        Helper function to convert the data of the secondary canonical triangle to a JSON string. The data includes
        the first object of the secondary canonical trianle encoded as a triple of lists, as well as an integer triple
        consisting of the degrees of the three line bundles.

        \return str A JSON string representation of the spherical twist triangle data
        """

        secondary_canonical_triangle = self.secondary_canonical_triangle()

        object1 = {
            "shift_vector" : secondary_canonical_triangle.object1.shift_vector,
            "dimension_vector" : secondary_canonical_triangle.object1.dimension_vector
        }

        secondary_complex_data = {
            "object1" : object1,
            "degrees" : [self.line_bundle_1.degree,
                        self.line_bundle_2.degree,
                        self.line_bundle_3.degree]
        }

        return json.dumps(secondary_complex_data)    
            
            

        
    def is_semistable(self, *args, logging=False, log_file=None):
        r"""!
        Method to check if the double spherical twist is semistable. The double spherical twist is semistable
        if the Harder-Narasimhan factorization is trivial.

        \param tuple args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.
        \param bool logging A boolean flag to indicate whether to log the Harder-Narasimhan factors that caused the object to be unstable
        \param str log_file The file to log the Harder-Narasimhan factors that caused the object to be unstable

        \return bool True if the double spherical twist is semistable, False otherwise

        \throws TypeError If the args are not of the correct type
        \throws ValueError If the number of args is incorrect
        \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
        # elif self.catagory == 'P2':
        #     if len(args) != 2:
        #         raise ValueError("Central charge of P2 requires two real number parameters")
        #     if not all(isinstance(x, (float, int)) for x in args):
        #         raise TypeError("P2 objects should have two real number parameters")
        elif self.catagory == 'K3':
            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        

        try:
            return len(self.get_HN_factors(*args)) == 1
        except HarderNarasimhanError as e:
            
            if logging and log_file:
                with open(log_file, 'a') as log_file:
                    msg_str = e.message + f"@ {e.stability_parameters}"
                    log_file.write(msg_str)
            elif logging:
                msg_str = e.message + f"@ {e.stability_parameters}"
                print(msg_str)
            
            return False
    



        

       
        

    
    def get_HN_factors(self, *args):
        r"""!
        This is the crux of the DoubleSphericalTwist class, where we compute the Harder-Narasimhan factors of the
        double spherical twist. A signficiant assumption (CONJECTURAL) that we make is that the Harder-Narasimhan filtration
        must arise from the defining triangle or the secondary canonical triangle. This is not always the case, but we
        have not yet implemented a general method to compute the HN factors in all cases.

        This method works by first examining the two edge cases for the defining triangle: if the largest phase of the subobject
        is less than the smallest phase of the quotient, then we assume the object is stable (CONJECTURAL). If the smallest phase
        of the subobject is larger than the largest phase of the quotient, then we know for a fact (by BDL) that the Harder-Narasimhan
        filtration can be computed by concatenating the HN factors of the subobject and quotient. If neither of these cases hold, we move
        on to the second canonical triangle and apply a similar logic.

        \param tuple args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For
                        P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return list A list of tuples where the first element is a DerivedCategoryObject and the second element is a float
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
        # elif self.catagory == 'P2':
        #     if len(args) != 2:
        #         raise ValueError("Central charge of P2 requires two real number parameters")
        #     if not all(isinstance(x, (float, int)) for x in args):
        #         raise TypeError("P2 objects should have two real number parameters")
        elif self.catagory == 'K3':
            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        

        
        ###########
        # First check the defining triangle
        ###########

        # Write triangle as Tw_b O(c) -> Tw_a Tw b O(c) -> B + B[shift] + ...
        subobject = SphericalTwist(self.line_bundle_2, self.line_bundle_3, self.degree)

        modified_defining_triangle = self.defining_triangle.rotateLeft()
        quotient_complex = modified_defining_triangle.object3 


        if subobject.is_semistable(*args):
            # get the phase of the single twist
            left_side_phase = cmath.phase(subobject.central_charge(*args)) / math.pi
            right_side_min_phase = quotient_complex.get_smallest_phase(*args)
            right_side_max_phase = quotient_complex.get_largest_phase(*args)

            if left_side_phase <= right_side_min_phase:
                potential_phase = cmath.phase(self.central_charge(*args))/math.pi
                for n in range(-4,4):
                    if left_side_phase <= potential_phase + n and potential_phase + n <= right_side_min_phase:
                        return [(self, potential_phase + n)]
                    
                raise HarderNarasimhanError(message=f"{subobject} is semistable and its phase is smaller than {quotient_complex}, but cannot correctly find the phase", 
                                            stability_parameters=args)
                    

                
            elif left_side_phase > right_side_max_phase:
                # By BDL20, the HN factors of the subobject and quotient concatenate to make
                # the HN factors of the twist
                if len(quotient_complex.dimension_vector) == 1:
                    return [(subobject, left_side_phase),
                             (quotient_complex, right_side_max_phase)]
                else:
                    if len(quotient_complex) != 2:
                        raise ValueError("The Hom object is not concentrated in 1 or 2 degrees")

                    cplx_summand_0 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[0]],
                                                shift_vector=[quotient_complex.shift_vector[0]],
                                                    dimension_vector=[quotient_complex.dimension_vector[0]])
                    cplx_summand_1 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[1]],
                                                shift_vector=[quotient_complex.shift_vector[1]],
                                                    dimension_vector=[quotient_complex.dimension_vector[1]])
                    if right_side_min_phase == quotient_complex.sheaf_vector[0].phase(*args) + quotient_complex.shift_vector[0]:
                        return [(subobject, left_side_phase),
                                (cplx_summand_1, right_side_max_phase),
                                (cplx_summand_0, right_side_min_phase)]
                    else:
                        return [(subobject, left_side_phase),
                                (cplx_summand_0, right_side_max_phase),
                                (cplx_summand_1, right_side_min_phase)]
                
        else:
            # Subobject (i.e. Tw_b O(c)) is not semistable, so we must
            # first look at its HN filtration to see if anything can be salvaged
            HN_factors_subobject = None
            try:
                HN_factors_subobject = subobject.get_HN_factors(*args)
            except HarderNarasimhanError as e:
                ## Add context to the raised error message
                raise HarderNarasimhanError(message=f"Subobject {subobject} in defining triangle of {self} is not semistable, cannot find its HN factors:\n{e.message}", 
                                            stability_parameters=args)

            right_side_min_phase = quotient_complex.get_smallest_phase(*args)
            right_side_max_phase = quotient_complex.get_largest_phase(*args)

            if HN_factors_subobject[0][1] <= right_side_min_phase:
                # largest phase of subobject is less than smallest phase of quotient
                potential_phase = cmath.phase(self.central_charge(*args))/math.pi
                for n in range(-3,3):
                    if HN_factors_subobject[0][1] <= potential_phase + n and potential_phase + n <= right_side_min_phase:
                        return [(self, potential_phase + n)]
                    
            elif HN_factors_subobject[-1][1] > right_side_max_phase:
                # smallest phase of subobject is greater than largest phase of quotient
                if len(quotient_complex) == 1:
                    return HN_factors_subobject + [(quotient_complex, right_side_max_phase)]
                elif len(quotient_complex) != 2:
                    raise ValueError(f"The Hom object is not concentrated in 2 degrees; currently\n{quotient_complex}")

                cplx_summand_0 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[0]],
                                            shift_vector=[quotient_complex.shift_vector[0]],
                                                dimension_vector=[quotient_complex.dimension_vector[0]])
                cplx_summand_1 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[1]],
                                            shift_vector=[quotient_complex.shift_vector[1]],
                                                dimension_vector=[quotient_complex.dimension_vector[1]])
                if right_side_min_phase == quotient_complex.sheaf_vector[0].phase(*args) + quotient_complex.shift_vector[0]:
                    return HN_factors_subobject + [(cplx_summand_1, right_side_max_phase),
                                                   (cplx_summand_0, right_side_min_phase)]
                else:
                    return HN_factors_subobject + [(cplx_summand_0, right_side_max_phase),
                                                   (cplx_summand_1, right_side_min_phase)]

        

        ###########
        # next check secondary canonical triangle
        ###########

        secondary_canonical_triangle = self.secondary_canonical_triangle().rotateLeft()
        

        first_twist = SphericalTwist(secondary_canonical_triangle.object1.line_bundle_pairs_vector[0][0],
                                     secondary_canonical_triangle.object1.line_bundle_pairs_vector[0][1],
                                    self.degree)
        
        HN_factors_first_term = None
        try:
            HN_factors_first_term = first_twist.get_HN_factors(*args)
        except HarderNarasimhanError as e:
            raise HarderNarasimhanError(message=f"Subobject {first_twist} in secondary canonical triangle of {self} is not semistable, cannot find its HN factors:\n{e.message}", 
                                        stability_parameters=args)
        HN_factors_last_term = None
        try:
            HN_factors_last_term = secondary_canonical_triangle.object3.get_HN_factors_ordered(*args)
        except HarderNarasimhanError as e:
            raise HarderNarasimhanError(message=f"Quotient {secondary_canonical_triangle.object3} in secondary canonical triangle of {self} is not semistable, cannot find its HN factors:\n{e.message}", 
                                        stability_parameters=args)
        

        smallest_HN_phase_first = HN_factors_first_term[-1][1]
        largest_HN_phase_first = HN_factors_first_term[0][1]
        smallest_HN_phase_last = HN_factors_last_term[-1][1]
        largest_HN_phase_last = HN_factors_last_term[0][1]

        if largest_HN_phase_first <= smallest_HN_phase_last:
            potential_phase = cmath.phase(self.central_charge(*args))/math.pi
            for n in range(-3,3):
                if largest_HN_phase_first <= potential_phase + n and potential_phase + n <= right_side_min_phase:
                    return [(self, potential_phase + n)]
                
        elif smallest_HN_phase_first > largest_HN_phase_last:
            return HN_factors_first_term + HN_factors_last_term
        

        raise HarderNarasimhanError(message=f"The Harder-Narasimhan factors of both quotients and both subobjects do not match any of 4 known scenarios; their HN factors necessarily intertwine.",
                                     stability_parameters=args)
        









        



            

        


###################################################################
#                        Main                                     #
###################################################################


if __name__ == "__main__":
    # lb1 = LineBundle(-1, catagory='P2')
    # lb2 = LineBundle(0, catagory='P2')
    # lb3 = LineBundle(5, catagory='P2')

    # from LocalP2 import LePotier
    # DLP = LePotier(granularity=3, width=5)

    # sph = SphericalTwist(line_bundle_1=lb1, line_bundle_2=lb3)

    # x_vals = np.linspace(-5, 5, 150)  # X values from -2 to 2

    # # Generate y values satisfying y > x^2
    # y_vals = []
    # for x in x_vals:
    #     y_min =DLP.curve_estimate(x)
    #     y_max = 12
    #     y_range = np.linspace(y_min, y_max, 100)  # 50 points per x value
    #     y_vals.append(y_range)

    # # Convert to numpy array
    # y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # # Repeat x values to match the shape of y
    # x_vals = np.repeat(x_vals, 160)  # Each x value repeats 10 times

    # masses = np.array([sph.mass(x, y) for x, y in zip(x_vals, y_vals)])

    # # Plot the surface
    # fig = go.Figure(data=[go.Scatter3d(z=masses, x=x_vals, y=y_vals,
    #                                 mode='markers', marker=dict(size=3, color=masses, colorscale='viridis'))])

   
    # fig.show()



    # lb4 = LineBundle(-5, catagory='P1')
    # lb5 = LineBundle(2, catagory='P1')

    # sph2 = SphericalTwist(line_bundle_1=lb4, line_bundle_2=lb5, degree=1)

    # print(sph2)


    # lb6 = LineBundle(1, catagory='K3')
    # lb7 = LineBundle(2, catagory='K3')
    # lb8 = LineBundle(4, catagory='K3')

    # sph_sum = SphericalTwistSum([(lb6, lb7), (lb8, lb9), (lb6, lb9), (lb6, lb7)], [1, 5, 0, 10], [-1, 2, -5, -1], degree=1)
    # print(sph_sum)
    # print(sph_sum.chernCharacter())

    


    # sph3 = SphericalTwist(lb7, lb8, degree=2)




    # print(sph3.shift(3))
    # print("Is semistable: ", sph3.is_semistable(1, 2, 1))

    # K3_deg = 1
    

    # sph4 = DoubleSphericalTwist(lb6, lb7, lb8, degree=K3_deg)

    # # print(sph4.mass(2, 3, K3_deg))

    # x_vals = np.linspace(-5, 11.10, 200)  # X values from -2 to 2

    # # Generate y values satisfying y > x^2
    # y_vals = []
    # for x in x_vals:
    #     y_range = np.linspace(0, 5, 100)  # 50 points per x value
    #     y_vals.append(y_range)

    # # Convert to numpy array
    # y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # # Repeat x values to match the shape of y
    # x_vals = np.repeat(x_vals, 100)  # Each x value repeats 10 times

    # masses = np.array([sph4.mass(x, y, K3_deg) for x, y in zip(x_vals, y_vals)])

    # # Plot the surface
    # fig = go.Figure(data=[go.Scatter3d(z=masses, x=x_vals, y=y_vals,
    #                                 mode='markers', marker=dict(size=3, color=masses, colorscale='viridis'))])

    # fig.update_layout(scene = dict(xaxis = dict(nticks=4, range=[-5,5])))
   
    # fig.show()

    print(f"\n\n\nRHom(O(1), O(2)):\t {_compute_rhom_helper([2, 1], 'K3')}\n")
    print("--------------------------")
    print(f"This gives triangle:\n\n\t {_compute_rhom_helper([2, 1], 'K3')[0]}xO(1) ---> O(2) ---> Tw_1 O(2)\n")
    print(f"\nRHom(O(-1), O(1)):\t {_compute_rhom_helper([1, -1], 'K3')}\n")
    print(f"\nRHom(O(-1), O(2)):\t {_compute_rhom_helper([2, -1], 'K3')}\n")
    print(f"Putting everything together,\n\nRHom(O(-1), Tw_1 O(2): {_compute_rhom_helper([2, 1, -1] ,'K3')}\n")
    # print(_compute_rhom_helper([2, 1], 'K3'))