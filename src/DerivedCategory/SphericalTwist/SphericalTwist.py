from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.DistinguishedTriangle import DistinguishedTriangle


import json
from typing import List, Union
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

    def __new__(cls, line_bundle_vector : List[LineBundle]):

        key = tuple(line_bundle_vector)
        if key not in cls._instances:
            instance = super().__new__(cls)
            cls._instances[key] = instance

        return cls._instances[key]

    
    def __init__(self, line_bundle_vector : List[LineBundle]):
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
        
        if not all(isinstance(obj, LineBundle) for obj in line_bundle_vector):
            raise TypeError("All elements of line_bundle_vector must be instances of LineBundle.")
        

        first_geometry_context = line_bundle_vector[0].geometry_context
        if not all(lb.geometry_context == first_geometry_context for lb in line_bundle_vector):
            raise TypeError("All elements of line_bundle_vector must arise from the same underlying category")

        #######
        # Set remaining member variables
        #######
        
        self.line_bundle_vector = line_bundle_vector ## The vector of line bundles in the Hom space

        self.geometry_context = first_geometry_context ## The geometric context of the spherical twist

        self._initialized = True ## Flag to indicate that the object has been initialized



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

        last_line_bundle = LineBundle(self.line_bundle_vector[-1].divisor, self.geometry_context)

        second_triangle_object = None
        if len(self.line_bundle_vector) == 2:
            second_triangle_object = LineBundle(self.line_bundle_vector[0].divisor, self.geometry_context)
        else:
            second_triangle_object = SphericalTwistComposition(line_bundle_vector=self.line_bundle_vector[:-1])


        ## Delay import to avoid circular import issues
        from src.DerivedCategory.RHom import RHom

        ## Use RHom package to create dictionary of shifts -> dimensions
        RHom_dict_first_line_bundle = RHom(last_line_bundle, second_triangle_object)


        # First extract the keys (shifts) and values (dimensions) from the previous output
        shift_vector, dimension_vector = zip(*RHom_dict_first_line_bundle.items())
        shift_vector = list(shift_vector)
        dimension_vector = list(dimension_vector)


        ## Convert the dictionary data into a GradedCoproductObject
        first_triangle_object = GradedCoproductObject(object_vector=[last_line_bundle] * len(shift_vector) ,
                                        shift_vector=shift_vector,
                                        dimension_vector=dimension_vector)
        

        return DistinguishedTriangle(derived_object1=first_triangle_object,
                                    derived_object2=second_triangle_object,
                                    derived_object3=self)
    
    @cached_property
    def canonical_triangles(self) -> List[DistinguishedTriangle]:

        return _get_canonical_triangles_helper(self.line_bundle_vector)


    

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
    

    def __len__(self) -> int:

        return len(self.line_bundle_vector) - 1

    
    
    def chernCharacter(self):
        r"""!
        Method to compute the Chern Character of the spherical twist. The Chern Character of the
        spherical twist is the Chern Character of the third object in the distinguished triangle.

        \return ChernCharacter The Chern Character of the spherical twist
        """

        defining_triangle = self.defining_triangle

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

        return GradedCoproductObject(sph_twists_vector=[self], dimension_vector=[1], shift_vector=[n])
    

    











###################################################################
#                       Helper Methods                            #
###################################################################



def ApplySphericalTwist(
        target: Union[DerivedCategoryObject, DistinguishedTriangle],
        line_bundle : LineBundle
) -> DerivedCategoryObject:
    r"""!
    Helper function to apply a spherical twist to an object in the derived category. 


    \return DerivedCategoryObject The target after applying the spherical twist
    """

    if not isinstance(target, DerivedCategoryObject) and not isinstance(target, DistinguishedTriangle):
        raise TypeError("object must be an instance of DerivedCategoryObject or DistinguishedTriangle")

    if not isinstance(line_bundle, LineBundle):
        raise TypeError("twist must be an instance of LineBundle")
    
    if target.geometry_context != line_bundle.geometry_context:
        raise ValueError("target and line_bundle must be defined on the same variety / GeometricContext.")


    if isinstance(target, SphericalTwistComposition):
        
        new_list = target.line_bundle_vector.copy()
        new_list.append(line_bundle)

        return SphericalTwistComposition(line_bundle_vector=new_list)
    
    elif isinstance(target, LineBundle):
        return SphericalTwistComposition(line_bundle_vector=[target, line_bundle])
    
    elif isinstance(target, GradedCoproductObject):
        
        sph_twist_vector = []

        for graded_obj in target.object_vector:
            sph_twist_vector.append( ApplySphericalTwist(target=graded_obj, line_bundle=line_bundle) )
            

        return GradedCoproductObject(object_vector=sph_twist_vector,
                                dimension_vector=target.dimension_vector,
                                shift_vector=target.shift_vector)
    
    
    elif isinstance(target, DistinguishedTriangle):
        return DistinguishedTriangle(derived_object1=ApplySphericalTwist(target=target.object1, line_bundle=line_bundle),
                                    derived_object2=ApplySphericalTwist(target=target.object2, line_bundle=line_bundle),
                                    derived_object3=ApplySphericalTwist(target=target.object3, line_bundle=line_bundle))
    
    else:
        raise NotImplementedError(f"ApplySphericalTwist not implemented for type {type(target)};\n\t\t\t(only implemented for DistinguishedTriangles, SphericalTwistComposition, LineBundle, and GradedCoproductObject where the underlying objects are SphericalTwistComposition or LineBundle)")

        






def _get_canonical_triangles_helper(line_bundles : List[LineBundle]) -> List[DistinguishedTriangle]:

    if len(line_bundles) < 2:
        raise ValueError("It does not make sense to compute canonical triangles with less than 1 object; the sequence must have at least 2 elements")

    elif len(line_bundles) == 2:
        return [ SphericalTwistComposition(line_bundle_vector=line_bundles).defining_triangle ]
    
    else:

        key = tuple(line_bundles)

        if key not in SphericalTwistComposition._canonical_triangle_cache:

            previous_step_list = _get_canonical_triangles_helper(line_bundles[:-1] )
            new_list = [SphericalTwistComposition(line_bundle_vector=line_bundles).defining_triangle]

            for triangle in previous_step_list:
                new_list.append( ApplySphericalTwist(line_bundle=line_bundles[-1], target=triangle) )

            SphericalTwistComposition._canonical_triangle_cache[key] = new_list

        return SphericalTwistComposition._canonical_triangle_cache[key]

        






# class DoubleSphericalTwist(DerivedCategoryObject):
#     """!
#     A class to represent the composition of successive spherical twists applied to a line bundle. The double
#     spherical twist is given as a distinguished triangle similar to the case of the single spherical twist; however,
#     for higher numbers of spherical twists, there are often multiple triangles that the object fits into. The added
#     functionality that this class provides is the ability to account for both triangles when computing the Harder-
#     Narasimhan filtration of the object. Specifically, one must account for the defining triangle

#         RHom(O(a), Tw_b O(c)) ⊗ O(a) ---->  Tw_b O(c) ----> Tw_a Tw_b O(c)

#     as well as what we refer to as the 'secondary canonical triangle' given by

#         Tw_a (RHom(O(b), O(c)) ⊗ O(b)) ----> Tw_a O(c) ----> Tw_a Tw_b O(c)

#     The Harder-Narasimhan factors of the double spherical twist are computed by first computing the Harder-Narasimhan
#     factors of the defining triangle, and then the secondary canonical triangle. Unlike the single SphericalTwist class,
#     we do not actually provide the Harder-Narasimhan filtration in all cases; there are edge cases where nothing can 
#     currently be said and we must return an empty list leading to a mass of 0. 
#     """
    

#     def __init__(self, line_bundle_1, line_bundle_2, line_bundle_3, degree=1):
#         r"""!
#         Initialize an instance of DoubleSphericalTwist with the specified line bundles. The spherical twist
#         is defined as the cone of the evaluation morphism 

#                 RHom(O(a), Tw_b O(c)) ⊗ O(a) ---->  Tw_b O(c) ----> Tw_a Tw_b O(c)

#         where O(a), O(b), and O(c) denote line bundles. 
#         The double spherical twist is represented as a distinguished triangle in the derived category of coherent
#         sheaves. 

#         Several helper methods are used to compute the dimensions of the RHom spaces between the pushforwards
#         of the line bundles, and then to construct the distinguished triangle.

#         \param LineBundleline_bundle_1 The last line bundle twisted around; i.e. O(a) where we are computing Tw_a Tw_b O(c)
#         \param LineBundle line_bundle_2 The first line bundle twisted around; i.e. O(b) where we are computing Tw_a Tw_b O(c)
#         \param LineBundle line_bundle_3 The line bundle we are applying the spherical twist to; i.e. O(c) where we are computing Tw_a Tw_b O(c)
#         \param int degree An optional argument for the degree of the variety, which is relevant to the dimension of the
#                        derived RHom space for K3 surfaces of picard rank 1. This does not affect the P1 or P2 implementations.
#                        Default is 1.

#         \throws TypeError If line_bundle_1, line_bundle_2, or line_bundle_3 are not instances of LineBundle
#         \throws ValueError If the line bundles are not defined on the same catagory
#         """

#         if not isinstance(line_bundle_1, LineBundle):
#             raise TypeError("line_bundle_1 must be an instance of LineBundleP1.")
#         if not isinstance(line_bundle_2, LineBundle):
#             raise TypeError("line_bundle_2 must be an instance of LineBundleP1.")
#         if not isinstance(line_bundle_3, LineBundle):
#             raise TypeError("line_bundle_3 must be an instance of LineBundleP1.")
        
#         if line_bundle_1.catagory != line_bundle_2.catagory or line_bundle_1.catagory != line_bundle_3.catagory:
#             raise ValueError("Line bundles must be defined on the same variety")

#         self.line_bundle_1 = line_bundle_1 ## The last line bundle twisted around; i.e. O(a) where we are computing Tw_a Tw_b O(c)

#         self.line_bundle_2 = line_bundle_2 ## The first line bundle twisted around; i.e. O(b) where we are computing Tw_a Tw_b O(c)

#         self.line_bundle_3 = line_bundle_3 ## The line bundle we are applying the spherical twist to; i.e. O(c) where we are computing Tw_a Tw_b O(c)

#         self.catagory = line_bundle_1.catagory ## The catagory of the line bundles, e.g. P1, P2, K3
        
#         self.degree = degree ## An optional argument for the degree of the K3 surface, if applicable

#         self.defining_triangle = self._sph_twist_DoubleLineBundles(line_bundle_1, line_bundle_2, line_bundle_3) ## The distinguished triangle of the double spherical twist


#         # UPDATE AS MORE CATAGORIES ARE IMPLEMENTED
#         if self.catagory not in __CURRENT_DOUBLE_TWIST_IMPLEMENTED__:
#             raise ValueError(f"Double spherical twists are currently only implemented for {__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ }")
        

#     def _sph_twist_DoubleLineBundles(self, line_bundle_1, line_bundle_2, line_bundle_3):
#         r"""!
#         Helper method to compute the distinguished triangle of the double spherical twist. The distinguished triangle
#         is given by the cone of the evaluation morphism 

#             RHom(O(a), Tw_b O(c)) ⊗ O(a) ---->  Tw_b O(c) ----> Tw_a Tw_b O(c)

#         \param LineBundle line_bundle_1 The last line bundle twisted around; i.e. O(a) where we are computing Tw_a Tw_b O(c)
#         \param LineBundle line_bundle_2 The first line bundle twisted around; i.e. O(b) where we are computing Tw_a Tw_b O(c)
#         \param LineBundle line_bundle_3 The line bundle we are applying the spherical twist to; i.e. O(c) where we are computing Tw_a Tw_b O(c)

#         \return DistinguishedTriangle The distinguished triangle of the double spherical twist

#         \throws TypeError If line_bundle_1, line_bundle_2, or line_bundle_3 are not instances of LineBundle
#         """

#         if not isinstance(line_bundle_1, LineBundle):
#             raise TypeError("line_bundle_1 must be an instance of LineBundle.")
#         if not isinstance(line_bundle_2, LineBundle):
#             raise TypeError("line_bundle_2 must be an instance of LineBundle.")
#         if not isinstance(line_bundle_3, LineBundle):
#             raise TypeError("line_bundle_3 must be an instance of LineBundle.")
            

#         homDims = _dimHom_Line_and_SingleTwist(line_bundle_1, line_bundle_2, line_bundle_3, self.degree)  

#         bundle_vector = []
#         dimension_vector = [] 
#         shift_vector = []

#         # create the necessary lists for the ChainComplex constructor 
#         for i in range(len(homDims)):
#             if homDims[i] == 0:
#                 continue
#             dimension_vector.append(homDims[i])
#             # The homDims tuple represents the degrees (-1, 0, 1, 2, 3). This is represented by cohomological
#             # shifts [1], [0], [-1], [-2], and [-3] respectively.
#             shift_vector.append(1 - i)
#             bundle_vector.append(LineBundle(line_bundle_1.degree, self.catagory))

#         if not bundle_vector:
#             print(f"line_bundle_1 = {line_bundle_1}\nline_bundle_2={line_bundle_2}\nline_bundle_3={line_bundle_3}\n_dimHom_Line_and_SingleTwist={homDims}")

#         object1 = ChainComplex(sheaf_vector=bundle_vector, shift_vector=shift_vector, dimension_vector=dimension_vector)

#         object2 = SphericalTwistCoproduct([(line_bundle_2, line_bundle_3)], dimension_vector=[1], shift_vector=[0], degree=self.degree)
#         object3 = DerivedCategoryObject(string=f"Tw_{line_bundle_1.degree} Tw_{line_bundle_2.degree} O({line_bundle_3.degree})", catagory=self.catagory)

#         return DistinguishedTriangle(object1, object2, object3)
    
#     def chernCharacter(self):
#         r"""!
#         Method to compute the Chern character of the double spherical twist. The Chern character of the double
#         spherical twist is the Chern character of the third object in the distinguished triangle.

#         \return ChernCharacter The Chern Character of the double spherical twist
#         """


#         return self.defining_triangle.object3.chernCharacter()
    
#     def central_charge(self, *args):
#         r"""!
#         Method to compute the central charge of the spherical twist. The central charge of the spherical
#         twist is the central charge of the third object in the distinguished triangle.

#         \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
#                         For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
#                         and an integer representing the degree of the K3 surface.

#         \return complex The central charge of the spherical twist as a complex number

#         \throws TypeError If the args are not of the correct type
#         \throws ValueError If the number of args is incorrect
#         \throws NotImplementedError If the catagory of the object is not P1, P2, or K3

#         """

#         if self.catagory == 'P1':
#             if len(args) != 1:
#                 raise ValueError("Central charge of P1 requires single complex number parameter")
#             if not isinstance(args[0], complex):
#                 raise TypeError("P1 objects should have a single complex parameter")
            
#             ch = self.chernCharacter()
            
#             return -1*ch[1] + args[0]*ch[0]

#         # elif self.catagory == 'P2':
#         #     if len(args) != 2:
#         #         raise ValueError("Central charge of P2 requires two real number parameters")
#         #     if not all(isinstance(x, (float, int)) for x in args):
#         #         raise TypeError("P2 objects should have two real number parameters")
            
#         #     ch = self.chernCharacter()
            
#         #     return complex(-1*ch[2] +
#         #                     args[1] * ch[0],
#         #                       ch[1] - args[0] * ch[0])
            

#         elif self.catagory == 'K3':

#             if len(args) != 3:
#                 raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
#             if not all(isinstance(x, (float, int)) for x in args):
#                 raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
#             if not isinstance(args[2], int):
#                 raise TypeError("The degree of the K3 surface must be an integer")

#             alpha = args[0]
#             beta = args[1]
#             d = args[2]

#             ch = self.chernCharacter()
            
#             return complex(2*d*alpha * ch[1] - ch[2] - ch[0] + (beta**2 - alpha**2)*d*ch[0], 
#                            2*d*ch[1] - 2*d*alpha*beta*ch[0])
#         else:
#             raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")

        
    
#     def __str__(self):
#         r"""!
#         Returns a string representation of the spherical twist by printing the defining triangle

#         \return str A string representation of the spherical twist
#         """

#         return str(self.defining_triangle.object3)
    

#     def defining_triangle_to_json(self):
#         r"""!
#         Helper function to convert the data of the canonical triangle for a double spherical twist to a JSON string.
#         The data includes the dimensions, shifts, and line bundles from the first object of the triangle (i.e. the SphericalTwistSum
#         object), as well as a triple of the line bundle integers.

#         \return str A JSON string representation of the spherical twist triangle data
#         """

#         object1 = {
#             "shift_vector" : self.defining_triangle.object1.shift_vector,
#             "dimension_vector" : self.defining_triangle.object1.dimension_vector
#         }

#         chain_complex_data = {
#             "object1" : object1,
#             "degrees" : [self.line_bundle_1.degree,
#                         self.line_bundle_2.degree,
#                         self.line_bundle_3.degree]
#         }

#         return json.dumps(chain_complex_data)  



    

    


#     def secondary_triangle_to_json(self):
#         r"""!
#         Helper function to convert the data of the secondary canonical triangle to a JSON string. The data includes
#         the first object of the secondary canonical trianle encoded as a triple of lists, as well as an integer triple
#         consisting of the degrees of the three line bundles.

#         \return str A JSON string representation of the spherical twist triangle data
#         """

#         secondary_canonical_triangle = self.secondary_canonical_triangle()

#         object1 = {
#             "shift_vector" : secondary_canonical_triangle.object1.shift_vector,
#             "dimension_vector" : secondary_canonical_triangle.object1.dimension_vector
#         }

#         secondary_complex_data = {
#             "object1" : object1,
#             "degrees" : [self.line_bundle_1.degree,
#                         self.line_bundle_2.degree,
#                         self.line_bundle_3.degree]
#         }

#         return json.dumps(secondary_complex_data)    
            
            

        
#     def is_semistable(self, *args, logging=False, log_file=None):
#         r"""!
#         Method to check if the double spherical twist is semistable. The double spherical twist is semistable
#         if the Harder-Narasimhan factorization is trivial.

#         \param tuple args The parameters for the stability condition. The number of parameters depends on the catagory of the object
#                         For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
#                         and an integer representing the degree of the K3 surface.
#         \param bool logging A boolean flag to indicate whether to log the Harder-Narasimhan factors that caused the object to be unstable
#         \param str log_file The file to log the Harder-Narasimhan factors that caused the object to be unstable

#         \return bool True if the double spherical twist is semistable, False otherwise

#         \throws TypeError If the args are not of the correct type
#         \throws ValueError If the number of args is incorrect
#         \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
#         """

#         if self.catagory == 'P1':
#             if len(args) != 1:
#                 raise ValueError("Central charge of P1 requires single complex number parameter")
#             if not isinstance(args[0], complex):
#                 raise TypeError("P1 objects should have a single complex parameter")
#         # elif self.catagory == 'P2':
#         #     if len(args) != 2:
#         #         raise ValueError("Central charge of P2 requires two real number parameters")
#         #     if not all(isinstance(x, (float, int)) for x in args):
#         #         raise TypeError("P2 objects should have two real number parameters")
#         elif self.catagory == 'K3':
#             if len(args) != 3:
#                 raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
#             if not all(isinstance(x, (float, int)) for x in args):
#                 raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
#             if not isinstance(args[2], int):
#                 raise TypeError("The degree of the K3 surface must be an integer")
#         else:
#             raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        

#         try:
#             return len(self.get_HN_factors(*args)) == 1
#         except HarderNarasimhanError as e:
            
#             if logging and log_file:
#                 with open(log_file, 'a') as log_file:
#                     msg_str = e.message + f"@ {e.stability_parameters}"
#                     log_file.write(msg_str)
#             elif logging:
#                 msg_str = e.message + f"@ {e.stability_parameters}"
#                 print(msg_str)
            
#             return False
    



        

       
        

    
#     def get_HN_factors(self, *args):
#         r"""!
#         This is the crux of the DoubleSphericalTwist class, where we compute the Harder-Narasimhan factors of the
#         double spherical twist. A signficiant assumption (CONJECTURAL) that we make is that the Harder-Narasimhan filtration
#         must arise from the defining triangle or the secondary canonical triangle. This is not always the case, but we
#         have not yet implemented a general method to compute the HN factors in all cases.

#         This method works by first examining the two edge cases for the defining triangle: if the largest phase of the subobject
#         is less than the smallest phase of the quotient, then we assume the object is stable (CONJECTURAL). If the smallest phase
#         of the subobject is larger than the largest phase of the quotient, then we know for a fact (by BDL) that the Harder-Narasimhan
#         filtration can be computed by concatenating the HN factors of the subobject and quotient. If neither of these cases hold, we move
#         on to the second canonical triangle and apply a similar logic.

#         \param tuple args The parameters for the stability condition. The number of parameters depends on the catagory of the object
#                         For P1, this is a single complex number. For
#                         P2, this is two real numbers. For K3, this is two real numbers
#                         and an integer representing the degree of the K3 surface.

#         \return list A list of tuples where the first element is a DerivedCategoryObject and the second element is a float
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
#         # elif self.catagory == 'P2':
#         #     if len(args) != 2:
#         #         raise ValueError("Central charge of P2 requires two real number parameters")
#         #     if not all(isinstance(x, (float, int)) for x in args):
#         #         raise TypeError("P2 objects should have two real number parameters")
#         elif self.catagory == 'K3':
#             if len(args) != 3:
#                 raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
#             if not all(isinstance(x, (float, int)) for x in args):
#                 raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
#             if not isinstance(args[2], int):
#                 raise TypeError("The degree of the K3 surface must be an integer")
#         else:
#             raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        

        
#         ###########
#         # First check the defining triangle
#         ###########

#         # Write triangle as Tw_b O(c) -> Tw_a Tw b O(c) -> B + B[shift] + ...
#         subobject = SphericalTwist(self.line_bundle_2, self.line_bundle_3, self.degree)

#         modified_defining_triangle = self.defining_triangle.rotateLeft()
#         quotient_complex = modified_defining_triangle.object3 


#         if subobject.is_semistable(*args):
#             # get the phase of the single twist
#             left_side_phase = cmath.phase(subobject.central_charge(*args)) / math.pi
#             right_side_min_phase = quotient_complex.get_smallest_phase(*args)
#             right_side_max_phase = quotient_complex.get_largest_phase(*args)

#             if left_side_phase <= right_side_min_phase:
#                 potential_phase = cmath.phase(self.central_charge(*args))/math.pi
#                 for n in range(-4,4):
#                     if left_side_phase <= potential_phase + n and potential_phase + n <= right_side_min_phase:
#                         return [(self, potential_phase + n)]
                    
#                 raise HarderNarasimhanError(message=f"{subobject} is semistable and its phase is smaller than {quotient_complex}, but cannot correctly find the phase", 
#                                             stability_parameters=args)
                    

                
#             elif left_side_phase > right_side_max_phase:
#                 # By BDL20, the HN factors of the subobject and quotient concatenate to make
#                 # the HN factors of the twist
#                 if len(quotient_complex.dimension_vector) == 1:
#                     return [(subobject, left_side_phase),
#                              (quotient_complex, right_side_max_phase)]
#                 else:
#                     if len(quotient_complex) != 2:
#                         raise ValueError("The Hom object is not concentrated in 1 or 2 degrees")

#                     cplx_summand_0 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[0]],
#                                                 shift_vector=[quotient_complex.shift_vector[0]],
#                                                     dimension_vector=[quotient_complex.dimension_vector[0]])
#                     cplx_summand_1 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[1]],
#                                                 shift_vector=[quotient_complex.shift_vector[1]],
#                                                     dimension_vector=[quotient_complex.dimension_vector[1]])
#                     if right_side_min_phase == quotient_complex.sheaf_vector[0].phase(*args) + quotient_complex.shift_vector[0]:
#                         return [(subobject, left_side_phase),
#                                 (cplx_summand_1, right_side_max_phase),
#                                 (cplx_summand_0, right_side_min_phase)]
#                     else:
#                         return [(subobject, left_side_phase),
#                                 (cplx_summand_0, right_side_max_phase),
#                                 (cplx_summand_1, right_side_min_phase)]
                
#         else:
#             # Subobject (i.e. Tw_b O(c)) is not semistable, so we must
#             # first look at its HN filtration to see if anything can be salvaged
#             HN_factors_subobject = None
#             try:
#                 HN_factors_subobject = subobject.get_HN_factors(*args)
#             except HarderNarasimhanError as e:
#                 ## Add context to the raised error message
#                 raise HarderNarasimhanError(message=f"Subobject {subobject} in defining triangle of {self} is not semistable, cannot find its HN factors:\n{e.message}", 
#                                             stability_parameters=args)

#             right_side_min_phase = quotient_complex.get_smallest_phase(*args)
#             right_side_max_phase = quotient_complex.get_largest_phase(*args)

#             if HN_factors_subobject[0][1] <= right_side_min_phase:
#                 # largest phase of subobject is less than smallest phase of quotient
#                 potential_phase = cmath.phase(self.central_charge(*args))/math.pi
#                 for n in range(-3,3):
#                     if HN_factors_subobject[0][1] <= potential_phase + n and potential_phase + n <= right_side_min_phase:
#                         return [(self, potential_phase + n)]
                    
#             elif HN_factors_subobject[-1][1] > right_side_max_phase:
#                 # smallest phase of subobject is greater than largest phase of quotient
#                 if len(quotient_complex) == 1:
#                     return HN_factors_subobject + [(quotient_complex, right_side_max_phase)]
#                 elif len(quotient_complex) != 2:
#                     raise ValueError(f"The Hom object is not concentrated in 2 degrees; currently\n{quotient_complex}")

#                 cplx_summand_0 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[0]],
#                                             shift_vector=[quotient_complex.shift_vector[0]],
#                                                 dimension_vector=[quotient_complex.dimension_vector[0]])
#                 cplx_summand_1 = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[1]],
#                                             shift_vector=[quotient_complex.shift_vector[1]],
#                                                 dimension_vector=[quotient_complex.dimension_vector[1]])
#                 if right_side_min_phase == quotient_complex.sheaf_vector[0].phase(*args) + quotient_complex.shift_vector[0]:
#                     return HN_factors_subobject + [(cplx_summand_1, right_side_max_phase),
#                                                    (cplx_summand_0, right_side_min_phase)]
#                 else:
#                     return HN_factors_subobject + [(cplx_summand_0, right_side_max_phase),
#                                                    (cplx_summand_1, right_side_min_phase)]

        

#         ###########
#         # next check secondary canonical triangle
#         ###########

#         secondary_canonical_triangle = self.secondary_canonical_triangle().rotateLeft()
        

#         first_twist = SphericalTwist(secondary_canonical_triangle.object1.line_bundle_pairs_vector[0][0],
#                                      secondary_canonical_triangle.object1.line_bundle_pairs_vector[0][1],
#                                     self.degree)
        
#         HN_factors_first_term = None
#         try:
#             HN_factors_first_term = first_twist.get_HN_factors(*args)
#         except HarderNarasimhanError as e:
#             raise HarderNarasimhanError(message=f"Subobject {first_twist} in secondary canonical triangle of {self} is not semistable, cannot find its HN factors:\n{e.message}", 
#                                         stability_parameters=args)
#         HN_factors_last_term = None
#         try:
#             HN_factors_last_term = secondary_canonical_triangle.object3.get_HN_factors_ordered(*args)
#         except HarderNarasimhanError as e:
#             raise HarderNarasimhanError(message=f"Quotient {secondary_canonical_triangle.object3} in secondary canonical triangle of {self} is not semistable, cannot find its HN factors:\n{e.message}", 
#                                         stability_parameters=args)
        

#         smallest_HN_phase_first = HN_factors_first_term[-1][1]
#         largest_HN_phase_first = HN_factors_first_term[0][1]
#         smallest_HN_phase_last = HN_factors_last_term[-1][1]
#         largest_HN_phase_last = HN_factors_last_term[0][1]

#         if largest_HN_phase_first <= smallest_HN_phase_last:
#             potential_phase = cmath.phase(self.central_charge(*args))/math.pi
#             for n in range(-3,3):
#                 if largest_HN_phase_first <= potential_phase + n and potential_phase + n <= right_side_min_phase:
#                     return [(self, potential_phase + n)]
                
#         elif smallest_HN_phase_first > largest_HN_phase_last:
#             return HN_factors_first_term + HN_factors_last_term
        

#         raise HarderNarasimhanError(message=f"The Harder-Narasimhan factors of both quotients and both subobjects do not match any of 4 known scenarios; their HN factors necessarily intertwine.",
#                                      stability_parameters=args)
        









        



            

        


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