from CoherentSheaf import LineBundle
from DerivedCategoryObject import DerivedCategoryObject
from DistinguishedTriangle import DistinguishedTriangle
from ChainComplex import ChainComplex
from ChernCharacter import ChernCharacter

import math
import cmath
import numpy as np
import json

import plotly.graph_objects as go
from plotly.graph_objs import *


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




class HarderNarasimhanError(Exception):
    r"""!
    Exception raised when the correct Harder-Narasimhan filtration cannot be found 
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args)

        self.message = kwargs.get('message') ## The error message

        self.stability_parameters = kwargs.get('stability_parameters') ## The stability parameters used to compute the Harder-Narasimhan factors




class SphericalTwist(DerivedCategoryObject):
    r"""!
    A spherical twist is a non-standard autoequivalence of the derived category of coherent sheaves, in the 
    sense that it does not arise from any composition of (1) standard autoequivalences on the variety (2) 
    tensoring by line bundles and (3) (co)homological shifts. Such autoequivalences typically only arise in the
    case of Calabi-Yau categories and toric varieties (where (-2)-curves can exist in the Fano setting as well, 
    e.g. P^2 blown up at 2 points), but often control the structure of the stability manifold. They are also 
    relevant to homological mirror symmetry since they are the derived equivalent of Dehn twists in the Fukaya
    category in the symplectic setting. They are explicitly defined as the cone of the derived evaluation morphism

            Hom(A, B) ⊗ A ---->  B ----> Tw_A B

    where A is a spherical object in the sense that RHom(A,A) is isomorphic as a graded-vector space to the 
    singular cohomology of an n-sphere.


    Currently, we only consider spherical twists around line bundles in the derived category of coherent sheaves
    since they always yield examples of spherical objects for Local P^n and K3 surfaces. On K3 surfaces, it is not
    true that the spherical twists account for all spherical objects, but they are still provide a rich source of
    examples to help predict mass growth.
    """
    
    def __init__(self, line_bundle_1, line_bundle_2, degree=1):
        r"""!
        Initialize an instance of SphericalTwist with the specified line bundles. The spherical twist
        is defined as the cone of the evaluation morphism 

                Hom(O(a), O(b)) ⊗ O(a) ---->  O(b) ----> Tw_a O(b)

        where O(a) and O(b) denote either line bundles or pushforwards of line bundles along the inclusion of the zero section. 
        The spherical twist is represented as a distinguished triangle in the derived category of coherent
        sheaves. 

        Several helper methods are used to compute the dimensions of the Hom spaces between the pushforwards
        of the line bundles, and then to construct the distinguished triangle.

        \param LineBundle line_bundle_1 The first line bundle in the Hom space

        \param LineBundle line_bundle_2: The second line bundle in the Hom space

        \param int degree An optional argument for the degree of the variety, which is relevant to the dimension of the
                       derived RHom space for K3 surfaces of picard rank 1. This does not affect the P1 or P2 implementations.
                       Default is 1.

        \throws TypeError If line_bundle_1 is not an instance of LineBundle
        \throws TypeError If line_bundle_2 is not an instance of LineBundle
        \throws ValueError If the line bundles are not defined on the same variety
        """

        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundle.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundle.")
        if line_bundle_1.catagory != line_bundle_2.catagory:
            raise ValueError("Line bundles must be defined on the same variety")

        self.line_bundle_1 = line_bundle_1 ## The first line bundle in the Hom space, i.e. O(a) in Hom(O(a), O(b))

        self.line_bundle_2 = line_bundle_2 ## The second line bundle in the Hom space, i.e. O(b) in Hom(O(a), O(b))

        self.catagory = line_bundle_1.catagory ## The catagory of the object, i.e. P1, P2, or K3

        self.degree = degree ## The degree of the K3 surface, if applicable

        self.defining_triangle = self._sph_twist_LineBundles(line_bundle_1, line_bundle_2) ## The distinguished triangle representing the spherical twist


    def __str__(self):
        r"""!
        Returns a string representation of the spherical twist by printing the defining triangle

        \return str A string representation of the spherical twist
        """

        return str(self.defining_triangle.object3)
    

    def defining_triangle_to_json(self):
        r"""!
        Method to convert the spherical twist to a JSON string. This is used to pass the chain complex
        data of the spherical twist in the browser. Since a spherical twist is defined by its defining
        triangle, 
        
                RHom(O(a), O(b)) ⊗ O(a) ----> O(b) ----> Tw_O(a) O(b)

        and Tw_O(a) O(b) is merely a symbol, the only data we really need to pass is what the first RHom
        object looks like (dimension, shifts, etc) and the degrees [a,b].

        \return str A JSON string representation of the chain complex. The first object contains the information
                    of the first RHom object, and the degrees of the line bundles are stored in the degrees key.
        """

        object1 = {
            "shift_vector" : self.defining_triangle.object1.shift_vector,
            "dimension_vector" : self.defining_triangle.object1.dimension_vector
        }

        chain_complex_data = {
            "object1" : object1,
            "degrees" : [self.line_bundle_1.degree,
                        self.line_bundle_2.degree]
        }

        return json.dumps(chain_complex_data)
    

    
    def _sph_twist_LineBundles(self, line_bundle_1, line_bundle_2):
        r"""!
        Helper function which uses the __dimHom_LineBundlesP2 method to compute the defining triangle for a 
        single spherical twist of a line bundle around another line bundle. The twist is given by

        Tw_lb1 lb2 = Cone(  Hom(lb_1, lb_2) ⊗ lb_1 ---->  lb_2 )


        \param LineBundle line_bundle_1: The first line bundle in the Hom space
        \param LineBundle line_bundle_2: The second line bundle in the Hom space

        \return DistinguishedTriangle The distinguished triangle representing the spherical twist
        """

        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundle.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundle.")

        homDims = _dimHom_LineBundles(line_bundle_1, line_bundle_2, self.degree)

        bundle_vector = []
        dimension_vector = [] 
        shift_vector = []

        # create the necessary lists for the ChainComplex constructor 
        for i in range(len(homDims)):
            if homDims[i] == 0:
                continue
            dimension_vector.append(homDims[i])
            shift_vector.append(-1*i)
            bundle_vector.append(LineBundle(line_bundle_1.degree, self.catagory))

        object1 = ChainComplex(sheaf_vector=bundle_vector, shift_vector=shift_vector, dimension_vector=dimension_vector)
        # While we technically can simply case object2 = line_bundle_2, we will use the ChainComplex constructor
        # since several successive methods call sheaf_vector[i] and shift_vector[i]. A significant overhaul would
        # be needed to change this.
        object2 = ChainComplex(sheaf_vector=[LineBundle(line_bundle_2.degree, self.catagory)], shift_vector=[0], dimension_vector=[1])
        object3 = DerivedCategoryObject(string=f"Tw_{line_bundle_1.degree} O({line_bundle_2.degree})", catagory=self.catagory)

        return DistinguishedTriangle(object1, object2, object3)
    
    def chernCharacter(self):
        r"""!
        Method to compute the Chern Character of the spherical twist. The Chern Character of the
        spherical twist is the Chern Character of the third object in the distinguished triangle.

        \return ChernCharacter The Chern Character of the spherical twist
        """

        return self.defining_triangle.object3.chernCharacter()
    
    def shift(self, n):
        r"""!
        Method to shift the spherical twist by n units. As a spherical twist is initially only
        specified as a string until its defining triangle is computed, the shift method simply
        relies on the implementation in the parent class DerivedCategoryObject.

        \param int n: The number of units to shift the object by

        \return DerivedCategoryObject The shifted object
        """

        if not isinstance(n, int):
            raise TypeError("Shift must be an integer")

        return self.defining_triangle.object3.shift(n)
    


    
    def central_charge(self, *args):
        r"""!
        Method to compute the central charge of the spherical twist. The central charge of the spherical
        twist is the central charge of the third object in the distinguished triangle.

        \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return The central charge of the spherical twist as a complex number 

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

        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
            
            ch = self.chernCharacter()
            
            return complex(-1*ch[2] +
                            args[1] * ch[0],
                              ch[1] - args[0] * ch[0])
            

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

        
    


    
    def is_semistable(self, *args):
        r"""!
        Method to check if the spherical twist is stable. The spherical twist is stable if the Harder-Narasimhan
        filtration is trivial, i.e. just the object itself.

        \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return True if the spherical twist is stable, False otherwise

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
                raise ValueError("Central charge of P2 requires two real number parameters. Currently {} parameters given: {}".format(len(args), args))
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

        try:
            return len(self.get_HN_factors(*args)) == 1
        except HarderNarasimhanError as e:
            print(f"Could not determine if {self} is semistable at {e.stability_parameters}: {e.message}")
            return False
    



    
    def mass(self, *args):
        r"""!
        Computes the mass of an object in the derived catagory. The mass of a stable object is simply the modulus
        of its central charge. For a non-stable object, the mass is the sum of the masses of the Harder-Narasimhan
        factors of the object. The notion of the mass of an object is derived from string theory, where BPS states
        are characterized as objects which satisfy the BPS bound M = |Z|. For non-BPS states, one simply has 
        |Z| < M. The mass of a Bridgeland stability condition is given by the sum of its semistable factors, which 
        corresponds to the decomposition of an object in the derived category into Harder-Narasimhan factors.
        As a consequence, this method heavily relies on the get_HN_factors method to compute the Harder-Narasimhan
        factors of the object.

        \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return The mass of the object, as a non-negative real number 

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
        

        try :
            HN_filtration = self.get_HN_factors(*args)

            mass = 0
            for (derived_cat_obj, _) in HN_filtration:
                mass += abs(derived_cat_obj.central_charge(*args))

            return mass
        except HarderNarasimhanError as e:
            print(f"Could not determine mass of {self} at {e.stability_parameters}: {e.message}")
            return -1
        
    

        
       



                
    def get_HN_factors(self, *args):
        r"""!
        This method is the main workhorse of the SphericalTwist class. It computes the Harder-Narasimhan factors
        of the spherical twist object. It is generally assumed that for a single spherical twist, the only way
        that an object can destabilize is when an element of the last term of the defining triangle

                             O(a) -----> Tw_O(a) O(b) -------->  O(b)[n] ⊕ O(b)[n+1]

        has larger phase than O(a). In this case, the Harder-Narasimhan factors of the spherical twist depend on 
        which object it is that has larger phase. For example, if the minimum shift has larger phase, then we assume
        that the object must be strictly stable - THIS IS A CONJECTURE. If the maximum shift has smaller phase, then
        the triangle above leads to a Harder-Narasimhan filtration, so that the individual line bundle sums are precisely
        the HN factors; this is a result of Bapat-Deopurkar-Licata (2020).
         
        The most difficult case is when the smaller phase O(b)[n] is smaller than O(a) is smaller than O(b)[n+1]. In this 
        case, some homological algebra is required to show that the cone of the composed map 

                              Tw_O(a) O(b) -------> O(b)[n+1]

        fits into a distinguished triangle O(a) ----> Cone ----> O(b)[n]. 

        Instead of returning the objects alone, the method returns a list of tuples, where the first element is the semistable
        factor and the second element is the phase of the object. This is done to make computing largest and smallest semistable
        factors easier in the DoubleSphericalTwist class.

        The list is always returned in reverse order of the phase, so that the smallest phase HN factor is last and the largest
        is first.



        * It should be noted that several assumptions in this have not been verified outside of the quiver
        case


        \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.

        \return A list of tuples where the first element is a DerivedCategoryObject and the second element is a float
                    representing the phase of the object. The list is always returned in such a way that the largest phase
                    HN factor is first and smallest is last.

        \throws TypeError If the args are not of the correct type
        \throws ValueError If the number of args is incorrect
        \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
        \throws HarderNarasimhanError If the spherical twist is stable but the phase cannot be found

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
        
        
        modified_defining_triangle = self.defining_triangle.rotateLeft()
        subobject = modified_defining_triangle.object1.sheaf_vector[0]

        quotient_complex = modified_defining_triangle.object3

        if subobject.phase(*args) <= quotient_complex.get_smallest_phase(*args):
            # The object is (ASSUMED TO BE --- CONJECTURE) stable
            potential_phase = cmath.phase(self.central_charge(*args)) / math.pi

            # Attempt to find the phase of the object; ideally this value of n should be unique

            # TODO: This is a temporary fix to the problem of finding the phase of the spherical twist object
            #       when the object is stable. We really shouldnt be considering odd dimensional shifts, but
            #       we run into an error when the phase of the twist is larger than both the subobject an quotient;
            #       this is a temporary fix to this occurse when the subobject and quotient differ by phase > 1 so 
            #       that they no longer lie in the same heart. In particular, this causes a discontinuity for the
            #       algebraic regions of the stability manifold.
            for n in range(-3,3):
                if subobject.phase(*args) <= potential_phase + n and potential_phase + n <= quotient_complex.get_largest_phase(*args):
                    return [(self, potential_phase + n)]
            
            
            raise HarderNarasimhanError(message=f"{self} should theoretically be stable, but could not find phase",
                                        stability_parameters=args)

            


        elif len(quotient_complex.dimension_vector) == 1:
            # The quotient object has only one term / is concentrated in a single degree and
            # its phase is smaller than the subobject.

            # The defining triangle O(a) -> Tw -> O(b)[shift] should in fact be the 
            # Harder-Narasimhan filtration in this case
            
            return [(modified_defining_triangle.object1, subobject.phase(*args)),
                     (quotient_complex, quotient_complex.get_smallest_phase(*args))]
            
        else:
            # Twist is unstable and the hom space is concentrated in more than one degree
            if len(quotient_complex.dimension_vector) != 2:
                raise ValueError("The Hom object is not concentrated in 1 or 2 degrees")
            
            phase0 = quotient_complex.sheaf_vector[0].phase(*args) + quotient_complex.shift_vector[0]
            phase1 = quotient_complex.sheaf_vector[1].phase(*args) + quotient_complex.shift_vector[1]

            largest_phase = max(phase0, phase1)

            # CASE 1: phi(subobj) > largest phase(quotient)
            if subobject.phase(*args) > largest_phase:
                # By BDL20, the HN factors of the subobject and quotient concatenate to make 
                # the HN factors of the twist
                if largest_phase == phase0:
                    return [(modified_defining_triangle.object1, subobject.phase(*args)), 
                            (ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[0]], shift_vector=[quotient_complex.shift_vector[0]], dimension_vector=[quotient_complex.dimension_vector[0]]), phase0),
                            (ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[1]], shift_vector=[quotient_complex.shift_vector[1]], dimension_vector=[quotient_complex.dimension_vector[1]]), phase1)]
                else:
                    return [(modified_defining_triangle.object1, subobject.phase(*args)), 
                            (ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[1]], shift_vector=[quotient_complex.shift_vector[1]], dimension_vector=[quotient_complex.dimension_vector[1]]), phase1),
                            (ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[0]], shift_vector=[quotient_complex.shift_vector[0]], dimension_vector=[quotient_complex.dimension_vector[0]]), phase0)]

            # CASE 2: smallest phase(Quotient) < phi(subobj) < largest phase(quotient)
            #         this is the most difficult case to handle since we must in fact consider
            #         the cone of the composed map Tw_O(a) O(b) ----> O(b)[shift]
            else:
                if largest_phase == phase0:
                    smaller_idx = 1
                    larger_idx = 0
                else:
                    smaller_idx = 0
                    larger_idx = 1
                    
                # isolate single element of larger shift in the quotient object
                smaller_phase_complex = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[smaller_idx]],
                                            shift_vector=[quotient_complex.shift_vector[smaller_idx]],
                                            dimension_vector=[quotient_complex.dimension_vector[smaller_idx]])
                larger_phase_complex = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[larger_idx]],
                                            shift_vector=[quotient_complex.shift_vector[larger_idx]],
                                            dimension_vector=[quotient_complex.dimension_vector[larger_idx]])

                phase_subobject = subobject.phase(*args)
                phase_larger_complex = quotient_complex.sheaf_vector[larger_idx].phase(*args) + quotient_complex.shift_vector[larger_idx]

                central_charge_cone = larger_phase_complex.central_charge(*args) + subobject.central_charge(*args)
                cone_object = DerivedCategoryObject(string="Cone", catagory=self.catagory, chern_character=None)
                cone_triangle = DistinguishedTriangle(modified_defining_triangle.object1, cone_object, larger_phase_complex)
                
                # Need to compute phase of cone to make a StableObject
                phase_cone = cmath.phase(central_charge_cone) / math.pi
                # TODO: This is a temporary fix to the problem of finding the phase of the spherical twist object
                #       when the object is stable. We really shouldnt be considering odd dimensional shifts, but
                #       we run into an error when the phase of the twist is larger than both the subobject an quotient;
                #       this is a temporary fix to this occurse when the subobject and quotient differ by phase > 1 so 
                #       that they no longer lie in the same heart. In particular, this causes a discontinuity for the
                #       algebraic regions of the stability manifold.
                for n in range(-3,3):
                    if phase_subobject <= phase_cone + n and phase_cone + n <= phase_larger_complex:
                        return [(cone_triangle.object2, phase_cone + n),
                                (smaller_phase_complex, smaller_phase_complex.get_smallest_phase(*args))]
                    


                raise HarderNarasimhanError(message=f"Could not find phase of cone {cone_object} in \n{cone_triangle}",
                                            stability_parameters=args)
                
                
                
                
            

        




class SphericalTwistSum(DerivedCategoryObject):
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



    def __str__(self):
        r"""!
        Returns a string representation of the spherical twist by printing the defining triangle. 
        The string representation is similar to that of the chain complex, where 2 lines are printed.
        The first line contains the number of times the spherical twist is applied, and the second line
        contains the actual twist. For example, the string representation of a spherical twist sum given
        by the data [(O(1), O(1)), [3], [-2]] would be

                ⊕3
        Tw_1 O(1)[-2]

        \return str A string representation of the spherical twist sum
        """

        # Initialize lists to hold each formatted component
        first_line_components = []
        second_line_components = []
        
        # Iterate over the zipped input lists
        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            # Format the second line component
            second_line_component = ""
            if s != 0:
                second_line_component = f"Tw_{lb1.degree} O({lb2.degree})[{s}]"
            else:
                second_line_component = f"Tw_{lb1.degree} O({lb2.degree})"

            second_line_components.append(second_line_component)
            
            # Calculate the position to place the '⊕n' above 'O'
            o_position = second_line_component.index(')')

            # Create a string with spaces up to the 'O' position, then add '⊕n'
            first_line_component = ""
            if n != 1:
                first_line_component = ' ' * o_position + f'⊕{n}'
            else:
                first_line_component = ' ' * o_position + '  '
            first_line_components.append(first_line_component)
            
        
        # Join all components with ' ⊕ ' separator
        first_line = ' '.join(first_line_components)
        second_line = ' ⊕ '.join(second_line_components)
        
        # Combine the two lines with a newline character
        if first_line.isspace():
            return second_line
        result = f"{first_line}\n{second_line}"
        
        return result
    

    def __len__(self):
        r"""!
        Returns the number of distinct spherical twists in the sum; we should generally expect that 
        the spherical twists are distinct since the constructor should hypothetically combine like terms.

        This method is primarily used in the DoubleSphericalTwist.get_HN_factors method to help determine
        edge cases.

        \return int The number of distinct spherical twists in the sum
        """

        return len(self.line_bundle_pairs_vector)
    
    def chernCharacter(self):
        r"""!
        Similar to the case of ChainComplex, since the Chern character is additive on exact sequences (i.e. factors
        through the Grothendieck group), we may always find the Chern character of an object obtained by sums and twists
        of known objects. In this case, we simply rely on the implementation of the above SphericalTwist class.

        \return ChernCharacter The Chern Character of the spherical twist sum
        """

        chern_character = ChernCharacter([0, 0, 0])
        
        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            sph_twist_chern = SphericalTwist(lb1, lb2, self.degree).chernCharacter()

            chern_character += int(n * (-1)**s) * sph_twist_chern
        
        return chern_character
    

    
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

        object2 = SphericalTwistSum([(line_bundle_2, line_bundle_3)], dimension_vector=[1], shift_vector=[0], degree=self.degree)
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



    
    def secondary_canonical_triangle(self):
        r"""!
        Method to compute the secondary canonical triangle of the spherical twist. The secondary canonical triangle
        is obtained by applying Tw_a to the defining triangle of Tw_b O(c); specifically, this gives

        Tw_a (Hom(O(b), O(c)) ⊗ O(b)) ----> Tw_a O(c) ----> Tw_a Tw_b O(c)

        \return DistinguishedTriangle The secondary canonical triangle of the spherical twist
        """

        # Get the object Tw_a O(c)
        object2 = SphericalTwistSum([(self.line_bundle_1, self.line_bundle_3)], dimension_vector=[1], shift_vector=[0], degree=self.degree)
    
        object3 = self.defining_triangle.object3

        single_twist_original = SphericalTwist(self.line_bundle_2, self.line_bundle_3, self.degree)
        hom_complex_single_twist = single_twist_original.defining_triangle.object1

        dimension_vec = []
        shift_vec = []
        line_bundle_pair_vec = []

        for i in range(len(hom_complex_single_twist.sheaf_vector)):
            if hom_complex_single_twist.dimension_vector[i] == 0:
                continue
            dimension_vec.append(hom_complex_single_twist.dimension_vector[i])
            shift_vec.append(hom_complex_single_twist.shift_vector[i])
            line_bundle_pair_vec.append((self.line_bundle_1, self.line_bundle_2))

        object1 = SphericalTwistSum(line_bundle_pair_vec, dimension_vec, shift_vec, self.degree)

        return DistinguishedTriangle(object1, object2, object3)
    


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
    
    def mass(self, *args, logging=False, log_file=None):
        r"""!
        The mass of the double spherical twist is the sum of the masses of the Harder-Narasimhan factors; the 
        Harder-Narasimhan factors are assumed to come from either the defining triangle or secondary canonical triangle.
        The mass of the double spherical twist is computed by first computing the Harder-Narasimhan factors of the
        defining triangle, and then the secondary canonical triangle.

        \param tuple args The parameters for the stability condition. The number of parameters depends on the catagory of the object
                        For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
                        and an integer representing the degree of the K3 surface.
        \param bool logging A boolean flag to indicate whether to log the Harder-Narasimhan factors that caused the object to be unstable
        \param str log_file The file to log the Harder-Narasimhan factors that caused the object to be unstable

        \return float The mass of the double spherical twist

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
            HN_factors = self.get_HN_factors(*args)

            mass = 0
            for (derived_cat_obj, _) in HN_factors:
                mass += abs(derived_cat_obj.central_charge(*args))

            return mass
        except HarderNarasimhanError as e:
            if logging and log_file:
                with open(log_file, 'a') as log_file:
                    msg_str = e.message + f"@ {e.stability_parameters}"
                    log_file.write(msg_str)
            elif logging:
                msg_str = e.message + f"@ {e.stability_parameters}"
                print(msg_str)
            return -1


        

       
        

    
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
#                  Static Helper Methods                          #
###################################################################


def _dimHom_LineBundles(line_bundle_1, line_bundle_2, degree_K3 = 1):
    r"""!
    General helper method to compute the dimension of the Hom spaces between 
    two line bundles (or pushforwards of line bundles, in the local Projectiv
    space case). The method is a wrapper for the specific implementations
    for P1, P2, and K3 surfaces.

    \param LineBundle line_bundle_1 The first line bundle in the Hom space
    \param LineBundle line_bundle_2 The second line bundle in the Hom space
    \param int degree_K3: The degree of the K3 surface, as an optional argument. This only affects
                      the K3 implementation. Default is 1.

    \return tuple A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles
    
    \throws TypeError If line_bundle_1 is not an instance of LineBundle
    \throws TypeError If line_bundle_2 is not an instance of LineBundle
    \throws ValueError If the line bundles are not defined on the same catagory
    \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
    """

    if not isinstance(line_bundle_1, LineBundle):
        raise TypeError("line_bundle_1 must be an instance of LineBundle.")
    if not isinstance(line_bundle_2, LineBundle):
        raise TypeError("line_bundle_2 must be an instance of LineBundle.")
    
    catagory = line_bundle_1.catagory

    if line_bundle_1.catagory != line_bundle_2.catagory:
        raise ValueError("Line bundles must be defined on the same variety")
    
    if catagory == 'P1':
        return _dimHom_LineBundlesP1(line_bundle_1, line_bundle_2)
    elif catagory == 'P2':
        return _dimHom_LineBundlesP2(line_bundle_1, line_bundle_2)
    elif catagory == 'K3':
        return _dimHom_LineBundlesK3(line_bundle_1, line_bundle_2, degree_K3)

    else:
        raise NotImplementedError("Only P1, P2 and K3 catagories are implemented")


def _dimHom_LineBundlesP1(line_bundle_1, line_bundle_2):
        r"""!
        Helper method which computes the dimension of the hom spaces between the pushforwards of the
        line bundles O(a) and O(b). The dimensions of the pushforwards are computed using the triangle

        i^* i_* E -> E -> E x O(2)[2]

        and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        the Hom spaces between line bundles on P1, the computation reduces to a case-by-case
        combinatorial problem. Since the homological index of the hom-space on P1 is bounded between
        0 and 1, the hom-space for local P1 is concentrated between degrees 0 and 2. Thus, we return
        a tuple of the form (a,b,c)

        \param LineBundle line_bundle_1 The first line bundle in the Hom space
        \param LineBundle line_bundle_2 The second line bundle in the Hom space

        \return tuple A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles

        \throws TypeError If line_bundle_1 is not an instance of LineBundle
        \throws TypeError If line_bundle_2 is not an instance of LineBundle
        """

        degree_dif = line_bundle_2.degree - line_bundle_1.degree

        if degree_dif == 0:
            return (1, 0, 1)
        if degree_dif >= 2:
            return (degree_dif + 1, degree_dif - 1, 0)
        elif degree_dif == 1:
            return (2, 0, 0)
        elif degree_dif == -1:
            return (0, 0, 2)
        else:
            return (0, -1*degree_dif - 1, -1*degree_dif + 1)



def _dimHom_LineBundlesP2(line_bundle_1, line_bundle_2):
        r"""!
        Helper method which computes the dimension of the hom spaces between the pushforwards of the
        line bundles O(a) and O(b). The dimensions of the pushforwards are computed using the triangle

        i^* i_* E -> E -> E x O(3)[2]

        and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        the Hom spaces between line bundles on P^2, the computation reduces to a case-by-case
        combinatorial problem. Since the homological index of the hom-space on P^2 is bounded between
        0 and 2, the hom-space for local P2 is concentrated between degrees 0 and 3. Thus, we return
        a tuple of the form (a,b,c,d)

        \param LineBundle line_bundle_1 The first line bundle in the Hom space
        \param LineBundle line_bundle_2 The second line bundle in the Hom space

        \return tuple A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles

        \throws TypeError If line_bundle_1 is not an instance of LineBundle
        \throws TypeError If line_bundle_2 is not an instance of LineBundle
        """

        degree_dif = line_bundle_2.degree - line_bundle_1.degree

        if degree_dif == 0:
            return (1, 0, 0, 1)
        elif degree_dif > -3 and degree_dif < 0:
            rank3 = math.comb(line_bundle_1.degree - line_bundle_2.degree + 2, 2)
            return (0, 0, 0, rank3)
        elif degree_dif > 0 and degree_dif < 3:
            rank0 = math.comb(degree_dif + 2, 2)
            return (rank0, 0, 0, 0)
        elif degree_dif >= 3:
            rank0 = math.comb(degree_dif + 2, 2)
            rank1 = math.comb(degree_dif - 1, 2)
            return (rank0, rank1, 0, 0)
        else:
            rank2 = math.comb(line_bundle_1.degree - line_bundle_2.degree -1, 2)
            rank3 = math.comb(line_bundle_1.degree - line_bundle_2.degree + 2, 2)    
            return (0, 0, rank2, rank3)
        
def _dimHom_LineBundlesK3(line_bundle_1, line_bundle_2, degree_K3):
    r"""!
    Helper method which computes the dimension of the hom spaces between line bundles on a 
    degree d K3 surface. Since the only K3s we consider are Picard rank 1, if O(a) and O(b) are
    distinct then either O(b-a) or O(a-b) must be ample; by Serre duality and a result of Huybrechts,
    one may always argue that Ext1(O(a), O(b)) = 0. In fact, Serre duality implies that the complex must
    be concentrated in a single degree, which is either 0 or 2 and corresponds to the cases that b > a and
    a < b, respectively. 

    We return a tuple of the form (n0,n1,n2) indicating the dimension of the graded RHom space. Aside from the
    instance that a=b, this will only include one nonzero element. Hirzebruch-Riemann-Roch shows that

    dim RHom(O(a), O(b)) = 1 + d(b-a)^2 

    \param LineBundle line_bundle_1 The first line bundle in the Hom space
    \param LineBundle line_bundle_2 The second line bundle in the Hom space
    \param int degree_K3 The degree of the K3 surface

    \return tuple A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles

    \throws TypeError If line_bundle_1 is not an instance of LineBundle
    \throws TypeError If line_bundle_2 is not an instance of LineBundle
    \throws TypeError If degree_K3 is not an integer
    """

    if not isinstance(line_bundle_1, LineBundle):
        raise TypeError("line_bundle_1 must be an instance of LineBundle.")
    if not isinstance(line_bundle_2, LineBundle):    
        raise TypeError("line_bundle_2 must be an instance of LineBundle.")
    if not isinstance(degree_K3, int):
        raise TypeError("The degree of the K3 surface must be an integer")

    degree_dif = line_bundle_2.degree - line_bundle_1.degree

    if degree_dif == 0:
        return (1, 0, 1)
    elif degree_dif > 0:
        return (degree_K3 * degree_dif**2 + 2, 0, 0)
    else:
        return (0, 0, degree_K3 * degree_dif**2 + 2)

        

def _dimHom_Line_and_SingleTwist(line_bundle_1, line_bundle_2, line_bundle_3, degree_K3):
    r"""!
    Helper method which computes the dimension of the hom spaces between a line bundle and 
    a single spherical twist of line bundles. This is a wrapper for the specific implementations
    for P1, P2, and K3 surfaces.

    \param LineBundle line_bundle_1 The line bundle to twist around: i.e. O(a) where we are computing Tw_a O(b)
    \param LineBundle line_bundle_2 The line bundle which the spherical twist is being applied to: i.e. O(b) where we are computing Tw_a O(b)
    \param LineBundle line_bundle_3: The line bundle which the spherical twist is being applied to: i.e. O(c) where we are computing Tw_a O(c)
    \param int degree_K3: The degree of the K3 surface, as an optional argument. This only affects
                      the K3 implementation. Default is 1.

    \return tuple A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles

    \throws TypeError If line_bundle_1 is not an instance of LineBundle
    \throws TypeError If line_bundle_2 is not an instance of LineBundle
    \throws TypeError If line_bundle_3 is not an instance of LineBundle
    \throws ValueError If the line bundles are not defined on the same catagory
    \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
    """

    if not isinstance(line_bundle_1, LineBundle):
        raise TypeError("line_bundle_1 must be an instance of LineBundle.")
    if not isinstance(line_bundle_2, LineBundle):
        raise TypeError("line_bundle_2 must be an instance of LineBundle.")
    if not isinstance(line_bundle_3, LineBundle):
        raise TypeError("line_bundle_3 must be an instance of LineBundle.")
    

    if line_bundle_1.catagory != line_bundle_2.catagory or line_bundle_1.catagory != line_bundle_3.catagory:
        raise ValueError("Line bundles must be defined on the same variety")
    

    if line_bundle_1.catagory == 'K3':
        return _dimHom_Line_and_SingleTwistK3(line_bundle_1, line_bundle_2, line_bundle_3, degree_K3)
    else:
        raise NotImplementedError("Double spherical twist is currently only implemented for K3 surfaces.")





def _dimHom_Line_and_SingleTwistK3(line_bundle_1, line_bundle_2, line_bundle_3, degree_K3):
    r"""!
    Helper function which computes the dimensions of the derived hom spaces between a line bundle 
    and a spherical twist of line bundles; that is, this function computes the dimensions of 
    
               RHom(O(a), Tw_{O(b)} O(c))
    
    which is necessarily concentrated in cohomological degrees -1 through 3. As a consequence, this 
    function returns a tuple of five integers (hom-1, hom0, hom1, hom2, hom3) corresponding to the 
    dimensions of the higher-ext groups. 

    \param LineBundle line_bundle_1 The line bundle to twist around: i.e. O(a) where we are computing Tw_a O(b)
    \param LineBundle line_bundle_2 The line bundle which the spherical twist is being applied to: i.e. O(b) where we are computing Tw_a O(b)
    \param LineBundle line_bundle_3 The line bundle which the spherical twist is being applied to: i.e. O(c) where we are computing Tw_a O(c)
    \param int degree_K3: The degree of the K3 surface, as an optional argument. This only affects
                      the K3 implementation. Default is 1.

    \return tuple A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles

    \throws TypeError If line_bundle_1 is not an instance of LineBundle
    \throws TypeError If line_bundle_2 is not an instance of LineBundle
    \throws TypeError If line_bundle_3 is not an instance of LineBundle
    \throws ValueError If the line bundles are not defined on the same catagory
    \throws NotImplementedError If the catagory of the object is not P1, P2, or K3

    """


    if not isinstance(degree_K3, int):
        raise TypeError("The degree of the K3 surface must be an integer")
    if not degree_K3 > 0:
        raise ValueError("The degree of the K3 surface must be a positive integer")
    if not isinstance(line_bundle_1, LineBundle):
        raise TypeError("line_bundle_1 must be an instance of LineBundleP2.")
    if not isinstance(line_bundle_2, LineBundle):
        raise TypeError("line_bundle_2 must be an instance of LineBundleP2.")
    if not isinstance(line_bundle_3, LineBundle):
        raise TypeError("line_bundle_3 must be an instance of LineBundleP2.")

    dif_12 = degree_K3 * (line_bundle_2.degree - line_bundle_1.degree)**2 + 2
    dif_23 = degree_K3 * (line_bundle_3.degree - line_bundle_2.degree)**2 + 2
    dif_13 = degree_K3 * (line_bundle_3.degree - line_bundle_1.degree)**2 + 2

    if line_bundle_2.degree < line_bundle_3.degree:

        if line_bundle_1.degree < line_bundle_2.degree:
            if dif_12 * dif_23 >= dif_13:
                return (dif_12 * dif_23 - dif_13, 0, 0, 0, 0)
            else:
                # There is an edge case when (b-c)(b-a) = -1; this reduces to Tw_{n+1} Tw_n O(n-1)
                # or Tw_{n-1} Tw_n O(n+1). In this case, the Hom space is concentrated in degree 0
                return (0, dif_13 - dif_12*dif_23, 0, 0, 0)
        elif line_bundle_1.degree == line_bundle_2.degree:
            return (0, 0, dif_23, 0, 0)
        elif line_bundle_1.degree > line_bundle_2.degree and line_bundle_1.degree < line_bundle_3.degree:
            return (0, dif_13, dif_12*dif_23, 0, 0)
        elif line_bundle_1.degree == line_bundle_3.degree:
            return (0, 1, dif_23 - 1, 0, 0)
        else:
            # line_bundle_1.deg > line_bundle_3.deg
            if dif_12 * dif_23 >= dif_13:
                return (0, 0, dif_23*dif_12-dif_13, 0, 0)
            else:
                # There is an edge case when (b-c)(b-a) = -1; this reduces to Tw_{n+1} Tw_n O(n-1)
                # or Tw_{n-1} Tw_n O(n+1). In this case, the Hom space is concentrated in degree 0
                # However, this theoretically should not occur here since we are considering the case
                # where a > c > b
                return (0, 0, 0, dif_13 - dif_12*dif_23, 0)
            
    elif line_bundle_2.degree == line_bundle_3.degree:

        # If b = c, then Tw_b O(c) = O(c)[-1]. Thus,
        #          Hom(O(a), Tw_b O(c)) = Hom(O(a), O(c)[-1]) = Hom(O(a), O(c))[-1] 

        (a,b,c) = _dimHom_LineBundlesK3(line_bundle_1, line_bundle_2, degree_K3)
        return (0, 0, a, b, c)
    
    else:
        # line_bundle_2.deg > line_bundle_3.deg
        if line_bundle_1.degree > line_bundle_2.degree:
            return (0, 0, 0, dif_13, dif_23*dif_13)
        elif line_bundle_1.degree == line_bundle_2.degree:
            return (0, 0, 0, 0, dif_23)
        elif line_bundle_1.degree < line_bundle_2.degree and line_bundle_1.degree > line_bundle_3.degree:
            if dif_12 * dif_23 >= dif_13:
                return (0, 0, dif_23*dif_12 - dif_13, 0, 0)
            else:
                return (0, 0, 0, dif_13 - dif_12*dif_23, 0)
        elif line_bundle_1.degree == line_bundle_3.degree:
            return (0, 1, dif_12*dif_23 - 1, 0, 0)
        else:
            # line_bundle_1.deg < line_bundle_3.deg
            return (0, dif_13, dif_12*dif_23, 0, 0)
        



            

        


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


    lb6 = LineBundle(1, catagory='K3')
    lb7 = LineBundle(2, catagory='K3')
    lb8 = LineBundle(4, catagory='K3')

    # sph_sum = SphericalTwistSum([(lb6, lb7), (lb8, lb9), (lb6, lb9), (lb6, lb7)], [1, 5, 0, 10], [-1, 2, -5, -1], degree=1)
    # print(sph_sum)
    # print(sph_sum.chernCharacter())

    


    # sph3 = SphericalTwist(lb7, lb8, degree=2)




    # print(sph3.shift(3))
    # print("Is semistable: ", sph3.is_semistable(1, 2, 1))

    K3_deg = 1
    

    sph4 = DoubleSphericalTwist(lb6, lb7, lb8, degree=K3_deg)

    # print(sph4.mass(2, 3, K3_deg))

    x_vals = np.linspace(-5, 11.10, 200)  # X values from -2 to 2

    # Generate y values satisfying y > x^2
    y_vals = []
    for x in x_vals:
        y_range = np.linspace(0, 5, 100)  # 50 points per x value
        y_vals.append(y_range)

    # Convert to numpy array
    y_vals = np.array(y_vals).flatten()  # Flatten the y array

    # Repeat x values to match the shape of y
    x_vals = np.repeat(x_vals, 100)  # Each x value repeats 10 times

    masses = np.array([sph4.mass(x, y, K3_deg) for x, y in zip(x_vals, y_vals)])

    # Plot the surface
    fig = go.Figure(data=[go.Scatter3d(z=masses, x=x_vals, y=y_vals,
                                    mode='markers', marker=dict(size=3, color=masses, colorscale='viridis'))])

    fig.update_layout(scene = dict(xaxis = dict(nticks=4, range=[-5,5])))
   
    fig.show()
