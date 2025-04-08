from GeometryContext import GeometryContext
from DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject
from GeometryContext import GeometryContext
from ChernCharacter import ChernCharacter
from CoherentSheaf import LineBundle
from HarderNarasimhanFiltration import HarderNarasimhanError, HarderNarasimhanFiltration
from SphericalTwist import SphericalTwistComposition

from typing import Dict

from dotenv import load_dotenv
import os

import math
import cmath

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']





class StabilityCondition():


    def __init__(self, geometry_context : GeometryContext, *parameters):


        if not isinstance(geometry_context, GeometryContext):
            raise TypeError("Geometry context must be of type GeometryContext")
        self.geometry_context = geometry_context

        if len(parameters) == 0:
            raise ValueError("No parameters given for stability condition")

        match geometry_context.catagory:
            case 'P1' | "LocalP1":
                if len(parameters) != 1:
                    raise ValueError("P1 stability condition requires a single complex number parameter")
                if not isinstance(parameters[0], complex):
                    raise TypeError("P1 stability condition requires a single complex number parameter")
            case 'P2' | "LocalP2":
                if len(parameters) != 2:
                    raise ValueError("P2 stability condition requires two real number parameters")
                if not all(isinstance(x, (float, int)) for x in parameters):
                    raise TypeError("P2 stability condition requires two real number parameters")
            case 'K3':
                if len(parameters) != 2:
                    raise ValueError("K3 stability condition requires three real number parameters: alpha, beta, and the degree")
                if not all(isinstance(x, (float, int)) for x in parameters):
                    raise TypeError("K3 stability condition requires three real number parameters: alpha, beta, and the degree")


            case _:
                raise NotImplementedError(f"Stability condition not implemented for {geometry_context.catagory}")
            
        self.parameters = parameters


    def centralCharge(self, derived_obj : DerivedCategoryObject) -> complex:
        r"""!
        
        """
        
        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Derived category object must be of type DerivedCategoryObject")

        chern_character = derived_obj.chernCharacter()
        polarization = self.geometry_context.polarization

        if self.geometry_context.catagory == 'P1' or self.geometry_context.catagory == 'LocalP1':
            ## TODO : Error checking?

            param_character = ChernCharacter(expr= -1 + self.parameters[0] * polarization, basis=[polarization], dimension=1)

            return self.geometry_context.divisor_data.evaluate((chern_character*param_character)[1])


        elif self.geometry_context.catagory == 'P2' or self.geometry_context.catagory == 'LocalP2':

            param_character = ChernCharacter(expr=-1 + 1j*polarization + (self.parameters[1] - self.parameters[1]*1j) * polarization**2, 
                                             basis=[polarization], dimension=2)

            return self.geometry_context.divisor_data.evaluate((chern_character * param_character)[2])
            
            ## TODO: figure out how to represent the central charge in the (s,q)-plane as a HRR evaluation
            
        elif self.geometry_context.catagory == 'K3':

            param_character = ChernCharacter.exp(  complex(self.parameters[0], self.parameters[1])*polarization )

            return self.geometry_context.divisor_data.evaluate( (param_character*chern_character)[2] )
        
        else :
            raise NotImplementedError(f"Central charge not implemented for {self.geometry_context.catagory}")


        #< a , b> * <ch0, ch1> = a ch1 + b ch0

    def phase(self, derived_obj : DerivedCategoryObject) -> float:
        r"""!
        Compute the phase of the derived object. The phase is computed as the argument of the central charge divided by pi.

        \param derived_obj The derived object to compute the phase for

        \return The phase of the derived object as a float

        \throws TypeError If the derived object is not of type DerivedCategoryObject
        """

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Derived category object must be of type DerivedCategoryObject")

        return cmath.phase(self.centralCharge(derived_obj)) / math.pi
    
    def is_semistable(self, derived_obj : DerivedCategoryObject) -> bool:
        r"""!
        Method to check if the spherical twist is stable. The spherical twist is stable if the Harder-Narasimhan
        filtration is trivial, i.e. just the object itself.

        \param derived_obj The derived object to check for stability

        \return True if the spherical twist is stable, False otherwise

        \throws TypeError If the derived object is not of type DerivedCategoryObject
        """

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Derived category object must be of type DerivedCategoryObject")

        return len(self.get_HN_factors(derived_obj)) == 1
        

    def get_HN_factors(self, derived_obj : DerivedCategoryObject) -> Dict[float, DerivedCategoryObject]:
        r"""!
        
        """

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Derived category object must be of type DerivedCategoryObject")
        
        if isinstance(derived_obj, LineBundle):
            ## Usually these will be stable for dim < 3 as long as the Picard rank is 1

            if derived_obj.geometry_context.divisor_data.variety_dimension <=2 and \
                len(derived_obj.geometry_context.divisor_data.basis) == 1:
                return HarderNarasimhanFiltration(stable_objects=[derived_obj],
                                                  phase_vector=[self.phase(derived_obj=derived_obj)])
            else:
                raise HarderNarasimhanError("Currently we can only confirm line bundles are semistable for Picard rank 1")
            
        elif isinstance(derived_obj, GradedCoproductObject):

            hn_filt = HarderNarasimhanFiltration([], [])
            for obj, shift, dim in derived_obj:

                obj_HN = self.get_HN_factors(obj)*dim
                obj_HN = obj_HN.shift(shift)
                hn_filt += obj_HN

        elif isinstance(derived_obj, SphericalTwistComposition):

            





    # def is_semistable(self, *args):
    #     r"""!
    #     Method to check if the spherical twist is stable. The spherical twist is stable if the Harder-Narasimhan
    #     filtration is trivial, i.e. just the object itself.

    #     \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
    #                     For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
    #                     and an integer representing the degree of the K3 surface.

    #     \return True if the spherical twist is stable, False otherwise

    #     \throws TypeError If the args are not of the correct type
    #     \throws ValueError If the number of args is incorrect
    #     \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
    #     """

    #     if self.catagory == 'P1':
    #         if len(args) != 1:
    #             raise ValueError("Central charge of P1 requires single complex number parameter")
    #         if not isinstance(args[0], complex):
    #             raise TypeError("P1 objects should have a single complex parameter")
    #     elif self.catagory == 'P2':
    #         if len(args) != 2:
    #             raise ValueError("Central charge of P2 requires two real number parameters. Currently {} parameters given: {}".format(len(args), args))
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("P2 objects should have two real number parameters")
    #     elif self.catagory == 'K3':

    #         if len(args) != 3:
    #             raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
    #         if not isinstance(args[2], int):
    #             raise TypeError("The degree of the K3 surface must be an integer")
    #     else:
    #         raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")

    #     try:
    #         return len(self.get_HN_factors(*args)) == 1
    #     except HarderNarasimhanError as e:
    #         print(f"Could not determine if {self} is semistable at {e.stability_parameters}: {e.message}")
    #         return False
    



    
    # def mass(self, *args):
    #     r"""!
    #     Computes the mass of an object in the derived catagory. The mass of a stable object is simply the modulus
    #     of its central charge. For a non-stable object, the mass is the sum of the masses of the Harder-Narasimhan
    #     factors of the object. The notion of the mass of an object is derived from string theory, where BPS states
    #     are characterized as objects which satisfy the BPS bound M = |Z|. For non-BPS states, one simply has 
    #     |Z| < M. The mass of a Bridgeland stability condition is given by the sum of its semistable factors, which 
    #     corresponds to the decomposition of an object in the derived category into Harder-Narasimhan factors.
    #     As a consequence, this method heavily relies on the get_HN_factors method to compute the Harder-Narasimhan
    #     factors of the object.

    #     \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
    #                     For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
    #                     and an integer representing the degree of the K3 surface.

    #     \return The mass of the object, as a non-negative real number 

    #     \throws TypeError If the args are not of the correct type
    #     \throws ValueError If the number of args is incorrect
    #     \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
    #     """


    #     if self.catagory == 'P1':
    #         if len(args) != 1:
    #             raise ValueError("Central charge of P1 requires single complex number parameter")
    #         if not isinstance(args[0], complex):
    #             raise TypeError("P1 objects should have a single complex parameter")
    #     elif self.catagory == 'P2':
    #         if len(args) != 2:
    #             raise ValueError("Central charge of P2 requires two real number parameters")
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("P2 objects should have two real number parameters")
    #     elif self.catagory == 'K3':

    #         if len(args) != 3:
    #             raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
    #         if not isinstance(args[2], int):
    #             raise TypeError("The degree of the K3 surface must be an integer")
    #     else:
    #         raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        

    #     try :
    #         HN_filtration = self.get_HN_factors(*args)

    #         mass = 0
    #         for (derived_cat_obj, _) in HN_filtration:
    #             mass += abs(derived_cat_obj.central_charge(*args))

    #         return mass
    #     except HarderNarasimhanError as e:
    #         print(f"Could not determine mass of {self} at {e.stability_parameters}: {e.message}")
    #         return -1
        
    

        
       



                
    # def get_HN_factors(self, *args):
    #     r"""!
    #     This method is the main workhorse of the SphericalTwist class. It computes the Harder-Narasimhan factors
    #     of the spherical twist object. It is generally assumed that for a single spherical twist, the only way
    #     that an object can destabilize is when an element of the last term of the defining triangle

    #                          O(a) -----> Tw_O(a) O(b) -------->  O(b)[n] ⊕ O(b)[n+1]

    #     has larger phase than O(a). In this case, the Harder-Narasimhan factors of the spherical twist depend on 
    #     which object it is that has larger phase. For example, if the minimum shift has larger phase, then we assume
    #     that the object must be strictly stable - THIS IS A CONJECTURE. If the maximum shift has smaller phase, then
    #     the triangle above leads to a Harder-Narasimhan filtration, so that the individual line bundle sums are precisely
    #     the HN factors; this is a result of Bapat-Deopurkar-Licata (2020).
         
    #     The most difficult case is when the smaller phase O(b)[n] is smaller than O(a) is smaller than O(b)[n+1]. In this 
    #     case, some homological algebra is required to show that the cone of the composed map 

    #                           Tw_O(a) O(b) -------> O(b)[n+1]

    #     fits into a distinguished triangle O(a) ----> Cone ----> O(b)[n]. 

    #     Instead of returning the objects alone, the method returns a list of tuples, where the first element is the semistable
    #     factor and the second element is the phase of the object. This is done to make computing largest and smallest semistable
    #     factors easier in the DoubleSphericalTwist class.

    #     The list is always returned in reverse order of the phase, so that the smallest phase HN factor is last and the largest
    #     is first.



    #     * It should be noted that several assumptions in this have not been verified outside of the quiver
    #     case


    #     \param args The parameters for the stability condition. The number of parameters depends on the catagory of the object
    #                     For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
    #                     and an integer representing the degree of the K3 surface.

    #     \return A list of tuples where the first element is a DerivedCategoryObject and the second element is a float
    #                 representing the phase of the object. The list is always returned in such a way that the largest phase
    #                 HN factor is first and smallest is last.

    #     \throws TypeError If the args are not of the correct type
    #     \throws ValueError If the number of args is incorrect
    #     \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
    #     \throws HarderNarasimhanError If the spherical twist is stable but the phase cannot be found

    #     """

    #     if self.catagory == 'P1':
    #         if len(args) != 1:
    #             raise ValueError("Central charge of P1 requires single complex number parameter")
    #         if not isinstance(args[0], complex):
    #             raise TypeError("P1 objects should have a single complex parameter")
    #     elif self.catagory == 'P2':
    #         if len(args) != 2:
    #             raise ValueError("Central charge of P2 requires two real number parameters")
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("P2 objects should have two real number parameters")
    #     elif self.catagory == 'K3':
    #         if len(args) != 3:
    #             raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
    #         if not isinstance(args[2], int):
    #             raise TypeError("The degree of the K3 surface must be an integer")
    #     else:
    #         raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        
        
    #     modified_defining_triangle = self.defining_triangle.rotateLeft()
    #     subobject = modified_defining_triangle.object1.sheaf_vector[0]

    #     quotient_complex = modified_defining_triangle.object3

    #     if subobject.phase(*args) <= quotient_complex.get_smallest_phase(*args):
    #         # The object is (ASSUMED TO BE --- CONJECTURE) stable
    #         potential_phase = cmath.phase(self.central_charge(*args)) / math.pi

    #         # Attempt to find the phase of the object; ideally this value of n should be unique

    #         # TODO: This is a temporary fix to the problem of finding the phase of the spherical twist object
    #         #       when the object is stable. We really shouldnt be considering odd dimensional shifts, but
    #         #       we run into an error when the phase of the twist is larger than both the subobject an quotient;
    #         #       this is a temporary fix to this occurse when the subobject and quotient differ by phase > 1 so 
    #         #       that they no longer lie in the same heart. In particular, this causes a discontinuity for the
    #         #       algebraic regions of the stability manifold.
    #         for n in range(-3,3):
    #             if subobject.phase(*args) <= potential_phase + n and potential_phase + n <= quotient_complex.get_largest_phase(*args):
    #                 return [(self, potential_phase + n)]
            
            
    #         raise HarderNarasimhanError(message=f"{self} should theoretically be stable, but could not find phase",
    #                                     stability_parameters=args)

            


    #     elif len(quotient_complex.dimension_vector) == 1:
    #         # The quotient object has only one term / is concentrated in a single degree and
    #         # its phase is smaller than the subobject.

    #         # The defining triangle O(a) -> Tw -> O(b)[shift] should in fact be the 
    #         # Harder-Narasimhan filtration in this case
            
    #         return [(modified_defining_triangle.object1, subobject.phase(*args)),
    #                  (quotient_complex, quotient_complex.get_smallest_phase(*args))]
            
    #     else:
    #         # Twist is unstable and the hom space is concentrated in more than one degree
    #         if len(quotient_complex.dimension_vector) != 2:
    #             raise ValueError("The Hom object is not concentrated in 1 or 2 degrees")
            
    #         phase0 = quotient_complex.sheaf_vector[0].phase(*args) + quotient_complex.shift_vector[0]
    #         phase1 = quotient_complex.sheaf_vector[1].phase(*args) + quotient_complex.shift_vector[1]

    #         largest_phase = max(phase0, phase1)

    #         # CASE 1: phi(subobj) > largest phase(quotient)
    #         if subobject.phase(*args) > largest_phase:
    #             # By BDL20, the HN factors of the subobject and quotient concatenate to make 
    #             # the HN factors of the twist
    #             if largest_phase == phase0:
    #                 return [(modified_defining_triangle.object1, subobject.phase(*args)), 
    #                         (CoherentSheafCoproduct(sheaf_vector=[quotient_complex.sheaf_vector[0]],
    #                                                 shift_vector=[quotient_complex.shift_vector[0]], 
    #                                                 dimension_vector=[quotient_complex.dimension_vector[0]]), phase0),
    #                         (CoherentSheafCoproduct(sheaf_vector=[quotient_complex.sheaf_vector[1]], 
    #                                                 shift_vector=[quotient_complex.shift_vector[1]], 
    #                                                 dimension_vector=[quotient_complex.dimension_vector[1]]), phase1)]
    #             else:
    #                 return [(modified_defining_triangle.object1, subobject.phase(*args)), 
    #                         (CoherentSheafCoproduct(sheaf_vector=[quotient_complex.sheaf_vector[1]], 
    #                                                 shift_vector=[quotient_complex.shift_vector[1]], 
    #                                                 dimension_vector=[quotient_complex.dimension_vector[1]]), phase1),
    #                         (CoherentSheafCoproduct(sheaf_vector=[quotient_complex.sheaf_vector[0]], 
    #                                                 shift_vector=[quotient_complex.shift_vector[0]], 
    #                                                 dimension_vector=[quotient_complex.dimension_vector[0]]), phase0)]

    #         # CASE 2: smallest phase(Quotient) < phi(subobj) < largest phase(quotient)
    #         #         this is the most difficult case to handle since we must in fact consider
    #         #         the cone of the composed map Tw_O(a) O(b) ----> O(b)[shift]
    #         else:
    #             if largest_phase == phase0:
    #                 smaller_idx = 1
    #                 larger_idx = 0
    #             else:
    #                 smaller_idx = 0
    #                 larger_idx = 1
                    
    #             # isolate single element of larger shift in the quotient object
    #             smaller_phase_complex = CoherentSheafCoproduct(sheaf_vector=[quotient_complex.sheaf_vector[smaller_idx]],
    #                                         shift_vector=[quotient_complex.shift_vector[smaller_idx]],
    #                                         dimension_vector=[quotient_complex.dimension_vector[smaller_idx]])
    #             larger_phase_complex = CoherentSheafCoproduct(sheaf_vector=[quotient_complex.sheaf_vector[larger_idx]],
    #                                         shift_vector=[quotient_complex.shift_vector[larger_idx]],
    #                                         dimension_vector=[quotient_complex.dimension_vector[larger_idx]])

    #             phase_subobject = subobject.phase(*args)
    #             phase_larger_complex = quotient_complex.sheaf_vector[larger_idx].phase(*args) + quotient_complex.shift_vector[larger_idx]

    #             central_charge_cone = larger_phase_complex.central_charge(*args) + subobject.central_charge(*args)
    #             cone_object = DerivedCategoryObject(string="Cone", catagory=self.catagory, chern_character=None)
    #             cone_triangle = DistinguishedTriangle(modified_defining_triangle.object1, cone_object, larger_phase_complex)
                
    #             # Need to compute phase of cone to make a StableObject
    #             phase_cone = cmath.phase(central_charge_cone) / math.pi
    #             # TODO: This is a temporary fix to the problem of finding the phase of the spherical twist object
    #             #       when the object is stable. We really shouldnt be considering odd dimensional shifts, but
    #             #       we run into an error when the phase of the twist is larger than both the subobject an quotient;
    #             #       this is a temporary fix to this occurse when the subobject and quotient differ by phase > 1 so 
    #             #       that they no longer lie in the same heart. In particular, this causes a discontinuity for the
    #             #       algebraic regions of the stability manifold.
    #             for n in range(-3,3):
    #                 if phase_subobject <= phase_cone + n and phase_cone + n <= phase_larger_complex:
    #                     return [(cone_triangle.object2, phase_cone + n),
    #                             (smaller_phase_complex, smaller_phase_complex.get_smallest_phase(*args))]
                    


    #             raise HarderNarasimhanError(message=f"Could not find phase of cone {cone_object} in \n{cone_triangle}",
    #                                         stability_parameters=args)



    #     def central_charge(self, *args) -> complex:
    #     r"""!
    #     Compute the central charge of an object in the derived category of coherent sheaves. For all the current categories
    #     implemented, the only stability conditions considered are numerical stability conditions; in particular, they only
    #     depend on the Chern character of the object. Since DerivedCategoryObjects are the highest level objects which 
    #     have a chernCharacter method, this method will be implemented here (the implementation does not change for any of 
    #     the children classes).

    #     \param tuple args The arguments required to compute the central charge. The number of arguments and the type
    #                  of arguments will depend on the catagory of the sheaves in the complex. For P1, the central
    #                  charge requires a single complex number. For P2, the central charge requires two floating-point
    #                  numbers. For K3, the central charge requires two floating-point numbers and an integer.

    #     \return complex The central charge of the chain complex as a complex number

    #     \throws ValueError If the number of arguments is incorrect
    #     \throws TypeError If the type of the arguments is incorrect
    #     \throws NotImplementedError If the catagory of the sheaves in the complex is not implemented
    #     """

        
    #     if self.catagory == 'P1':
    #         if len(args) != 1:
    #             raise ValueError("Central charge for P1 requires exactly one argument.")
    #         if not isinstance(args[0], complex):
    #             raise TypeError("Central charge for P1 requires a complex number as an argument.")

    #         ch = self.chernCharacter()
    #         return complex(-1*ch[1] + args[0] * ch[0])


    #     elif self.catagory == 'P2':
    #         if len(args) != 2:
    #             raise ValueError("Central charge for P2 requires exactly two arguments.")
    #         if not all(isinstance(arg, (float,int)) for arg in args):
    #             raise TypeError("Central charge for P2 requires two floating-point numbers as arguments.")
            
    #         ch = self.chernCharacter()
    #         return complex(-1*ch[2] + args[1] * ch[0],
    #                         ch[1] - args[0] * ch[0])
        
    #     elif self.catagory == 'K3':
    #         if len(args) != 3:
    #             raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
    #         if not isinstance(args[2], int):
    #             raise TypeError("The degree of the K3 surface must be an integer")

    #         alpha = args[0]
    #         beta = args[1]
    #         d = args[2]
    #         ch = self.chernCharacter()
            
    #         return complex(2*d*alpha * ch[1] - ch[2] - ch[0] + (beta**2 - alpha**2)*d*ch[0], 
    #                        2*d*ch[1] - 2*d*alpha*beta*ch[0])

    #     else:
    #         raise NotImplementedError("Central charge not implemented for this variety.")
    

    # def phase(self, *args) -> float:
    #     r"""!
    #     Computes the phase of the central charge of the coherent sheaf. The central charge
    #     is an element of the dual of the numerical Grothendieck group; in other words, a 
    #     funtction

    #     Z : K -> C

    #     where K is the numerical Grothendieck group, and C is the complex numbers. The phase
    #     of the central charge is the argument of this complex number.

    #     \param *args: float or int
    #         The parameters of the central charge. The number of parameters should be equal
    #         to the number of parameters required by the central charge for the given catagory.
    #         For example, a P1 object requires a single complex number parameter, while a P2
    #         object requires two real number parameters.

    #     \return float The phase of the central charge of the coherent sheaf, in units of pi
    #     """

    #     return cmath.phase(self.central_charge(*args)) / math.pi
    
        
    # @abstractmethod
    # def is_semistable(self, *args) -> bool:
    #     r"""!
    #     Method to determine if the derived category object is semistable with respect to a given stability condition. 
    #     This will simply act as a wrapper for the central charge method, which should be implemented in child classes.

    #     \param tuple args 
    #         The parameters of the stability condition. The number of parameters will depend on the catagory of the object.
    #         For P1, this will be a single complex number. For P2, this will be two real numbers. For K3, this will be
    #         two real numbers and one integer.

    #     \return bool True if the object is semistable with respect to the stability condition, False otherwise

    #     \throws ValueError
    #         If the DerivedCategoryObject is not initialized
    #         If the number of parameters is incorrect for the catagory
    #     \throws TypeError
    #         If the parameters are not of the correct type

    #     """

    #     pass

    # def mass(self, *args, logging : bool =False, log_file : bool =None) -> float:
    #     r"""!
    #     The mass of the double spherical twist is the sum of the masses of the Harder-Narasimhan factors; the 
    #     Harder-Narasimhan factors are assumed to come from either the defining triangle or secondary canonical triangle.
    #     The mass of the double spherical twist is computed by first computing the Harder-Narasimhan factors of the
    #     defining triangle, and then the secondary canonical triangle.

    #     \param tuple args The parameters for the stability condition. The number of parameters depends on the catagory of the object
    #                     For P1, this is a single complex number. For P2, this is two real numbers. For K3, this is two real numbers
    #                     and an integer representing the degree of the K3 surface.
    #     \param bool logging A boolean flag to indicate whether to log the Harder-Narasimhan factors that caused the object to be unstable
    #     \param str log_file The file to log the Harder-Narasimhan factors that caused the object to be unstable

    #     \return float The mass of the double spherical twist

    #     \throws TypeError If the args are not of the correct type
    #     \throws ValueError If the number of args is incorrect
    #     \throws NotImplementedError If the catagory of the object is not P1, P2, or K3
    #     """

    #     if self.catagory == 'P1':
    #         if len(args) != 1:
    #             raise ValueError("Central charge of P1 requires single complex number parameter")
    #         if not isinstance(args[0], complex):
    #             raise TypeError("P1 objects should have a single complex parameter")
    #     # elif self.catagory == 'P2':
    #     #     if len(args) != 2:
    #     #         raise ValueError("Central charge of P2 requires two real number parameters")
    #     #     if not all(isinstance(x, (float, int)) for x in args):
    #     #         raise TypeError("P2 objects should have two real number parameters")
    #     elif self.catagory == 'K3':
    #         if len(args) != 3:
    #             raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
    #         if not all(isinstance(x, (float, int)) for x in args):
    #             raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
    #         if not isinstance(args[2], int):
    #             raise TypeError("The degree of the K3 surface must be an integer")
    #     else:
    #         raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        
        

    #     try:
    #         HN_factors = self.get_HN_factors(*args)

    #         mass = 0
    #         for (derived_cat_obj, _) in HN_factors:
    #             mass += abs(derived_cat_obj.central_charge(*args))

    #         return mass
    #     except HarderNarasimhanError as e:
    #         if logging and log_file:
    #             with open(log_file, 'a') as log_file:
    #                 msg_str = e.message + f"@ {e.stability_parameters}"
    #                 log_file.write(msg_str)
    #         elif logging:
    #             msg_str = e.message + f"@ {e.stability_parameters}"
    #             print(msg_str)
    #         return -1
        