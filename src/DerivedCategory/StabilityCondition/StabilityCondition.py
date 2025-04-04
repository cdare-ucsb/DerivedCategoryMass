

from dotenv import load_dotenv
import os

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']



class StabilityCondition():


    def __init__(self, *parameters):
        


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