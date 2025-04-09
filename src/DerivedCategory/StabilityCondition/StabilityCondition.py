from src.DerivedCategory.GeometryContext import GeometryContext
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject, NumericalObject
from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.CoherentSheaf import LineBundle
from .HarderNarasimhanFiltration import HarderNarasimhanError, HarderNarasimhanFiltration
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition

from dotenv import load_dotenv
import os

import math
import cmath

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']





class StabilityCondition():

    _hn_filt_cache = {}


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

            param_character = ChernCharacter(expr= -1 + self.parameters[0] * polarization, basis=[polarization], dimension=1)

            return self.geometry_context.divisor_data.evaluate((chern_character*param_character)[1])


        elif self.geometry_context.catagory == 'P2' or self.geometry_context.catagory == 'LocalP2':

            ## We can represent the (s,q)-parameterization of Chunyi Li by (-ch2 + q ch0) +i (ch1 - s ch0)
            ## This is equivalent to simply evaluating the Chern character with the term
            ##            
            ##                    ∫ (-1 + iH + (q - is) H^2) * ch(E)
            ##
            ## Which is given below.
            param_character = ChernCharacter(expr=-1 + 1j*polarization + (self.parameters[1] - self.parameters[1]*1j) * polarization**2, 
                                             basis=[polarization], dimension=2)

            return self.geometry_context.divisor_data.evaluate((chern_character * param_character)[2])
            
            
        elif self.geometry_context.catagory == 'K3':

            ## The Mukai vector is simply <ch0, ch1, ch2 + ch0>. In addition, the mukai pairing is given by
            ## <r1, d1H, s1> * <r2, d2H, s2> = d1*d2 **H**2 - s1r2 - s2r1 . In order to represent this as a 
            ## normal evaluation of the chern character, we notice that we can change the middle term d1H to -d1H
            ## and then invert the whole expression. Thus,
            ##
            ##             ∫ exp((B + iω )H) * ch(E) td(X) = -1 *( exp( -1 * (b+iω) * H ) (ch0(E), ch1(E), ch2(E)) + ch0(E)   )
            ##                                             = <1, -bH - i ωH, -b^2/2 - iωbH - iω^2/2 + b^2/2> * v(E)

            param_character = ChernCharacter.exp( linear_expr=-1*complex(self.parameters[0], self.parameters[1])*polarization,
                                                 dimension=2 )

            purely_chern_comp = self.geometry_context.divisor_data.evaluate( (param_character*chern_character)[2] )

            return -1*(purely_chern_comp + chern_character[0])

            
        
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
        
        #############
        # TODO: Double check that this works as intended for graded coproducts
        #############

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
        

    def get_HN_factors(self, derived_obj : DerivedCategoryObject) -> HarderNarasimhanFiltration:
        r"""!
        
        """

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Derived category object must be of type DerivedCategoryObject")
        
        if derived_obj in self._hn_filt_cache:
            return self._hn_filt_cache[derived_obj]

        ## Otherwise, the object is not in the cache, so we need to compute it
        
        if isinstance(derived_obj, LineBundle):
            ## Usually these will be stable for dim < 3 as long as the Picard rank is 1

            if derived_obj.geometry_context.divisor_data.variety_dimension <=2 and \
                len(derived_obj.geometry_context.divisor_data.basis) == 1:
                self._hn_filt_cache[derived_obj] = HarderNarasimhanFiltration(stable_objects=[derived_obj],
                                                  phase_vector=[self.phase(derived_obj=derived_obj)])
            else:
                raise HarderNarasimhanError("Currently we can only confirm line bundles are semistable for Picard rank 1")
            
        elif isinstance(derived_obj, GradedCoproductObject):

            hn_filt = HarderNarasimhanFiltration([], [])
            for obj, shift, dim in derived_obj:

                obj_HN = self.get_HN_factors(obj)*dim
                obj_HN = obj_HN.shift(shift) ## Shift should also affect the phases of the objects as well
                hn_filt += obj_HN
            
            self._hn_filt_cache[derived_obj] = hn_filt

        elif isinstance(derived_obj, SphericalTwistComposition):

            stable = True
            for triangle in derived_obj.canonical_triangles:

                triangle = triangle.rotateLeft()
                subobject = triangle.object1
                quotient = triangle.object3

                if max(self.get_HN_factors(subobject)).phase < min(self.get_HN_factors(quotient)).phase:
                    continue
                elif min(self.get_HN_factors(subobject)).phase > max(self.get_HN_factors(quotient)).phase:
                    stable = False
                    self._hn_filt_cache[derived_obj] = self.get_HN_factors(subobject) + self.get_HN_factors(quotient)
                    break
                else:

                    if len(self.get_HN_factors(subobject)) > 1:
                        raise HarderNarasimhanError(f"The subobject {subobject} is not semistable, and the Harder-Narasimhan filtrations of {subobject} and {quotient} intertwine.")

                    ## At this point, we can assume that subobject is stable.
                    phi = self.phase(subobject)

                    upper_HN, lower_HN = self.get_HN_factors(quotient).splitAtPhase(phi)

                    new_quotient = upper_HN.toGradedCoproductObject()

                    ######################
                    # TODO: We are explicitly assuming that new_factor is semistable since 
                    #       this one particular triangle that it fits into does not destabilize it.
                    #       In reality, this is probably not the case outside of basic quivery 
                    #       catagories with only finitely many simple objects.
                    ######################

                    new_factor = NumericalObject(chern_char= (subobject.chernCharacter() + new_quotient.chernCharacter()),
                                                 geometry_context=derived_obj.geometry_context )
                    
                    new_factor_phase = self.phase(new_factor)
                    shift = 0
                    while new_factor_phase + shift < phi:
                        ###############
                        # TODO: This may need to be improved
                        ###############
                        shift += 1

                    self._hn_filt_cache[derived_obj] = lower_HN + HarderNarasimhanFiltration(stable_objects=[new_factor],
                                                  phase_vector=[new_factor_phase + shift])

            if stable:
                # The spherical twist is stable, so we can just return the HN factors of the object
                self._hn_filt_cache[derived_obj] = HarderNarasimhanFiltration(stable_objects=[derived_obj],
                                                  phase_vector=[self.phase(derived_obj=derived_obj)])
                    
            
        else:
            raise NotImplementedError(f"Harder-Narasimhan factors not implemented for {type(derived_obj)}")
        
        return self._hn_filt_cache[derived_obj]



            
    def mass(self, derived_obj : DerivedCategoryObject) -> float:

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Expected a DerivedCategoryObject.")

        return sum(
            abs(self.centralCharge(factor.obj)) * factor.multiplicity
            for factor in self.get_HN_factors(derived_obj)
        )


