from src.DerivedCategory.GeometryContext import GeometryContext
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject, NumericalObject
from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.CoherentSheaf import LineBundle, CoherentSheaf
from .HarderNarasimhanFiltration import HarderNarasimhanError, HarderNarasimhanFiltration
from .SlopeStability import SlopeStability
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition


from dotenv import load_dotenv
import os
from typing import Iterator

from sympy import Add
from itertools import product

from functools import reduce

import math
import cmath
import numpy as np

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['P1', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']


class UnstableInstanceException(Exception):
    r"""!
    Exception raised when a function meant for stable objects is called on an unstable object.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args)

        self.message = kwargs.get('message') ## The error message



class StabilityCondition():

    _hn_filt_cache = {}
    _slope_hn_filt_cache = {}


    def __init__(self, geometry_context : GeometryContext, *parameters, P2_sq_param : bool = False):

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
                self.P2_sq_param = P2_sq_param
            case 'K3':
                basis = self.geometry_context.divisor_data.basis
                expected_len = len(basis) + 1
                if len(parameters) != expected_len:
                    raise ValueError(f"K3 stability condition requires {expected_len} real parameters: one for each B-field component and one ω scalar")
                if not all(isinstance(x, (float, int)) for x in parameters):
                    raise TypeError("K3 stability parameters must be real numbers")
                
                if not parameters[-1] > 0:
                    raise ValueError("The last parameter (ω) must be positive")


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

            if self.P2_sq_param:
                param_character = ChernCharacter(expr=-1 + 1j*polarization + (self.parameters[1] - self.parameters[1]*1j) * polarization**2, 
                                             basis=[polarization], dimension=2)
            else:
                param_character = ChernCharacter.exp(linear_expr= (self.parameters[0] + 1j* self.parameters[1]) * polarization,
                                                     basis=[polarization], dimension=2)

            return self.geometry_context.divisor_data.evaluate((chern_character * param_character)[2])
            
            
        elif self.geometry_context.catagory == 'K3':

            ## The Mukai vector is simply <ch0, ch1, ch2 + ch0>. In addition, the mukai pairing is given by
            ## <r1, d1H, s1> * <r2, d2H, s2> = d1*d2 **H**2 - s1r2 - s2r1 . In order to represent this as a 
            ## normal evaluation of the chern character, we notice that we can change the middle term d1H to -d1H
            ## and then invert the whole expression. Thus,
            ##
            ##             ∫ exp((B + iω )H) * ch(E) td(X) = -1 *( exp( -1 * (B + i ωH) ).(ch0(E), ch1(E), ch2(E)) + ch0(E) )
            ##                                             = <1, -bH - i ωH, -b^2/2 - iωbH - iω^2/2 + b^2/2> * v(E)

            from sympy import I

            basis = self.geometry_context.divisor_data.basis
            b_coeffs = self.parameters[:-1]
            omega = self.parameters[-1]

            B = sum(bi * Di for bi, Di in zip(b_coeffs, basis))  # B-field
            W = omega * self.geometry_context.polarization       # ω H
            twist = -1 * (B + I * W)


            param_character = ChernCharacter.exp( linear_expr=twist, basis=basis,
                                                 dimension=2 )
            

            purely_chern_comp = self.geometry_context.divisor_data.evaluate( (param_character*chern_character)[2] )


            z_sym = -1*(purely_chern_comp + chern_character[0])
            return complex(z_sym.evalf())

            
        
        else :
            raise NotImplementedError(f"Central charge not implemented for {self.geometry_context.catagory}")


        #< a , b> * <ch0, ch1> = a ch1 + b ch0



    def tiltedSlope(self, coh : CoherentSheaf) -> float:

        if not isinstance(coh, CoherentSheaf):
            raise TypeError("Coherent sheaf must be of type CoherentSheaf")
        if not coh.geometry_context == self.geometry_context:
            raise TypeError("Coherent sheaf must be of the same geometry context as the stability condition")
        if self.geometry_context.catagory == 'P1' or self.geometry_context.catagory == 'LocalP1':
            return coh.slope
        
        B = Add(*[coeff * basis_elem for coeff, basis_elem in zip(self.parameters[:-1], self.geometry_context.divisor_data.basis)] )
        W = self.parameters[-1] * self.geometry_context.polarization

        numerator =  self.geometry_context.divisor_data.evaluate( (coh.c1 - coh.rank * B), W)

        return numerator / coh.rank


    def phase(self, derived_obj : DerivedCategoryObject) -> float:
        r"""!
        Compute the phase of the derived object. The phase is computed as the argument of the central charge divided by pi.

        \param derived_obj The derived object to compute the phase for

        \return The phase of the derived object as a float

        \throws TypeError If the derived object is not of type DerivedCategoryObject
        """

        if isinstance(derived_obj, (LineBundle, CoherentSheaf, NumericalObject)):
            return cmath.phase(self.centralCharge(derived_obj)) / math.pi
        elif isinstance(derived_obj, GradedCoproductObject):

            if len(derived_obj) > 1:
                raise UnstableInstanceException("Graded coproduct whose objects are concentrated in multiple slopes are unstable; there is no coherent notion of phase for them.")
            else:
                return self.phase(derived_obj.object_vector[0]) + derived_obj.shift_vector[0]
            
        elif isinstance(derived_obj, SphericalTwistComposition):

            if not self.is_semistable(derived_obj):
                raise UnstableInstanceException("Spherical twist is unstable; cannot compute phase.")
            
            return self.get_HN_factors(derived_obj)[0].phase
                

        
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
        
        if not derived_obj.geometry_context == self.geometry_context:
            raise TypeError("Derived category object must be of the same geometry context as the stability condition")
        
        if derived_obj in self._hn_filt_cache:
            return self._hn_filt_cache[derived_obj]

        ## Otherwise, the object is not in the cache, so we need to compute it
        
        if isinstance(derived_obj, LineBundle):
            ## Usually these will be stable for dim < 3 as long as the Picard rank is 1

            if derived_obj.geometry_context.divisor_data.variety_dimension <=2 and \
                len(derived_obj.geometry_context.divisor_data.basis) == 1:
                ## For Picard rank 1 surfaces and curves, line bundles are always stable

                self._hn_filt_cache[derived_obj] = HarderNarasimhanFiltration(stable_objects=[derived_obj],
                                                  phase_vector=[self.phase(derived_obj=derived_obj)])
                

            elif self.geometry_context.catagory == 'K3':
                ## In higher picard rank K3 surfaces, line bundles CAN destabilize

                candidates = self._K3_line_bundle_destabilizing_candidates(derived_obj,
                                                                        r_max = 5,
                                                                        require_effective = True)
                maximal_destabilizing_subsheaf = max(candidates, key=lambda E: self.phase(E), default=None)

                
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


            if len(derived_obj) == 1:
                tri = derived_obj.defining_triangle

                if tri.object2 in tri.object1:
                    cone = tri.object1 - tri.object2 
                    cone=cone.shift(1)
                    return self.get_HN_factors(cone)
                elif tri.object1.get(0) is not None and tri.object1.get(1) is None:
                    ## The RHom complex is concentrated in degree 0.
                    ## We are working with a Kernel sheaf.

                    kernel_ch = tri.object1.get(0).chernCharacter() - tri.object2.chernCharacter()
                    kernel_sheaf = CoherentSheaf(kernel_ch, geometry_context=derived_obj.geometry_context)

                    return self.get_HN_factors(kernel_sheaf).shift(1)


            else:       
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
                        ## The Harder-Narasimhan filtrations of the subobject and quotient intertwine

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
                    ## TODO: fix phase here
                    self._hn_filt_cache[derived_obj] = HarderNarasimhanFiltration(stable_objects=[derived_obj],
                                                    phase_vector=[1])
                    
            
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



    def _K3_line_bundle_destabilizing_candidates(self, L : LineBundle, r_max : int, require_effective : bool = True) -> Iterator[LineBundle]:
        
        if not isinstance(L, LineBundle):
            raise TypeError("Expected a LineBundle.")
        
        basis = self.geometry_context.divisor_data.basis

        ###
        # Step 1: Convert lb to a list of ints representing the coefficients of the basis
        ###
        L_coeffs = [int(L.chernCharacter()[b]) for b in basis]

        ### Get all effective ranges
        ranges = [range(di + 1) for di in L_coeffs]  # assume effective basis and E <= D




        slope_stab = SlopeStability(self.geometry_context)

        print("About to call SlopeStability.get_HN_factors")

        L_slope_HN = slope_stab.get_HN_factors(L)

        print("Done calling SlopeStability.get_HN_factors")


        destabilizing_objects = []

        

        if all(self.tiltedSlope(factor.obj) > 0 for factor in L_slope_HN):
            #####################################
            # CASE 1: The twisted slope of L is positive, so L lies in the torision part T_B,ω =
            #           sheaves whose Harder–Narasimhan factors all have slope μB,ω>0
            #####################################

            ## We can now use some homological-algebra to simplify the problem. Foremost, we 
            ## only care about stable objects in the heart A_{B, ω} that can destabilize L. Such 
            ## objects A_{B, ω} are are either shifts by [1] of slope-stable objects with μB,ω <=0
            ## or slope-stable objects with μB,ω > 0.
            ##
            ## In the first case, there should be no morphisms from E[1] -> L, since that would imply
            ## RHom(E[1], L) = Ext^{-1}(E, L) is non-zero. However, the Ext of coherent sheaves is always
            ## Concentrated in degree 1. 
            ##
            ## In the latter case, we are simply looking at a morphism E -> L where E is a slope stable 
            ## sheaf with slope <= 0. In order for there to exist a non-zero morphism E -> L, we must
            ## have H^0(X, L - E) is non-zero. When E is a line bundle, this simply requires some numerics 
            ## such as chi( L - E) > 0 and E <= L. When E is a vector bundle, we must ensure that E->L actually
            ## gives a monomorphism. This happens when the of f:E->L is in the heart A_{B, ω}. Now the cone(f)
            ## is a 2 term complex [A -> B] with H^{-1}(cone(f)) = ker(f) and H^0(cone(f)) = coker(f). This
            ## is in the heart A_{B, ω} if and only if ker(f) is in the torsion free part F_B,ω and coker(f) is in the
            ## torsion part T_B,ω. This is equivalent to the condition that the twisted slope of the kernel is 
            ## <= 0 and the twisted slope of the cokernel is > 0. Thankfully, since our target L is a line bundle,
            ## the cokernel is always a torsion sheaf so that it is in T_B,ω. Thus, we only need to check
            ## that the kernel is in the heart A_{B, ω} and has twisted slope <= 0. 

            print("Case 1: The twisted slope of L is positive, so L lies in the torision part T_B,ω")


            #####
            #
            # STEP (a) : Look for line bundles E in T_B,ω that may destabilize L
            #      ----------------------------------------
            #
            #     The product function actually makes this an n-fold loop; specifically, its growth rate
            #     is roughly O(n^r) where r is the Picard rank
            #####
            for coeffs in product(*ranges):
                E_coeffs = np.array(coeffs)
                if require_effective and np.any(E_coeffs < 0):
                    continue
                if np.all(E_coeffs == 0):
                    continue
                E_symb = Add(*[ei * Di for ei, Di in zip(E_coeffs, basis)])

                # if self.geometry_context.divisor_data.evaluate(E_symb**2) != -2:
                #     continue

                E = LineBundle(E_symb, geometry_context=L.geometry_context)

                if self.tiltedSlope(E) <= 0:
                    ## E is in the torsion-free part F_B,ω
                    continue

                if self.phase(E) <= self.phase(L):
                    ## Phase is smaller anyways, dont need it
                    continue

                if len(L_slope_HN) == 1 and E.slope >= L.slope:
                    ## If L is slope-stable and E has a larger slope, it
                    ## is a classical fact of slope stability that there 
                    ## cannot exist a morphism from E->L
                    continue

                diff_class = L.divisor - E_symb

                ## Compute to see if there indeed exists a morphism
                RR_calc = 2 + 0.5*float(L.geometry_context.divisor_data.evaluate(diff_class**2).evalf())
                if RR_calc > 0:
                    destabilizing_objects.append(E)
            

            B = Add(*[coeff * basis_elem for coeff, basis_elem in zip(self.parameters[:-1], basis)] )
            W = self.parameters[-1] * self.geometry_context.polarization

            print("Moving on to higher rank bundles")

            #####
            #
            # STEP (b) : Look for higher-rank vector bundles E in T_ B,ω that may destabilize L
            #      ----------------------------------------
            #
            #####
            for r in range(2, r_max+1):

                ## The product function actually makes this an n-fold loop; specifically, its growth rate
                ## is roughly O(n^r) where r is the Picard rank
                for coeffs in product(*ranges):

                    E_c1_coeffs = np.array(coeffs)
                    if require_effective and np.any(E_c1_coeffs < 0):
                        continue
                    if np.all(E_c1_coeffs == 0):
                        continue
                    E_c1 = Add(*[ei * Di for ei, Di in zip(E_c1_coeffs, basis)])

                    E_slope = float(self.geometry_context.divisor_data.evaluate( (E_c1 - r*B), W).evalf())/r
                    if E_slope <= 0:
                        ## E is in the torsion-free part F_B,ω, which means its really E[1]
                        continue


                    ## The cone is a 2 term complex [A -> B] with H^{-1}(cone(f)) = ker(f) and H^0(cone(f)) = coker(f). This
                    ## is in the heart A_{B,ω} if and only if ker(f) is in the torsion free part F_B,ω and coker(f) is in the
                    ## torsion part T_B,ω.
                    ker_c1 = E_c1 - L.divisor

                    ker_tilt_slope =  float(self.geometry_context.divisor_data.evaluate( (ker_c1 - (r-1) * B), W).evalf())/ (r-1)

                    if ker_tilt_slope > 0:
                        ## The kernel is in the torsion part T_B,ω but we only want to look at the 
                        ## torsion-free part <= 0
                        continue

                    ## We now know that E -> L is a monomorphism, so lastly want to figure out what the second chern character is
                    ## We can use a Bogomolov-Gieseker inequality to bound the second chern character. In fact, a 
                    ## slightly stricter bound guaranteeing non-emptyness of moduli space is <v, v> >= -2, which
                    ## evaluates as follows:

                    bound = (float(self.geometry_context.divisor_data.evaluate( E_c1**2 ).evalf() ) - 2*r**2 + 2) / 2*r

                    for ch2 in range(-1*int(bound), int(bound)+1):
                    
                        mukai_vector_coeffs = [r] + list(coeffs) + [ch2]
                        if reduce(math.gcd, mukai_vector_coeffs) != 1:
                            ## In order for the moduli space to be nonempty, the mukai vector must be
                            ## primitive
                            continue


                        ## Use Riemann-Roch for χ(X, E^v(L)) = h^0 - h^1 + h^2; 
                        chi = -1*float(self.geometry_context.divisor_data.evaluate( E_c1 * L.divisor ).evalf()) \
                            + r*float(self.geometry_context.divisor_data.evaluate(L.divisor, L.divisor).evalf())/2 + 2*r + ch2
                        
                        if chi <=0:
                            ## Riemann Roch tells us there are no morphisms
                            continue


                        H2 = self.geometry_context.polarization**2
                        ch2 = ch2 / float(self.geometry_context.divisor_data.evaluate(H2).evalf())

                        chern = ChernCharacter(expr=r + E_c1 + ch2 * H2,
                                            basis=basis,
                                            dimension=2)
                        destabilizing_sheaf = CoherentSheaf(chern, geometry_context=L.geometry_context)

                        if self.phase(destabilizing_sheaf) > self.phase(L):
                            destabilizing_objects.append(destabilizing_sheaf)






            
        elif all(self.tiltedSlope(factor.obj) <= 0 for factor in L_slope_HN):
            #####################################
            # CASE 2: The twisted slope of L is negative, so L lies in the torsion-free part F_B,ω =
            #           sheaves whose Harder–Narasimhan factors all have slope μB,ω<=0
            #####################################

            ## We can now use some homological-algebra to simplify the problem. Foremost, we 
            ## only care about stable objects in the heart A_{B, ω} that can destabilize L. Such 
            ## objects A_{B, ω} are are either shifts by [1] of slope-stable objects with μB,ω <=0
            ## or slope-stable objects with μB,ω > 0.
            ##
            ## This is much more difficult than the previous case since there can be morphisms E -> L[1], 
            ## and E[1] -> L[1]. The morphisms E[1] -> L[1] can be found the same way as in the previous case;
            ## any such morphism must come from a morphism E -> L, and in order for E[1] to be in the heart in 
            ## the first place it must also be a coherent sheaf with twisted slope μB,ω <= 0. However,
            ## the morphisms E -> L[1] are much more difficult to find; such a morphism E -> L[1] arises
            ## from divisors with nontrivial first cohomology (e.g. H^1(X, L - E) != 0). As pointed out in
            ## RHom.py, there is no way for us to precisely compute the dimensionality of H^1(X, L - E) using 
            ## numerical / divisor data alone.
            

            print("Case 2: The twisted slope of L is negative, so L lies in the torsion-free part F_B,ω")
            

            #####
            #
            # STEP (a) : Look for line bundles that may destabilize L
            #      ----------------------------------------
            #
            # The product function actually makes this an n-fold loop; specifically, its growth rate
            # is roughly O(n^r) where r is the Picard rank
            #####
            for coeffs in product(*ranges):
                E_coeffs = np.array(coeffs)
                if require_effective and np.any(E_coeffs < 0):
                    continue
                if np.all(E_coeffs == 0):
                    continue
                E_symb = Add(*[ei * Di for ei, Di in zip(E_coeffs, basis)])

                # if self.geometry_context.divisor_data.evaluate(E_symb**2) != -2:
                #     continue

                E = LineBundle(E_symb, geometry_context=L.geometry_context)

                if self.tiltedSlope(E) <= 0:
                    ## E is in the torsion-free part F_B,ω, which means
                    ## it is really E[1] that we are looking at 
                    

                    if self.phase(E) <= self.phase(L):
                        ## Wont destabilize, skip
                        continue


                    if len(L_slope_HN) == 1 and E.slope >= L.slope:
                        ## If L is slope-stable and E has a larger slope, it
                        ## is a classical fact of slope stability that there 
                        ## cannot exist a morphism from E->L
                        continue


                    ## Use Riemann-Roch to see if Hom(E, L) is non-zero
                    diff_class = L.divisor - E_symb
                    RR_calc = 2 + 0.5*float(L.geometry_context.divisor_data.evaluate(diff_class**2).evalf())

                    if RR_calc > 0:
                        destabilizing_objects.append(E)

                else:
                    ## E is in the torsion part T_B,ω; the only way it can 
                    ## destabilize L is if there exists a morphism E -> L[1], 
                    ## which is equivalent to saying that H^0(X, L - E) is non-zero.

                    if self.phase(E) <= self.phase(L) + 1:
                        ## E -> L[1] does not actually destabilize
                        continue

                    ## Use Riemann-Roch to see if Ext^1(E, L) is non-zero
                    diff_class = L.divisor - E_symb
                    RR_calc = 2 + 0.5*float(L.geometry_context.divisor_data.evaluate(diff_class**2).evalf())

                    if RR_calc <= 0:
                        destabilizing_objects.append(E)
            
            
            print("Moving on to higher rank bundles")


            B = Add(*[coeff * basis_elem for coeff, basis_elem in zip(self.parameters[:-1], basis)] )
            W = self.parameters[-1] * self.geometry_context.polarization

            #####
            #
            # STEP (b) : Look for higher-rank vector bundles that may destabilize L
            #      ----------------------------------------
            #
            # The product function actually makes this an n-fold loop; specifically, its growth rate
            # is roughly O(n^r) where r is the Picard rank
            #####
            for r in range(2, r_max+1):

                for coeffs in product(*ranges):

                    E_c1_coeffs = np.array(coeffs)
                    if require_effective and np.any(E_c1_coeffs < 0):
                        continue
                    if np.all(E_c1_coeffs == 0):
                        continue
                    E_c1 = Add(*[ei * Di for ei, Di in zip(E_c1_coeffs, basis)])

                    ### Since E is not a DerivedCategoryObject yet, we cant just call twistedSlope()
                    E_twisted_slope = float(self.geometry_context.divisor_data.evaluate( (E_c1 - r*B), W).evalf())/r
                    
                    if E_twisted_slope <= 0:
                        ## E is in the torsion-free part F_B,ω, which means we are really looking at E[1].
                        ## Thus, this reduces to finding a standard morphism E->L. We can do some basic checks
                        ## to see if this is even possible.

                        E_std_slope_numerator = self.geometry_context.divisor_data.evaluate(E_c1 * self.geometry_context.polarization )
                        E_std_slope_denominator = r * self.geometry_context.divisor_data.evaluate(self.geometry_context.polarization, self.geometry_context.polarization)
                        E_std_slope = float(E_std_slope_numerator.evalf()) / float(E_std_slope_denominator.evalf())

                        if len(L_slope_HN) == 1 and E_std_slope >= L.slope:
                                ## If L is slope-stable and E has a larger slope, it
                                ## is a classical fact of slope stability that there 
                                ## cannot exist a morphism from E->L
                                continue

                        

                        ## The cone is a 2 term complex [A -> B] with H^{-1}(cone(f)) = ker(f) and H^0(cone(f)) = coker(f). This
                        ## is in the heart A_{B,ω} if and only if ker(f) is in the torsion free part F_B,ω and coker(f) is in the
                        ## torsion part T_B,ω.
                        ker_c1 = E_c1 - L.divisor

                        ker_tilt_slope =  float(self.geometry_context.divisor_data.evaluate( (ker_c1 - (r-1) * B), W).evalf())/ (r-1)

                        if ker_tilt_slope > 0:
                            ## The kernel is not in the torsion-free part F_B,ω 
                            continue

                        ## We now know that E -> L is a monomorphism, so lastly want to figure out what the second chern character is
                        ## We can use a Bogomolov-Gieseker inequality to bound the second chern character. In fact, a 
                        ## slightly stricter bound guaranteeing non-emptyness of moduli space is <v, v> >= -2

                        bound = (float(self.geometry_context.divisor_data.evaluate( E_c1**2 ).evalf() ) - 2*r**2 + 2) / 2*r

                        for ch2 in range(-1*int(bound), int(bound)+1):
                            mukai_vector_coeffs = [r] + list(coeffs) + [ch2]
                            if reduce(math.gcd, mukai_vector_coeffs) != 1:
                                ## In order for the moduli space to be nonempty, the mukai vector must be
                                ## primitive
                                continue


                            chi = -1*float(self.geometry_context.divisor_data.evaluate( E_c1 * L.divisor ).evalf()) \
                                + r*float(self.geometry_context.divisor_data.evaluate(L.divisor, L.divisor).evalf())/2 + 2*r + ch2
                            
                            if chi <=0:
                                ## Riemann Roch indicates there are no morphisms
                                continue


                            H2 = self.geometry_context.polarization**2
                            ch2 = ch2 / float(self.geometry_context.divisor_data.evaluate(H2).evalf())

                            chern = ChernCharacter(expr=r + E_c1 + ch2 * H2,
                                                basis=basis,
                                                dimension=2)
                            
                            destabilizing_sheaf = CoherentSheaf(chern, geometry_context=L.geometry_context) 
                            if self.phase(destabilizing_sheaf) > self.phase(L):
                                destabilizing_objects.append( destabilizing_sheaf )

                    else:
                        ## E is in the torsion part T_B,ω, so we are now searching for morphisms f: E -> L[1]
                        ## Notice that morphism f: E-> L[1] correspond precisly to Yoneda extension classes
                        ##  0 -> L -> A -> E -> 0, since in the derived category this short exact sequence 
                        ## becomes L -> A -> E -> L[1] with the last map being precisely f. As a consequence
                        ## the cone of f is quasi-isomorphic to the shift of the extension class A, i.e.
                        ##
                        ## cone(f) = A[1]
                        ##
                        ## Again, in order for E -> L[1] to be a monomorphism, the cone(f) must be in the heart A_{B,ω}
                        ## which is the same as saying A[1] must be in the heart A_{B,ω}. However, a shift of a coherent
                        ## sheaf is in the heart if and only if its part of the torsion-free part F_B,ω consisting
                        ## of sheaves with slope μ_Β,ω <= 0. 


                        A_c1 = E_c1 + L.divisor
                        A_tilt_slope =  float(self.geometry_context.divisor_data.evaluate( (A_c1 - (r+1) * B), W).evalf())/ (r+1)

                        if A_tilt_slope > 0:
                            ## The cone is not in the heart A_{B,ω} 
                            continue

                        bound = (float(self.geometry_context.divisor_data.evaluate( E_c1**2 ).evalf() ) - 2*r**2 + 2) / 2*r

                        for ch2 in range(-1*bound, bound+1):
                            mukai_vector_coeffs = [r] + list(coeffs) + [ch2]
                            if reduce(math.gcd, mukai_vector_coeffs) != 1:
                                ## In order for the moduli space to be nonempty, the mukai vector must be
                                ## primitive
                                continue


                            chi = -1*float(self.geometry_context.divisor_data.evaluate( E_c1 * L.divisor ).evalf()) \
                                + r*float(self.geometry_context.divisor_data.evaluate(L.divisor, L.divisor).evalf())/2 + 2*r + ch2
                            
                            if chi > 0:
                                ## Riemann Roch indicates there are no morphisms
                                continue


                            H2 = self.geometry_context.polarization**2
                            ch2 = ch2 / float(self.geometry_context.divisor_data.evaluate(H2).evalf())

                            chern = ChernCharacter(expr=r + E_c1 + ch2 * H2,
                                                basis=basis,
                                                dimension=2)
                            
                            destabilizing_sheaf = CoherentSheaf(chern, geometry_context=L.geometry_context) 
                            if self.phase(destabilizing_sheaf) > self.phase(L):
                                destabilizing_objects.append( destabilizing_sheaf )

            





        else: 
            #####################################
            # CASE 3: The Harder-Narasimhan factors lie in both the torsion part T_B,ω and the
            #           torsion-free part F_B,ω
            #####################################

            raise NotImplementedError("Need to implement this since technicially every object should fit into 0 ->T[1]-> E-> F ->0")

                    


    


        




        






    def in_tilted_heart(self, derived_obj : DerivedCategoryObject) -> bool:

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Expected a DerivedCategoryObject.")
        if not derived_obj.geometry_context == self.geometry_context:
            raise TypeError("Derived category object must be of the same geometry context as the stability condition")
    


        if isinstance(derived_obj, CoherentSheaf):
            slope_stab = SlopeStability(self.geometry_context)

            slope_hn = slope_stab.get_HN_factors(derived_obj)

            return all( self.tiltedSlope(factor.obj) > 0 for factor in slope_hn )
            
        elif isinstance(derived_obj, GradedCoproductObject):

            if len(derived_obj) == 1:

                if derived_obj.shift_vector[0] == 0:
                    ## We are just dealing with a direct sum of the same object
                    return self.in_tilted_heart(derived_obj.object_vector[0])
                elif derived_obj.shift_vector[0] == 1 and isinstance(derived_obj.object_vector[0], CoherentSheaf):
                    
                    slope_hn = self.get_HN_factors(derived_obj.object_vector[0])
                    return all( self.tiltedSlope(factor.obj) <= 0 for factor in slope_hn )
                
            else:
                ####
                # TODO: Do something to adderss the fact that <F[1], T> is different from the union of F[1] and T;
                #         we should be checking extensions somehow
                ####
                return False


   





