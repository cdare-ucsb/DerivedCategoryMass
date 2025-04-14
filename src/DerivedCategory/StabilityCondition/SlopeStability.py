 
from src.DerivedCategory.CoherentSheaf import CoherentSheaf, LineBundle
from src.DerivedCategory.GeometryContext import GeometryContext
from src.DerivedCategory.StabilityCondition import HarderNarasimhanFiltration

from typing import List

import numpy as np
from sympy import Add

from itertools import product




class SlopeStability:

    def __init__(self, geometry_context):

        if not isinstance(geometry_context, GeometryContext):
            raise TypeError("Expected a GeometryContext.")
        self.geometry_context = geometry_context
        self._hn_filt_cache = {}


    

    def line_bundle_destabilizers(self, L : LineBundle, require_effective : bool = True) -> List[LineBundle]:

        if not isinstance(L, LineBundle):
            raise TypeError("Expected a LineBundle.")
        if not L.geometry_context == self.geometry_context:
            raise TypeError("Line bundle must be of the same geometry context as the stability condition")
        

        if self.geometry_context.catagory == "K3":

            basis = L.geometry_context.divisor_data.basis

            ###
            # Step 1: Convert lb to a list of ints representing the coefficients of the basis
            ###

            L_coeffs = [int(L.chernCharacter()[b]) for b in basis]

            ####
            # Step 2: Create bounding box to look for (-2)-candidates
            ####
            ranges = [range(di + 1) for di in L_coeffs]  # assume effective basis and E <= D

            destabilizers = []



            for coeffs in product(*ranges):
                E_coeffs = np.array(coeffs)
                
                
                if require_effective and np.any(E_coeffs < 0):
                    continue
                if np.all(E_coeffs == 0):
                    continue
                E_symb = Add(*[ei * Di for ei, Di in zip(E_coeffs, basis)])

                # if self.geometry_context.divisor_data.evaluate(E_symb**2) != -2:
                #     continue

                E = LineBundle(E_symb, geometry_context=self.geometry_context)

                

                if E.slope <=  L.slope:
                    continue

                diff_class = L.divisor - E_symb

                RR_calc = 2 + 0.5*float(self.geometry_context.divisor_data.evaluate(diff_class**2).evalf())
                if RR_calc > 0:
                    # print(f"{E} slope: {E.slope}, morph exist")
                    destabilizers.append(E)
                

            return destabilizers



    def phase(self, coh : CoherentSheaf) -> float:

        if not isinstance(coh, CoherentSheaf):
            raise TypeError("Expected a CoherentSheaf.")
        if not coh.geometry_context == self.geometry_context:
            raise TypeError("Coherent sheaf must be of the same geometry context as the stability condition")

        if coh.rank == 0:
            return 1
        else:
             
            n = self.geometry_context.divisor_data.variety_dimension
             
            degree = self.geometry_context.divisor_data.evaluate(coh.c1, (self.geometry_context.polarization)**(n-1))

            arg = np.arctan2(coh.rank, -1*float(degree.evalf()) )
            if arg <= 0:
                arg += 2 * np.pi  # shift to (0, 2π]
            
            return arg / np.pi  # normalize to (0, 2]




                
    def get_HN_factors(self, coh : CoherentSheaf) -> HarderNarasimhanFiltration:

        if not isinstance(coh, CoherentSheaf):
            raise TypeError("Expected a CoherentSheaf.")
        if not coh.geometry_context == self.geometry_context:
            raise TypeError("Coherent sheaf must be of the same geometry context as the stability condition")
        if coh in self._hn_filt_cache:
            return self._hn_filt_cache[coh]



        if isinstance(coh, LineBundle):

            destabilizers = list(self.line_bundle_destabilizers(coh, require_effective=True))
            
            if len(destabilizers) == 0:
                self._hn_filt_cache[coh] = HarderNarasimhanFiltration(stable_objects=[coh],
                                            phase_vector=[self.phase(coh)])
            else:
                ## We must find the maximal destabilizing subobject

                max_destabilizer = max(destabilizers, key=lambda E: E.slope)

                quotient = CoherentSheaf(chern_character= ( coh.chernCharacter() - max_destabilizer.chernCharacter()),
                                        geometry_context=self.geometry_context)

                quotient_phase = self.phase(quotient)

                self._hn_filt_cache[coh] = HarderNarasimhanFiltration(stable_objects=[max_destabilizer, quotient],
                                            phase_vector=[self.phase(coh=coh), quotient_phase])
                
            return self._hn_filt_cache[coh]
