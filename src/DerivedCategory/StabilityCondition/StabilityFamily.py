from src.DerivedCategory.GeometryContext import GeometryContext
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject, NumericalObject
from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.CoherentSheaf import LineBundle
from . import HarderNarasimhanError, HarderNarasimhanFiltration, StabilityCondition
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition

from dotenv import load_dotenv
import os

from sympy import symbols, lambdify, I, expand, exp

import numpy as np
import math
import cmath

# Load .env file
load_dotenv()
IMPLEMENTED_CATAGORIES = os.getenv("IMPLEMENTED_CATAGORIES").split(",") # ['LocalP1', 'P1', 'LocalP2', 'P2', 'K3']
__CURRENT_DOUBLE_TWIST_IMPLEMENTED__ = os.getenv("CURRENT_DOUBLE_TWIST_IMPLEMENTED").split(",") # ['K3']



class StabilityFamily():

    def __init__(self, geometry_context: GeometryContext, parameter_list : np.ndarray ):
        
        if not isinstance(geometry_context, GeometryContext):
            raise TypeError("geometry_context must be an instance of GeometryContext")
        if not isinstance(parameter_list, np.ndarray):
            raise TypeError("parameter_list must be a numpy array")
        
        # Ensure it's a 2D array: (N, num_parameters)
        if parameter_list.ndim != 2:
            raise ValueError("parameter_list must be a 2D numpy array of shape (N, d)")
        
        match geometry_context.catagory:
            case 'P1' | 'LocalP1':
                if parameter_list.shape[1] != 1:
                    raise ValueError("For P1, parameter_list must have shape (N, 1) for N complex parameters.")
                if not np.issubdtype(parameter_list.dtype, np.complexfloating):
                    raise TypeError("For P1, parameters must be complex numbers.")
                
                self._init_p1_central_charge_lambdify()

            case 'P2' | 'LocalP2':
                if parameter_list.shape[1] != 2:
                    raise ValueError("For P2, parameter_list must have shape (N, 2) for real parameters (s, q).")
                if not np.issubdtype(parameter_list.dtype, np.floating):
                    raise TypeError("For P2, parameters must be real numbers (float or int).")
                
                self._init_p2_central_charge_lambdify()

            case 'K3':
                if parameter_list.shape[1] != 2:
                    raise ValueError("For K3, parameter_list must have shape (N, 2) for (α, β).")
                if not np.issubdtype(parameter_list.dtype, np.floating):
                    raise TypeError("For K3, parameters must be real numbers (float or int).")
                
                self._init_k3_central_charge_lambdify()

            case _:
                raise NotImplementedError(f"Stability condition family not implemented for {geometry_context.catagory}")

        self.geometry_context = geometry_context
        self.parameter_list = parameter_list



    def centralCharge(self, derived_obj : DerivedCategoryObject) -> np.ndarray:

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Derived category object must be of type DerivedCategoryObject")

        ch = derived_obj.chernCharacter()
        polarization = self.geometry_context.polarization

        match self.geometry_context.catagory:
            case 'P1' | 'LocalP1':

                ch0 = float(ch[0])  # 1
                ch1 = float(ch[polarization])  # coefficient of H

                w_array = self.parameter_list[:, 0]  # shape (N,)
                return self._Z_p1_func(w_array, ch0, ch1)
            
            case 'P2' | 'LocalP2':
                ch0 = float(ch[0])  # rank
                ch1 = float(ch[polarization])  # linear term (e.g., ch1 H)
                ch2 = float(ch[polarization**2])  # reduce top term to scalar

                # Extract (s, q) parameters from parameter grid
                s_array = self.parameter_list[:, 0]
                q_array = self.parameter_list[:, 1]

                # Evaluate lambdified central charge
                return self._Z_p2_func(s_array, q_array, ch0, ch1, ch2)
            
            case 'K3':
                # Step 2: extract ch0, ch1_i, ch2
                ch0 = float(ch[0])  # scalar rank
                # Can represent ch2 as a scalar since its already a top-dimensional class
                ch2 = float(self.geometry_context.divisor_data.evaluate(ch[2]))  

                basis = self.geometry_context.divisor_data.basis
                ch1_vec = [float(ch[D]) for D in basis]  # ch1 = sum ch1_i * D_i

                # Step 3: extract parameters from grid (N, r+1): b0...br, omega
                param_array = self.parameter_list
                b_grid = param_array[:, :-1].T  # (r, N)
                omega_grid = param_array[:, -1]  # (N,)

                # Step 4: apply the lambdified function Z_k3_func
                # Inputs = (b0...br, omega, ch0, ch1_0...ch1_r, ch2)
                return self._Z_k3_func(*b_grid, omega_grid, ch0, *ch1_vec, ch2)
            
            case _:
                raise NotImplementedError(f"Stability condition family not implemented for {self.geometry_context.catagory}")


    













    ####################################################
    #                                                  #
    #          Central Charge Lambdifications          #
    #                                                  #
    ####################################################


    def _init_p1_central_charge_lambdify(self):
        """Precompute the symbolic expression and lambdified function for P1 central charge."""
        w = symbols('w')
        H = self.geometry_context.polarization
        ch0, ch1 = symbols("ch0 ch1")  # Assume Chern character in basis 1 + H
        # Central charge is ∫ (-1 + w H) · (ch0 + ch1 H) = -ch1 + w ch0
        Z_expr = -1 * ch1 + w * ch0

        # Lambdify Z(w) as a function of w, ch0, ch1
        self._Z_p1_func = lambdify((w, ch0, ch1), Z_expr, modules='numpy')


    def _init_p2_central_charge_lambdify(self):
        """
        Precompute a lambdified function to evaluate the central charge on P2:
        Z(s, q, ch0, ch1, ch2) = (-ch2 + q ch0) + i (ch1 - s ch0)
        """
        from sympy import symbols, I, lambdify

        # Define symbolic variables
        s, q = symbols("s q")
        ch0, ch1, ch2 = symbols("ch0 ch1 ch2")

        # Construct the symbolic central charge expression
        Z_expr = (-ch2 + q * ch0) + I * (ch1 - s * ch0)

        # Lambdify
        self._Z_p2_func = lambdify((s, q, ch0, ch1, ch2), Z_expr, modules="numpy")


    def _init_k3_central_charge_lambdify(self):
        basis = self.geometry_context.divisor_data.basis
        intersection = self.geometry_context.divisor_data.top_intersection_form

        r = len(basis)
        H = self.geometry_context.polarization

        # Parameters: b_i and omega
        b_syms = symbols(f"b0:{r}")
        omega_sym = symbols("omega")

        # B = sum b_i * D_i
        B = sum(b * D for b, D in zip(b_syms, basis))
        W = omega_sym * H
        twist = -1 * (B + I * W)

        # Build exp(-B + iωH)
        param_expr = expand(exp(twist).series(n=3).removeO())

        # Dummy ch(E) input symbols
        ch0 = symbols("ch0")
        ch1_syms = symbols(f"ch1_0:{r}")
        ch2 = symbols("ch2")

        # Build ch1 = sum ch1_i * D_i
        ch1_expr = sum(ci * Di for ci, Di in zip(ch1_syms, basis))
        ch_expr = expand(ch0 + ch1_expr + ch2 * H**2)  # H^2 placeholder, will integrate away

        # Pair: integrate degree 2 part of exp * ch
        product = expand(param_expr * ch_expr)


        deg2_expr = 0
        for term in product.as_ordered_terms():
            coeff, mon = term.as_coeff_Mul()

            # Match monomial part: look for degree-2 expressions like D_i * D_j
            if mon == 1:
                continue

            try:
                poly = mon.as_poly(*basis)
                if poly.total_degree() != 2:
                    continue

                for m in poly.monoms():  # monoms() → list of (deg_0, ..., deg_r)
                    if sum(m) != 2:
                        continue

                    # Identify which two divisors this monomial represents
                    idxs = [i for i, deg in enumerate(m) if deg > 0]
                    if len(idxs) == 1:
                        # diagonal term like D_i^2
                        i = idxs[0]
                        div1 = div2 = basis[i]
                        power = m[i]
                        if power != 2:
                            continue
                    elif len(idxs) == 2:
                        i, j = idxs
                        div1 = basis[i]
                        div2 = basis[j]
                    else:
                        continue  # should not happen for degree 2 monomials

                    # Use intersection form to reduce D_i D_j to a number
                    value = intersection.get((div1, div2), 0)
                    deg2_expr += coeff * poly.coeff_monomial(m) * value

            except Exception:
                continue


        # Final central charge expression: - (deg2 + ch0)
        Z_expr = -1 * (deg2_expr + ch0)

        # Lambdify
        all_syms = (*b_syms, omega_sym, ch0, *ch1_syms, ch2)
        self._Z_k3_func = lambdify(all_syms, Z_expr, modules="numpy")
