from src.DerivedCategory.GeometryContext import GeometryContext
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject, NumericalObject
from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.CoherentSheaf import LineBundle
from . import HarderNarasimhanError, HarderNarasimhanFiltration, StabilityCondition
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition

from dotenv import load_dotenv
import os

from sympy import symbols, lambdify

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
            case 'P2' | 'LocalP2':
                if parameter_list.shape[1] != 2:
                    raise ValueError("For P2, parameter_list must have shape (N, 2) for real parameters (s, q).")
                if not np.issubdtype(parameter_list.dtype, np.floating):
                    raise TypeError("For P2, parameters must be real numbers (float or int).")
            case 'K3':
                if parameter_list.shape[1] != 2:
                    raise ValueError("For K3, parameter_list must have shape (N, 2) for (α, β).")
                if not np.issubdtype(parameter_list.dtype, np.floating):
                    raise TypeError("For K3, parameters must be real numbers (float or int).")
            case _:
                raise NotImplementedError(f"Stability condition family not implemented for {geometry_context.catagory}")

        self.geometry_context = geometry_context
        self.parameter_list = parameter_list


    def centralCharge(self, derived_obj : DerivedCategoryObject) -> np.ndarray:

        if not isinstance(derived_obj, DerivedCategoryObject):
            raise TypeError("Derived category object must be of type DerivedCategoryObject")

        chern_character = derived_obj.chernCharacter()
        polarization = self.geometry_context.polarization

        if self.geometry_context.catagory == 'P1' or self.geometry_context.catagory == 'LocalP1':

            param_character = ChernCharacter(expr= -1 + self.parameters[0] * polarization, basis=[polarization], dimension=1)

            return self.geometry_context.divisor_data.evaluate((chern_character*param_character)[1])
    