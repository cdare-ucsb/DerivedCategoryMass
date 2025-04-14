
from src.DerivedCategory.DerivedCategoryObject import DerivedCategoryObject, GradedCoproductObject
from typing import List, Tuple
from dataclasses import dataclass


class HarderNarasimhanError(Exception):
    r"""!
    Exception raised when the correct Harder-Narasimhan filtration cannot be found 
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args)

        self.message = kwargs.get('message') ## The error message

        self.stability_parameters = kwargs.get('stability_parameters') ## The stability parameters used to compute the Harder-Narasimhan factors


@dataclass
class HNFactor:
    """
    Represents a single semistable factor in a Harder–Narasimhan filtration.

    \param phase: The phase (argument of central charge / π)
    \param obj: A semistable derived category object
    \param multiplicity: The number of times this object appears in the filtration
    """
    phase: float
    obj: DerivedCategoryObject
    multiplicity: int = 1
    original_shift: int = 0

    def __repr__(self):
        base = f"{self.obj}"
        if self.original_shift != 0:
            base += f"[{self.original_shift}]"
        
        if self.multiplicity > 1:
            base = f"({base})^{self.multiplicity}"
        
        return f"{base} ∈ P({self.phase:.3f})"
    
    def __lt__(self, other):
        if not isinstance(other, HNFactor):
            return NotImplemented
        return self.phase < other.phase

    def __le__(self, other):
        if not isinstance(other, HNFactor):
            return NotImplemented
        return self.phase <= other.phase

    def __gt__(self, other):
        if not isinstance(other, HNFactor):
            return NotImplemented
        return self.phase > other.phase

    def __ge__(self, other):
        if not isinstance(other, HNFactor):
            return NotImplemented
        return self.phase >= other.phase


class HarderNarasimhanFiltration():

    def __init__(self, stable_objects : List[DerivedCategoryObject], phase_vector : List[float], dimension_vector : List[int] = None, shift_vector : List[int] = None):

        
        ##################
        # Input Validation
        ##################


        # stable objects validation
        if not isinstance(stable_objects, list):
            raise TypeError("The stable objects must be a list of DerivedCategoryObject instances.")

        if stable_objects:
            first_obj_geometry_context = stable_objects[0].geometry_context
            if not all(obj.geometry_context == first_obj_geometry_context for obj in stable_objects):
                raise ValueError("All stable objects must have the same geometry context.")
            self.geometry_context = first_obj_geometry_context
        else:
            self.geometry_context = None  # Handle trivial filtration gracefully
        

        # phase vector validation
        if not isinstance(phase_vector, list):
            raise TypeError("The phase vector must be a list of floats.")
        if len(phase_vector) != len(stable_objects):
            raise ValueError("The phase vector must have the same length as the number of stable objects.")
        if not all(isinstance(p, (int, float)) for p in phase_vector):
            raise TypeError("The phase vector must contain only integers or floats.")
        
        # dimension vector validation
        if dimension_vector is None:
            dimension_vector = [1] * len(stable_objects)
        else:
            if not isinstance(dimension_vector, list):
                raise TypeError("The dimension vector must be a list of integers.")
                
            if len(dimension_vector) != len(stable_objects):
                raise ValueError("The dimension vector must have the same length as the number of stable objects.")
            if not all(isinstance(d, int) for d in dimension_vector):
                raise TypeError("The dimension vector must contain only integers.")
            if not all(d >= 0 for d in dimension_vector):
                raise ValueError("The dimension vector must contain only non-negative integers.")
            
        if shift_vector is None:
            shift_vector = [0] * len(stable_objects)
        else:
            if not isinstance(shift_vector, list):
                raise TypeError("The shift vector must be a list of integers.")
            if len(shift_vector) != len(stable_objects):
                raise ValueError("The shift vector must have the same length as the number of stable objects.")
            if not all(isinstance(s, int) for s in shift_vector):
                raise TypeError("The shift vector must contain only integers.")
            
        # assign the vectors to member variables

        self.factors = sorted(
            [HNFactor(phase=phase, obj=obj, multiplicity=dim, original_shift=shift) for obj, phase, dim, shift in zip(stable_objects, phase_vector, dimension_vector, shift_vector) if dim != 0],
            reverse=True
        )

    def shift(self, n: int) -> 'HarderNarasimhanFiltration':
        r"""!
        Returns a new Harder-Narasimhan filtration with all phases shifted by an integer amount.

        \param n: The integer amount to shift the phases by
        \return: A new HarderNarasimhanFiltration with shifted phases
        """
        if not isinstance(n, int):
            raise TypeError("Shift amount must be an integer.")

        return HarderNarasimhanFiltration(
            stable_objects=[f.obj for f in self.factors],
            phase_vector=[f.phase + n for f in self.factors],
            dimension_vector=[f.multiplicity for f in self.factors],
            shift_vector=[f.original_shift + n for f in self.factors]
        )


    def splitAtPhase(self, phase_cut: float) -> Tuple['HarderNarasimhanFiltration', 'HarderNarasimhanFiltration']:
        r"""!
        Splits the Harder–Narasimhan filtration at the specified phase value.

        \param phase_cut: The phase value to split at
        \return: A tuple (HN≥, HN<) where:
                - the first component contains all HNFactors with phase ≥ phase_cut
                - the second component contains all HNFactors with phase < phase_cut
        """
        if not isinstance(phase_cut, (float, int)):
            raise TypeError("Phase cut must be a float or int.")

        above_or_equal = [f for f in self.factors if f.phase >= phase_cut]
        below = [f for f in self.factors if f.phase < phase_cut]

        return (
            HarderNarasimhanFiltration(
                stable_objects=[f.obj for f in above_or_equal],
                phase_vector=[f.phase for f in above_or_equal],
                dimension_vector=[f.multiplicity for f in above_or_equal],
                shift_vector=[f.original_shift for f in above_or_equal],
            ),
            HarderNarasimhanFiltration(
                stable_objects=[f.obj for f in below],
                phase_vector=[f.phase for f in below],
                dimension_vector=[f.multiplicity for f in below],
                shift_vector=[f.original_shift for f in below],
            ),
        )

    def toGradedCoproductObject(self) -> GradedCoproductObject:

        stable_objects=[f.obj for f in self.factors],
        dimension_vector=[f.multiplicity for f in self.factors],
        shift_vector=[f.original_shift for f in self.factors]

        return GradedCoproductObject(
            object_vector=stable_objects,
            dimension_vector=dimension_vector,
            shift_vector=shift_vector)
    


    def __repr__(self):
        """
        Returns a multi-line string representation of the Harder-Narasimhan filtration,
        with each HN factor on its own line, separated by a horizontal rule.
        """
        if not self.factors:
            return "HarderNarasimhanFiltration([])"
        
        separator = "\n\n------------------\n\n"
        return separator.join(str(factor) for factor in self.factors)



    def __add__(self, other):
        if isinstance(other, HarderNarasimhanFiltration):
            from collections import defaultdict

            # Phase → {(obj, shift) → multiplicity}
            merged = defaultdict(lambda: defaultdict(int))

            for factor in self.factors + other.factors:
                merged[factor.phase][(factor.obj, factor.original_shift)] += factor.multiplicity

            combined = []
            for phase in sorted(merged.keys(), reverse=True):
                for (obj, shift), mult in merged[phase].items():
                    combined.append(HNFactor(phase, obj, mult, original_shift=shift))

            return HarderNarasimhanFiltration(
                stable_objects=[f.obj for f in combined],
                phase_vector=[f.phase for f in combined],
                dimension_vector=[f.multiplicity for f in combined],
                shift_vector=[f.original_shift for f in combined]
            )
        else:
            return NotImplemented

    def __mul__(self, n: int) -> 'HarderNarasimhanFiltration':
        """
        Scale the filtration by a positive integer, multiplying all multiplicities.

        \param n: The integer multiplier (must be ≥ 1)
        \return: A new HarderNarasimhanFiltration with scaled multiplicities
        """
        if not isinstance(n, int) or n < 1:
            raise ValueError("HN filtrations can only be multiplied by positive integers.")

        return HarderNarasimhanFiltration(
            stable_objects=[f.obj for f in self.factors],
            phase_vector=[f.phase for f in self.factors],
            dimension_vector=[f.multiplicity * n for f in self.factors],
            shift_vector=[f.original_shift for f in self.factors]
        )
    
    def __rmul__(self, n: int) -> 'HarderNarasimhanFiltration':
        r"""!
        Scale the filtration by a positive integer, multiplying all multiplicities.

        \param n: The integer multiplier (must be ≥ 1)
        \return: A new HarderNarasimhanFiltration with scaled multiplicities
        """
        return self.__mul__(n)


    def __len__(self):
        r"""!
        Returns the number of objects in the Harder-Narasimhan filtration.
        """
        return len(self.factors)
    
    
    def __iter__(self):
        r"""!
        Returns an iterator over the HN factors in descending phase order.
        """
        return iter(self.factors)
    
    def __getitem__(self, idx):
        if idx < 0:
            idx += len(self.factors)
        return self.factors[idx]
        
