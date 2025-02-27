from .Divisor import  Divisor
import numpy as np


class NeronSeveriGroup:

    def __init__(self, basis_divisors, intersection_matrix):
        """
        Initialize an instance of NeronSeveriGroup with the specified basis divisors and intersection matrix.

        Parameters:
        ----------
        basis_divisors : list of Divisor
            The basis divisors for the Neron-Severi group
        intersection_matrix : np.ndarray
            The intersection matrix of the divisors

        Raises:
        -------
        TypeError
            If any of the basis divisors are not instances of Divisor
        ValueError
            If the intersection matrix is not square or if the number of basis divisors does not match the size of the intersection matrix
        """

        if not all(isinstance(x, Divisor) for x in basis_divisors):
            raise TypeError("All basis divisors must be instances of Divisor")

        if intersection_matrix.shape[0] != intersection_matrix.shape[1]:
            raise ValueError("Intersection matrix must be square")

        if intersection_matrix.shape[0] != len(basis_divisors):
            raise ValueError("Number of basis divisors must match the size of the intersection matrix")

        self.basis_divisors = basis_divisors
        self.intersection_matrix = intersection_matrix