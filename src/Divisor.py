


class Divisor:
    
    def __init__(self, degree):
        """
        Initialize an instance of Divisor with the specified degree.
        
        Parameters:
        ----------
        degree : int
            The degree of the divisor
            
        Raises:
        -------
        TypeError
            If the degree is not an integer
        """
        
        if not isinstance(degree, int):
            raise TypeError("Degree must be an integer")
        
        self.degree = degree