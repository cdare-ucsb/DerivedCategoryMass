from ChernCharacter import ChernCharacterP2, ChernCharacterP1
import math
import cmath


###############################################################################
#                                                                             #
#                            Single-Degree Objects                            #
# ----------------------------------------------------------------------------#
#  These objects are used to represent coherent sheaves, vector bundles, and  #
#  line bundles. Line bundles (i.e. vector bundles of rank 1) will be our     #
#  building blocks for the sake of this project, since every coherent sheaf   #
#  in projective space  has a resolution by line bundles. Moreover, all       #         
#  spherical twists that we will consider can be represented as a composition #
#  of twists around (the pushforward of) line bundles.                        #    
#                                                                             #
###############################################################################





####################################
#      Generic Parent Classes      #
####################################

class CoherentSheaf():
    
    def __init__(self):
        pass

    def chernCharacter(self):
        pass


class VectorBundle(CoherentSheaf):

    def __init__(self):
        pass


class LineBundle(VectorBundle):

    def __init__(self):
        pass





class CoherentSheafP1(CoherentSheaf):
    """
    This class represents a general coherent sheaf on the projective line. 
    Coherent sheaves are generalizations of vector bundles since they can have torsion.
    Since we will assume that the base field is the complex numbers, we can assume that
    all coherent sheaves are simple represented by C[x,y]-modules that are localized in
    the sense that they agree on the overlap of the distinguished affine open sets (e.g.
    U_x = {x != 0}, U_y = {y != 0})

    While coherent sheaves are not uniquely determined by their characteristic classes, 
    we will use Chern classes as the minimal amount of information to represent a coherent sheaf.


    """

    def __init__(self, rank, c1):
        self.c0 = int(rank)
        self.c1 = int(c1)

    def chernCharacter(self):
        return ChernCharacterP1(self.c0, self.c1)
    
    def __str__(self):
        """
        String representation of the coherent sheaf class

        Returns:
        -------
        str
            A string representation of the coherent sheaf
        """
        return f'CoherentSheaf of rank {self.c0} and degree {self.c1}'
    
    def central_charge(self, w):
        """
        Method to compute the central charge of the coherent sheaf. The central charge is a complex-valued
        homomorphism on the Grothendieck group of coherent sheaves, and is defined for local P1 via

        Zw = -ch1 + w ch0

        The central charge is used to compute the phase of the central charge, which is used to
        determine the stability of the sheaf.

        Parameters:
        ----------
        w : complex
            The complex number parameterizing the central charge

        Returns:
        -------
        complex
            The central charge of the coherent sheaf as a complex number

        Raises:
        -------
        TypeError
            If w is not a complex number

        """
        
        if not isinstance(w, complex):
            raise TypeError("w must be a complex number")
    

        return -1*self.c1 + w*self.c0
    
    def phase(self, w):
        """
        The phase, in terms of stability conditions, is actually considered as the argument
        divided by pi. This is used to determine the stability of the sheaf with regard to its 
        subobjects

        Parameters:
        ----------
        w : complex
            The complex parameter for the central charge
        

        Returns:
        -------
        float
            The phase of the central charge, divided by pi

        Raises:
        -------
        TypeError
            If w is not a complex number
        """


        if not isinstance(w, complex):
            raise TypeError("w must be a complex number")
        
        return cmath.phase(self.central_charge(w)) / math.pi


class VectorBundleP1(VectorBundle, CoherentSheafP1):

    def __init__(self, rank, deg):
        self.c0 = int(rank)
        self.c1 = int(deg)
    
    def chernCharacter(self):
        return ChernCharacterP1(self.c0, self.c1)

    def __str__(self):
        """
        String representation of the vector bundle class

        Returns:
        -------
        str
            A string representation of the vector bundle
        """
        return f'VectorBundle of rank {self.c0} and degree {self.c1}'
    



class LineBundleP1(LineBundle, CoherentSheafP1):

    def __init__(self, deg):
        self.c0 = 1
        self.c1 = deg


    def chernCharacter(self):
        return ChernCharacterP1(self.c0, self.c1)
    
    def __str__(self):
        """
        String representation of the line bundle class; the standard notation of a line bundle
        of degree d is O(d), since every line bundle is a twist of the structure sheaf.

        Returns:
        -------
        str
            A string representation of the line bundle
        """
        return f'O({self.c1})'

    





class CoherentSheafP2(CoherentSheaf):
    '''
    This class represents a general coherent sheaf on a projective plane. 
    Coherent sheaves are generalizations of vector bundles since they can have torsion.
    Since we will assume that the base field is the complex numbers, we can assume that
    all coherent sheaves are simple represented by C[x,y,z]-modules that are localized 
    in the sense that they agree on the overlap of distinguished affine open sets (e.g.
    U_x = {x != 0}, U_y = {y != 0}, U_z = {z != 0}).

    While coherent sheaves are not uniquely determined by their characteristic classes, 
    we will use Chern classes as the minimal amount of information to represent a coherent sheaf.


    Attributes:
    ----------
    rank : int
        The rank of the coherent sheaf
    c1 : int
        The degree of the coherent sheaf
    c2 : float
        The second Chern class of the coherent sheaf
    '''

    def __init__(self, rank, c1, c2 = 0):
        """
        Initialize an instance of CoherentSheaf with the specified chern
        classes.

        Parameters:
        ----------
        rank : int
            The rank of the coherent sheaf
        c1 : int
            The degree of the coherent sheaf
        c2 : float
            The second Chern class of the coherent sheaf
        """
        self.c0 = int(rank)
        self.c1 = int(c1)
        self.c2 = c2
    
    def chernCharacter(self):
        """
        Method which returns a Chern Characteristic class of the coherent sheaf
        The chern character is defined on line bundles as follows:

        ch(L) = exp(c1(L)) = 1 + c1(L) + c1^2(L) / 2 + ...

        However, for vector bundles the chern character will generally be written as

        ch(E) = rank(E) + c1(E) + (c_1^2(E) - c2(E)) / 2 + (c_1^3(E) - 3c1(E)c2(E) + c3(E)) / 6 + ...

        The numerics of the Chern Character are handled by a separate class, so that this method
        only returns a specific chern character.

        Returns:
        -------
        ChernCharacter
            The Chern character of the coherent sheaf
        
        """
        ch2 = float( self.c1**2 - 2 * self.c2 ) / 2 # The second Chern class of the sheaf
        return ChernCharacterP2(self.c0, self.c1, ch2)

    def __str__(self):
        """
        String representation of the coherent sheaf class

        Returns:
        -------
        str
            A string representation of the coherent sheaf
        """
        return f'CoherentSheaf of rank {self.c0}, degree {self.c1}, and c2 {self.c2}'
    
    def central_charge(self, s, q):
        """
        Method to compute the central charge of the coherent sheaf. The central charge is a complex-valued
        homomorphism on the Grothendieck group of coherent sheaves, and is defined for local P2 via

        Z(s, q) = -ch2 + q * ch0 + i(ch1 - s * ch0)

        where s and q are real numbers (c.f. Chunyi Li, The Space of Stability Conditions on the Projective
        Plane). The central charge is used to compute the phase of the central charge, which is used to
        determine the stability of the sheaf.

        Parameters:
        ----------
        s : float
            The parameter controlling the imaginary part of the central charge
        q : float
            The parameter controlling the real part of the central charge

        Returns:
        -------
        complex
            The central charge of the coherent sheaf as a complex number

        Raises:
        -------
        TypeError
            If s or q are not floating-point decimals

        """
        
        if not isinstance(s, float) or not isinstance(q, float):
            raise TypeError("s and q must be floating-point decimals.")
    
        chern_char = self.chernCharacter()
        return complex(-chern_char.ch2 + q * chern_char.ch0, chern_char.ch1 - s * chern_char.ch0)

    def phase(self, s, q):
        """
        The phase, in terms of stability conditions, is actually considered as the argument
        divided by pi. This is used to determine the stability of the sheaf with regard to its 
        subobjects

        Parameters:
        ----------
        s : float
            The parameter controlling the imaginary part of the central charge
        q : float
            The parameter controlling the real part of the central charge

        Returns:
        -------
        float
            The phase of the central charge, divided by pi

        Raises:
        -------
        TypeError
            If s or q are not floating-point
        """


        if not isinstance(s, float) or not isinstance(q, float):
            raise TypeError("s and q must be floating-point decimals.")
        
        return cmath.phase(self.central_charge(s, q)) / math.pi
    


class VectorBundleP2(CoherentSheafP2, VectorBundle):
    """
    A vector bundle is a special type of coherent sheaf; in particular, it is a locally free
    sheaf of finite rank. These will ultimately be better generalizations of line bundles, and
    any direct sum of line bundles is an example of a vector bundle. Other notable vector bundles
    that we will use are the tangent and cotangent bundles, which are (somewhat) canonical rank
    2 vector bundles. There is not much additional information to add to the class, since the
    rank, degree, and second Chern class are sufficient to determine the Chern Character. Moreover,
    any additional information would be redundant since we are not keeping track of torsion.

    Attributes:
    ----------
    rank : int
        The rank of the vector bundle
    deg : int
        The degree of the vector bundle
    c2 : float
        The second Chern class of the vector bundle
    """

    def __init__(self, rank, deg, c2):
        """
        Initialize an instance of VectorBundle with the specified chern
        classes.

        Parameters:
        ----------
        rank : int
            The rank of the vector bundle
        deg : int
            The degree of the vector bundle
        c2 : float
            The second Chern class of the vector bundle
        """
        self.rank = int(rank)
        self.deg = int(deg)
        self.c2 = float(c2)

    def __str__(self):
        """
        String representation of the vector bundle class

        Returns:
        -------
        str
            A string representation of the vector bundle
        """
        return f'VectorBundle of rank {self.c0}, degree {self.c1}, and c2 {self.c2}'
    
    def chernCharacter(self):
        """
        Method which returns a Chern Characteristic class of the vector bundle
        The chern character is defined on line bundles as follows:

        ch(L) = exp(c1(L)) = 1 + c1(L) + c1^2(L) / 2 + ...

        However, for vector bundles the chern character will generally be written as

        ch(E) = rank(E) + c1(E) + (c_1^2(E) - c2(E)) / 2 + (c_1^3(E) - 3c1(E)c2(E) + c3(E)) / 6 + ...

        The numerics of the Chern Character are handled by a separate class, so that this method
        only returns a specific chern character.

        Returns:
        -------
        ChernCharacter
            The Chern character of the vector bundle
        
        """
        return super().chernCharacter()
    
    def isLineBundleP2Sum(self):
        """
        Helper function to determine if the vector bundle is a sum of line bundles. The function 
        simpily calls the isLineBundleSum() method of the Chern Character class.

        Returns:
        -------
        bool
            True if the vector bundle is a sum of line bundles, False otherwise
        """
        return self.chernCharacter().isLineBundleP2Sum()
    

    def isCotangentBundleP2Sum(self):
        """
        Helper function to determine if the vector bundle is a sum of cotangent bundles. The function
        simply calls the isCotangentBundleSum() method of the Chern Character class.

        Returns:
        -------
        bool
            True if the vector bundle is a sum of cotangent bundles, False otherwise
        """
        return self.chernCharacter().isCotangentBundleP2Sum()



class LineBundleP2(VectorBundleP2, LineBundle):
    """
    A line bundle is a special type of vector bundle; in particular, it is a rank 1 vector bundle.
    These will ultimately be our building blocks for the sake of this project, since every coherent
    sheaf on P^2 has a resolution by line bundles. Moreover, all spherical twists that we will consider
    can be represented as a composition of twists around (the pushforward of) line bundles.

    As the Picard group of P^2 is isomorphic to Z, we can represent any line bundle by its degree.

    Attributes:
    ----------
    c0 : int
        The rank of the line bundle will always be 1
    c1 : int
        The degree of the line bundle
    c2 : float
        The second Chern class of the line bundle, which will always be of the form d^2 / 2
    """
    
    def __init__(self, deg):
        """
        Initialize an instance of LineBundle with the specified degree.

        Parameters:
        ----------
        deg : int
            The degree of the line bundle
        """
        self.c0 = 1
        self.c1 = int(deg)
        self.c2 = 0

    def __str__(self):
        """
        String representation of the line bundle class; the standard notation of a line bundle
        of degree d is O(d), since every line bundle is a twist of the structure sheaf.

        Returns:
        -------
        str
            A string representation of the line bundle
        """
        return f'O({self.c1})'
    
    def is_semistable(self, s, q):
        """
        Helper method which helps the program know that line bundles are always stable objects

        Parameters:
        ----------
        s : float
            The parameter controlling the imaginary part of the central charge
        q : float
            The parameter controlling the real part of the central charge

        Returns:
        -------
        bool
            True --- line bundles should always be stable in the geometric chamber
        """
        return True
    


    
class CotangentBundleP2(VectorBundleP2):
    """
    An implementation of the canonical bundle of the complex projective plane. The canonical bundle
    is a rank 2 vector bundle, and is the dual of the tangent bundle; it is generated by the
    differentials of the homogeneous coordinates x, y, and z. Locally, these may be represented as
    d(y/z), d(x/z), and d(x/y) (respectively), which simplify to 

                    d(y/z) = dy/z - y dz/z^2    
                    d(x/z) = dx/z - x dz/z^2

    This gives an embedding of sheaves Ω --> O(-1)^3 whose cokernel is the structure sheaf. 

    Attributes:
    ----------
    c0 : int
        The rank of the cotangent bundle, which is 2
    c1 : int
        The first chern class of the cotangent bundle, which is of the form 2d-3
    c2 : float
        The second chern class of the cotangent bundle, which is of the form d^2 -3d + 3
    """
    def __init__(self, deg=0):
        """
        Initialize an instance of CotangentBundle with the specified degree.
        By default, the degree is 0 so that the class fits into the standard Euler Sequence

        Parameters:
        ----------
        deg : int
            The degree of the cotangent bundle
        """
        self.c0 = 2
        self.c1 = 2*int(deg) - 3
        self.c2 = deg**2 - 3*deg + 3

    def __str__(self):
        """
        String representation of the cotangent bundle class; the standard notation of a cotangent bundle
        of degree d is Ω(d).

        Returns:
        -------
        str
            A string representation of the cotangent bundle
        """
        return f'\u03a9({int((self.c1 + 3) / 2)})'
    
    def chernCharacter(self):
        """
        Method which returns a Chern Characteristic class of the cotangent bundle. The chern
        Character of the cotangent bundle is always of the form (2, 2d - 3, d^2 - 3d + 3).

        Returns:
        -------

        ChernCharacter
            The Chern character of the cotangent bundle
        """
        ch2 = float(self.c1**2 - 2 * self.c2) / 2
        return ChernCharacterP2(self.c0, self.c1, ch2)

