from .CoherentSheaf import LineBundle
from .DerivedCategoryObject import DerivedCategoryObject
from .DistinguishedTriangle import DistinguishedTriangle
from .ChainComplex import ChainComplex
import math




###############################################################################
#                                                                             #
#                           Spherical Twist                                   #
# ----------------------------------------------------------------------------#
#  This object is used to represent the composition of spherical twists in    #
#  in the derived category of local projective space  consisting of sheaves   #
#  supported on the zero divisor (the original projective space)              #
#  In particular, the only compositions of twists we consier are twists       #
#  around (pushforwards) of line bundles. Theoretically, since the derived    #
#  category of coherent sheaves on projective space is constructible, the     #
#  braid relations between spherical twists allow any object obtained as a    #
#  series of spherical twists to be represented as a sequence of spherical    #
#  twists specifically around line bundles.                                   #
#                                                                             #
#  In order to determine possible Harder-Narasimhan filtrations of the        #
#  spherical twists, we need to iteratively keep track of previous Harder-    #
#  Narasimhan filtrations. This is done by expanding the leaves of a binary   #
#  tree, where each node is a spherical twist, and the children of each node  #
#  are the spherical twists obtained by applying a successive twist.          #
#                                                                             #   
###############################################################################      




class SphericalTwist:
    
    def __init__(self, line_bundle_1, line_bundle_2):
        """
        Initialize an instance of SphericalTwist with the specified line bundles. The spherical twist
        is defined as the cone of the evaluation morphism 

                Hom(i*O(a), i*O(b)) ⊗ i*O(a) ---->  i*O(b) ----> Tw_a O(b)

        where i*O(a) is the pushforward of the line bundle O(a) and Tw_a O(b) is the spherical twist. 
        The spherical twist is represented as a distinguished triangle in the derived category of coherent
        sheaves on local P^1. 

        Several helper methods are used to compute the dimensions of the Hom spaces between the pushforwards
        of the line bundles, and then to construct the distinguished triangle.

        """
        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundleP1.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundleP1.")
        
        if line_bundle_1.catagory != line_bundle_2.catagory:
            raise ValueError("Line bundles must be defined on the same variety")

        self.line_bundle_1 = line_bundle_1
        self.line_bundle_2 = line_bundle_2
        self.catagory = line_bundle_1.catagory

        self.defining_triangle = self._sph_twist_LineBundles(line_bundle_1, line_bundle_2)


    def __str__(self):
        """
        Returns a string representation of the spherical twist by printing the defining triangle

        Returns:
        -------
        str
            A string representation of the spherical twist
        """
        return str(self.defining_triangle)
    

    
    def _sph_twist_LineBundles(self, line_bundle_1, line_bundle_2):
        """
        Helper function which uses the __dimHom_LineBundlesP2 method to compute the defining triangle for a 
        single spherical twist of a line bundle around another line bundle. The twist is given by

        Tw_lb1 lb2 = Cone(  Hom(lb_1, lb_2) ⊗ lb_1 ---->  lb_2 )


        Parameters:
        ----------
        line_bundle_1 : LineBundle
            The first line bundle in the Hom space
        line_bundle_2 : LineBundle
            The second line bundle in the Hom space

        Returns:
        -------
        DistinguishedTriangle
            The distinguished triangle representing the spherical twist

        Raises:
        -------
        TypeError
            If line_bundle_1 is not an instance of LineBundle
            If line_bundle_2 is not an instance of LineBundle
        """

        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundleP1.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundleP1.")

        homDims = _dimHom_LineBundles(line_bundle_1, line_bundle_2)

        bundle_vector = []
        dimension_vector = [] 
        shift_vector = []

        # create the necessary lists for the ChainComplex constructor 
        for i in range(len(homDims)):
            if homDims[i] == 0:
                continue
            dimension_vector.append(homDims[i])
            shift_vector.append(-1*i)
            bundle_vector.append(LineBundle(line_bundle_1.degree, self.catagory))

        object1 = ChainComplex(sheaf_vector=bundle_vector, shift_vector=shift_vector, dimension_vector=dimension_vector)
        object2 = ChainComplex(sheaf_vector=[LineBundle(line_bundle_2.degree, self.catagory)], shift_vector=[0], dimension_vector=[1])
        object3 = DerivedCategoryObject(string=f"Tw_{line_bundle_1.degree} O({line_bundle_2.degree})", catagory=self.catagory)

        return DistinguishedTriangle(object1, object2, object3)
    
    def chernCharacter(self):
        """
        Method to compute the Chern Character of the spherical twist. The Chern Character of the
        spherical twist is the Chern Character of the third object in the distinguished triangle.

        Returns:
        -------
        ChernCharacter
            The Chern Character of the spherical twist
        """

        return self.defining_triangle.object3.chernCharacter()
    


    
    def central_charge(self, *args):
        """
        Method to compute the central charge of the spherical twist. The central charge of the spherical
        twist is the central charge of the third object in the distinguished triangle.

        Parameters:
        ----------
        w : complex
            The complex parameter for the stability condition on local P1

        Returns:
        -------
        complex
            The central charge of the spherical twist as a complex number

        Raises:
        -------
        TypeError
            If w is not a complex number
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
        else:
            raise NotImplementedError("Only P1 and P2 catagories are implemented")

        return self.defining_triangle.object3.central_charge(*args)
    


    
    def is_semistable(self, *args):
        """
        Method to check if the spherical twist is stable. A spherical twist is stable if the
        Chern Character of the third object in the distinguished triangle is (0,0,0). This is
        equivalent to the condition that the Chern Character of the first object is equal to
        the Chern Character of the second object.


        Parameters:
        -------
        w : complex
            The complex parameter for the stability condition on local P1


        Returns:
        -------
        bool
            True if the spherical twist is stable, False otherwise
        """

        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters. Currently {} parameters given: {}".format(len(args), args))
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
        else:
            raise NotImplementedError("Only P1 and P2 catagories are implemented")


        # Write triangle as A -> Tw -> B + B[shift]
        modified_defining_triangle = self.defining_triangle.shiftLeft()
        subobject = modified_defining_triangle.object1.sheaf_vector[0]

        # phase(A)
        left_side_phase = subobject.phase(*args)

        
        quotient_complex = modified_defining_triangle.object3 

        right_side_phase = 0
        right_lb_base_phase = quotient_complex.sheaf_vector[0].phase(*args)


        if len(quotient_complex.shift_vector) == 1:
            right_side_shift = quotient_complex.shift_vector[0]
            
            right_side_phase = right_lb_base_phase + right_side_shift
        elif len(quotient_complex.shift_vector) == 2:
            # We may assume / Know that for spherical twists, the two line bundles
            # will be the same.w

            right_side_shift = 0 # should be minimum shift
            # Get minimum shift
            if quotient_complex.dimension_vector[0] > quotient_complex.dimension_vector[1]:
                right_side_shift = quotient_complex.shift_vector[1]
            else:
                right_side_shift = quotient_complex.shift_vector[0]

            right_side_phase = right_lb_base_phase + right_side_shift
        else:
            raise ValueError("For a single spherical twist the Hom object should not be concentrated in more than a single degree")
        

        return left_side_phase <= right_side_phase
    



    
    def mass(self, *args):


        if self.catagory == 'P1':
            if len(args) != 1:
                raise ValueError("Central charge of P1 requires single complex number parameter")
            if not isinstance(args[0], complex):
                raise TypeError("P1 objects should have a single complex parameter")
        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
        else:
            raise NotImplementedError("Only P1 and P2 catagories are implemented")
        


        if self.is_semistable(*args):
            return abs(self.central_charge(*args))
        
        modified_defining_triangle = self.defining_triangle.shiftLeft()
        subobject = modified_defining_triangle.object1.sheaf_vector[0]

        quotient_complex = modified_defining_triangle.object3

        mass = 0
        
        if len(quotient_complex.dimension_vector) == 1:

            # Write triangle as O(a) -> Tw -> O(b)[shift]
            
            mass = abs(subobject.central_charge(*args))
            mass += quotient_complex.dimension_vector[0] * abs(quotient_complex.sheaf_vector[0].central_charge(*args))
            
            return mass
            
        else:
            
            if len(quotient_complex.dimension_vector) != 2:
                raise ValueError("The Hom object is not concentrated in 1 or 2 degrees")
            
            phase0 = quotient_complex.sheaf_vector[0].phase(*args) + quotient_complex.shift_vector[0]
            phase1 = quotient_complex.sheaf_vector[1].phase(*args) + quotient_complex.shift_vector[1]

            largest_phase = max(phase0, phase1)

            # CASE 1: phi(subobj) > largest phase(quotient)
            if subobject.phase(*args) > largest_phase:
                mass = abs(subobject.central_charge(*args))
                # By BDL we have that the mass is the sum of masses of two objects
                for i in range(len(quotient_complex.sheaf_vector)):
                    dim = quotient_complex.dimension_vector[i]
                    mass += dim * abs(quotient_complex.sheaf_vector[i].central_charge(*args))
                
                return mass
            # CASE 2: smallest phase(Quotient) < phi(subobj) < largest phase(quotient)
            else:
                if largest_phase == phase0:
                    # smallest phase object first
                    mass = abs(quotient_complex.sheaf_vector[1].central_charge(*args))
                    # isolate single element of larger shift in the quotient object
                    new_complex = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[0]],
                                                shift_vector=[quotient_complex.shift_vector[0]],
                                                dimension_vector=[quotient_complex.dimension_vector[0]])
                    mass += abs(subobject.central_charge(*args) + new_complex.central_charge(*args))

                    return mass

                else:
                    # smallest phase object first
                    mass = abs(quotient_complex.sheaf_vector[0].central_charge(*args))
                    # isolate single element of larger shift in the quotient object
                    new_complex = ChainComplex(sheaf_vector=[quotient_complex.sheaf_vector[1]],
                                                shift_vector=[quotient_complex.shift_vector[1]],
                                                dimension_vector=[quotient_complex.dimension_vector[1]])
                    mass += abs(subobject.central_charge(*args) + new_complex.central_charge(*args))

                    return mass












###################################################################
#                  Static Helper Methods                          #
###################################################################


def _dimHom_LineBundles(line_bundle_1, line_bundle_2):

    if not isinstance(line_bundle_1, LineBundle):
        raise TypeError("line_bundle_1 must be an instance of LineBundle.")
    if not isinstance(line_bundle_2, LineBundle):
        raise TypeError("line_bundle_2 must be an instance of LineBundle.")
    
    catagory = line_bundle_1.catagory

    if line_bundle_1.catagory != line_bundle_2.catagory:
        raise ValueError("Line bundles must be defined on the same variety")
    
    if catagory == 'P1':
        return _dimHom_LineBundlesP1(line_bundle_1, line_bundle_2)
    elif catagory == 'P2':
        return _dimHom_LineBundlesP2(line_bundle_1, line_bundle_2)
    else:
        raise NotImplementedError("Only P1 and P2 catagories are implemented")


def _dimHom_LineBundlesP1(line_bundle_1, line_bundle_2):
        """
        Helper method which computes the dimension of the hom spaces between the pushforwards of the
        line bundles O(a) and O(b). The dimensions of the pushforwards are computed using the triangle

        i^* i_* E -> E -> E x O(2)[2]

        and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        the Hom spaces between line bundles on P1, the computation reduces to a case-by-case
        combinatorial problem. Since the homological index of the hom-space on P1 is bounded between
        0 and 1, the hom-space for local P1 is concentrated between degrees 0 and 2. Thus, we return
        a tuple of the form (a,b,c)

        Parameters:
        ----------
        line_bundle_1 : LineBundleP1
            The first line bundle in the Hom space
        line_bundle_2 : LineBundleP1
            The second line bundle in the Hom space

        Returns:
        -------
        tuple
            A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles

        Raises:
        -------
        TypeError
            If line_bundle_1 is not an instance of LineBundleP1
            If line_bundle_2 is not an instance of LineBundleP1
        """



        degree_dif = line_bundle_2.degree - line_bundle_1.degree

        if degree_dif == 0:
            return (1, 0, 1)
        if degree_dif >= 2:
            return (degree_dif + 1, degree_dif - 1, 0)
        elif degree_dif == 1:
            return (2, 0, 0)
        elif degree_dif == -1:
            return (0, 0, 2)
        else:
            return (0, -1*degree_dif - 1, -1*degree_dif + 1)



def _dimHom_LineBundlesP2(line_bundle_1, line_bundle_2):
        """
        Helper method which computes the dimension of the hom spaces between the pushforwards of the
        line bundles O(a) and O(b). The dimensions of the pushforwards are computed using the triangle

        i^* i_* E -> E -> E x O(3)[2]

        and applying Hom(-, O(b)) to obtain a long-exact sequence. Using the standard dimensions of
        the Hom spaces between line bundles on P^2, the computation reduces to a case-by-case
        combinatorial problem. Since the homological index of the hom-space on P^2 is bounded between
        0 and 2, the hom-space for local P2 is concentrated between degrees 0 and 3. Thus, we return
        a tuple of the form (a,b,c,d)

        Parameters:
        ----------
        line_bundle_1 : LineBundle
            The first line bundle in the Hom space
        line_bundle_2 : LineBundle
            The second line bundle in the Hom space

        Returns:
        -------
        tuple
            A tuple of the dimensions of the Hom spaces between the pushforwards of the line bundles

        Raises:
        -------
        TypeError
            If line_bundle_1 is not an instance of LineBundle
            If line_bundle_2 is not an instance of LineBundle
        """


        degree_dif = line_bundle_2.degree - line_bundle_1.degree

        if degree_dif == 0:
            return (1, 0, 0, 1)
        elif degree_dif > -3 and degree_dif < 0:
            rank3 = math.comb(line_bundle_1.degree - line_bundle_2.degree + 2, 2)
            return (0, 0, 0, rank3)
        elif degree_dif > 0 and degree_dif < 3:
            rank0 = math.comb(degree_dif + 2, 2)
            return (rank0, 0, 0, 0)
        elif degree_dif >= 3:
            rank0 = math.comb(degree_dif + 2, 2)
            rank1 = math.comb(degree_dif - 1, 2)
            return (rank0, rank1, 0, 0)
        else:
            rank2 = math.comb(line_bundle_1.degree - line_bundle_2.degree -1, 2)
            rank3 = math.comb(line_bundle_1.degree - line_bundle_2.degree + 2, 2)    
            return (0, 0, rank2, rank3)
        


# TODO: Implement the Hom spaces for the case of a single twist and a line bundle
def _dimHom_Line_and_SingleTwist(line_bundle_1, line_bundle_2, line_bundle_3):

    if not isinstance(line_bundle_1, LineBundleP2):
        raise TypeError("line_bundle_1 must be an instance of LineBundleP2.")
    if not isinstance(line_bundle_2, LineBundleP2):
        raise TypeError("line_bundle_2 must be an instance of LineBundleP2.")
    if not isinstance(line_bundle_3, LineBundleP2):
        raise TypeError("line_bundle_3 must be an instance of LineBundleP2.")

    degree_dif2 = line_bundle_3.degree - line_bundle_2.degree

    if degree_dif2 > 0 and degree_dif2 < 3:

        multiplier = math.comb(degree_dif2 + 2, 2)
        
        
        if line_bundle_2.c1 - line_bundle_1.c1 > 0:

            ## Requires better understanding of 
            #  Hom(O(-3), O)^{+ 6} --> Hom(O(-3), O(2)) ---> Hom(O(-3), Tw) 

            pass

        


###################################################################
#                        Main                                     #
###################################################################


if __name__ == "__main__":
    lb1 = LineBundle(-1, catagory='P2')
    lb2 = LineBundle(0, catagory='P2')
    lb3 = LineBundle(2, catagory='P2')

    sph = SphericalTwist(line_bundle_1=lb1, line_bundle_2=lb3)

    print(sph)


    lb4 = LineBundle(-5, catagory='P1')
    lb5 = LineBundle(2, catagory='P1')

    sph2 = SphericalTwist(line_bundle_1=lb4, line_bundle_2=lb5)

    print(sph2)