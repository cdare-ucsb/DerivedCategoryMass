from CoherentSheaf import LineBundle
from DerivedCategoryObject import DerivedCategoryObject
from DistinguishedTriangle import DistinguishedTriangle
from ChainComplex import ChainComplex
from ChernCharacter import ChernCharacter
import math
import cmath




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




class SphericalTwist(DerivedCategoryObject):
    
    def __init__(self, line_bundle_1, line_bundle_2, degree=1):
        """
        Initialize an instance of SphericalTwist with the specified line bundles. The spherical twist
        is defined as the cone of the evaluation morphism 

                Hom(i*O(a), i*O(b)) ⊗ i*O(a) ---->  i*O(b) ----> Tw_a O(b)

        where i*O(a) is the pushforward of the line bundle O(a) and Tw_a O(b) is the spherical twist. 
        The spherical twist is represented as a distinguished triangle in the derived category of coherent
        sheaves. 

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
        self.degree = degree

        self.defining_triangle = self._sph_twist_LineBundles(line_bundle_1, line_bundle_2)


    def __str__(self):
        """
        Returns a string representation of the spherical twist by printing the defining triangle

        Returns:
        -------
        str
            A string representation of the spherical twist
        """
        return str(self.defining_triangle.object3)
    

    
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

        homDims = _dimHom_LineBundles(line_bundle_1, line_bundle_2, self.degree)

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
    
    def shift(self, n):

        if not isinstance(n, int):
            raise TypeError("Shift must be an integer")
        
        new_object1 = self.defining_triangle.object1.shift(n)
        new_object2 = self.defining_triangle.object2.shift(n)
        new_object3 = self.defining_triangle.object3.shift(n)

        return DistinguishedTriangle(new_object1, new_object2, new_object3)
    


    
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
            
            ch = self.chernCharacter()
            
            return -1*ch[1] + args[0]*ch[0]

        elif self.catagory == 'P2':
            if len(args) != 2:
                raise ValueError("Central charge of P2 requires two real number parameters")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("P2 objects should have two real number parameters")
            
            ch = self.chernCharacter()
            
            return complex(-1*ch[2] +
                            args[1] * ch[0],
                              ch[1] - args[0] * ch[0])
            

        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")

            alpha = args[0]
            beta = args[1]
            d = args[2]

            ch = self.chernCharacter()
            
            return complex(2*d*alpha * ch[1] - ch[2] - ch[0] + (beta**2 - alpha**2)*d*ch[0], 
                           2*d*ch[1] - 2*d*alpha*beta*ch[0])
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")

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
        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")


        # Write triangle as A -> Tw -> B + B[shift]
        modified_defining_triangle = self.defining_triangle.rotateLeft()
        subobject = modified_defining_triangle.object1.sheaf_vector[0]

        # phase(A)
        left_side_phase = subobject.phase(*args)

        
        quotient_complex = modified_defining_triangle.object3 

        # right_side_phase = 0
        # right_lb_base_phase = quotient_complex.sheaf_vector[0].phase(*args)


        # if len(quotient_complex.shift_vector) == 1:
        #     right_side_shift = quotient_complex.shift_vector[0]
            
        #     right_side_phase = right_lb_base_phase + right_side_shift
        # elif len(quotient_complex.shift_vector) == 2:
        #     # We may assume / Know that for spherical twists, the two line bundles
        #     # will be the same.w

        #     right_side_shift = 0 # should be minimum shift
        #     # Get minimum shift
        #     if quotient_complex.dimension_vector[0] > quotient_complex.dimension_vector[1]:
        #         right_side_shift = quotient_complex.shift_vector[1]
        #     else:
        #         right_side_shift = quotient_complex.shift_vector[0]

        #     right_side_phase = right_lb_base_phase + right_side_shift
        # else:
        #     raise ValueError("For a single spherical twist the Hom object should not be concentrated in more than a single degree")
        

        return left_side_phase <= quotient_complex.get_smallest_phase(*args)
    



    
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
        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")
        


        if self.is_semistable(*args):
            return abs(self.central_charge(*args))
        
        modified_defining_triangle = self.defining_triangle.rotateLeft()
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
                

class SphericalTwistSum(DerivedCategoryObject):


    def __init__(self, line_bundle_pairs_vector, dimension_vector, shift_vector, degree=1):
        """
        Initialize an instance of SphericalTwistSum with the specified line bundles. This class is used
        to represent objects of the form

        Tw_{a_1} O(b_1)^{n_1}[s_1] ⊕ Tw_{a_2} O(b_2)^{n_2}[s_2]}

        where the spherical twists are defined as the cone of the evaluation morphism

                Hom(O(a), O(b)) ⊗ O(a) ---->  O(b) ----> Tw_a O(b)

        This class is primarily used for DoubleSphericalTwist, since the second canonical triangle
        will often be a sum of spherical twists.
        """
        if not all(isinstance(x, tuple) for x in line_bundle_pairs_vector):
            raise TypeError("line_bundle_pairs_vector must be a list of tuples")
        if not all(len(x) == 2 for x in line_bundle_pairs_vector):
            raise ValueError("line_bundle_pairs_vector must be a list of tuples with exactly two elements")
        if not all(isinstance(x[0], LineBundle) and isinstance(x[1], LineBundle) for x in line_bundle_pairs_vector):
            raise TypeError("line_bundle_pairs_vector must be a list of tuples where both objects are line bundles")
        if not all(x[0].catagory == x[1].catagory for x in line_bundle_pairs_vector):
            raise ValueError("Line bundles must be defined on the same variety; not all catagories currently match.")
        if not all(isinstance(x, int) for x in dimension_vector):
            raise TypeError("dimension_vector must be a list of integers")
        if not all( x >= 0 for x in dimension_vector):
            raise ValueError("dimension_vector must be a list of non-negative integers")
        if not all(isinstance(x, int) for x in shift_vector):
            raise TypeError("shift_vector must be a list of integers")
        
        
        if len(line_bundle_pairs_vector) != len(dimension_vector) or len(line_bundle_pairs_vector) != len(shift_vector):
            raise ValueError("line_bundle_pairs_vector, dimension_vector, and shift_vector must have the same length")
        
        if len(line_bundle_pairs_vector) == 0:
            raise ValueError("line_bundle_pairs_vector must have at least one element")
        
        self.line_bundle_pairs_vector = line_bundle_pairs_vector
        self.dimension_vector = dimension_vector
        self.shift_vector = shift_vector
        self.degree = degree
        self.catagory = line_bundle_pairs_vector[0][0].catagory

        self._remove_zeros_from_dimension_vector()
        self._combine_repeats()



    def __str__(self):
        """
        Returns a string representation of the spherical twist by printing the defining triangle

        Returns:
        -------
        str
            A string representation of the spherical twist
        """
        # Initialize lists to hold each formatted component
        first_line_components = []
        second_line_components = []
        
        # Iterate over the zipped input lists
        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            # Format the second line component
            second_line_component = ""
            if s != 0:
                second_line_component = f"Tw_{lb1.degree} O({lb2.degree})[{s}]"
            else:
                second_line_component = f"Tw_{lb1.degree} O({lb2.degree})"

            second_line_components.append(second_line_component)
            
            # Calculate the position to place the '⊕n' above 'O'
            o_position = second_line_component.index(')')

            # Create a string with spaces up to the 'O' position, then add '⊕n'
            first_line_component = ""
            if n != 1:
                first_line_component = ' ' * o_position + f'⊕{n}'
            else:
                first_line_component = ' ' * o_position + '  '
            first_line_components.append(first_line_component)
            
        
        # Join all components with ' ⊕ ' separator
        first_line = ' '.join(first_line_components)
        second_line = ' ⊕ '.join(second_line_components)
        
        # Combine the two lines with a newline character
        if first_line.isspace():
            return second_line
        result = f"{first_line}\n{second_line}"
        
        return result
    
    def chernCharacter(self):
        
        chern_character = ChernCharacter([0, 0, 0])
        
        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            sph_twist_chern = SphericalTwist(lb1, lb2, self.degree).chernCharacter()

            chern_character += int(n * (-1)**s) * sph_twist_chern
        
        return chern_character
    

    
    def shift(self, n):

        new_shift_vector = [x + n for x in self.shift_vector]

        return SphericalTwistSum(self.line_bundle_pairs_vector, self.dimension_vector, new_shift_vector, self.degree)
    
    def central_charge(self, *args):
        """
        Method to compute the central charge of the spherical twist sum. Since all stability conditions are numerical
        for our current implementations, it suffices to compute the central charge of the Chern character

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
        elif self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2, and K3 catagories are implemented")

        central_charge = 0

        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            sph_twist = SphericalTwist(lb1, lb2, self.degree)
            central_charge += n * (-1)**s * sph_twist.central_charge(*args)
        
        return central_charge
    
    def get_lowest_shift_component(self):
        
        min_shift = float('inf')
        min_lb_pair = None
        min_dimension = None

        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            if s < min_shift:
                min_shift = s
                min_lb_pair = (lb1, lb2)
                min_dimension = n
        
        return SphericalTwistSum([min_lb_pair], [min_dimension], [min_shift], self.degree)
            
    def get_highest_shift_component(self):
        
        max_shift = float('-inf')
        max_lb_pair = None
        max_dimension = None

        for (lb1, lb2), n, s in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            if s > max_shift:
                max_shift = s
                max_lb_pair = (lb1, lb2)
                max_dimension = n
        
        return SphericalTwistSum([max_lb_pair], [max_dimension], [max_shift], self.degree)
    
    def is_semistable(self, *args):

        if self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
            
            if len(self.dimension_vector) == 1:
                sph_twist = SphericalTwist(self.line_bundle_pairs_vector[0][0], self.line_bundle_pairs_vector[0][1], self.degree)
                return sph_twist.is_semistable(*args)
            else:
                return False
            
    def _remove_zeros_from_dimension_vector(self):
        """
        Helper function which iterates through the dimension vector, and if a certain Coherent sheaf
        is only included 0 times, we may effectively erase it.
        """
        for i in range(len(self.dimension_vector)):
            if i < len(self.dimension_vector) and self.dimension_vector[i] == 0:
                del self.line_bundle_pairs_vector[i]
                del self.shift_vector[i]
                del self.dimension_vector[i]

    def _combine_repeats(self):
        """
        Helper function to combine repeated sheaves in the complex. This is useful for simplifying
        the complex, as we can combine repeated sheaves into a single sheaf with a larger dimension.
        """

        # Dictionary to hold combined dimensions for each (sheaf, shift) pair
        combined = {}

        # Iterate over the tuples
        for (lb1,lb2), dim, shift in zip(self.line_bundle_pairs_vector, self.dimension_vector, self.shift_vector):
            key = (lb1, lb2, shift)
            if key in combined:
                combined[key] += dim
            else:
                combined[key] = dim

        # Unpack the combined dictionary back into the class variables
        self.line_bundle_pairs_vector = []
        self.dimension_vector = []
        self.shift_vector = []

        for (lb1, lb2, shift), dim in combined.items():
            self.line_bundle_pairs_vector.append((lb1, lb2))
            self.dimension_vector.append(dim)
            self.shift_vector.append(shift)


    

        






class DoubleSphericalTwist(DerivedCategoryObject):

    def __init__(self, line_bundle_1, line_bundle_2, line_bundle_3, degree=1):
        """
        Initialize an instance of DoubleSphericalTwist with the specified line bundles. The spherical twist
        is defined as the cone of the evaluation morphism 

                Hom(O(a), Tw_b O(c)) ⊗ O(a) ---->  Tw_b O(c) ----> Tw_a Tw_b O(c)

        where O(a), O(b), and O(c) denote line bundles. 
        The double spherical twist is represented as a distinguished triangle in the derived category of coherent
        sheaves. 

        Several helper methods are used to compute the dimensions of the Hom spaces between the pushforwards
        of the line bundles, and then to construct the distinguished triangle.

        Parameters:
        ----------
        line_bundle_1 : LineBundle
            The last line bundle twisted around; i.e. O(a) where we are computing Tw_a Tw_b O(c)
        line_bundle_2 : LineBundle
            The first line bundle twisted around; i.e. O(b) where we are computing Tw_a Tw_b O(c)
        line_bundle_3 : LineBundle
            The line bundle we are applying the spherical twist to; i.e. O(c) where we are computing Tw_a Tw_b O(c)
        degree : int
            Optional argument for degree of K3 surface

        """
        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundleP1.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundleP1.")
        if not isinstance(line_bundle_3, LineBundle):
            raise TypeError("line_bundle_3 must be an instance of LineBundleP1.")
        
        if line_bundle_1.catagory != line_bundle_2.catagory or line_bundle_1.catagory != line_bundle_3.catagory:
            raise ValueError("Line bundles must be defined on the same variety")

        self.line_bundle_1 = line_bundle_1
        self.line_bundle_2 = line_bundle_2
        self.line_bundle_3 = line_bundle_3
        self.catagory = line_bundle_1.catagory
        self.degree = degree

        self.defining_triangle = self._sph_twist_DoubleLineBundles(line_bundle_1, line_bundle_2, line_bundle_3)

        if self.catagory != 'K3':
            raise ValueError("Double spherical twists are currently only implemented for K3 surfaces")
        

    def _sph_twist_DoubleLineBundles(self, line_bundle_1, line_bundle_2, line_bundle_3):

        if not isinstance(line_bundle_1, LineBundle):
            raise TypeError("line_bundle_1 must be an instance of LineBundle.")
        if not isinstance(line_bundle_2, LineBundle):
            raise TypeError("line_bundle_2 must be an instance of LineBundle.")
        if not isinstance(line_bundle_3, LineBundle):
            raise TypeError("line_bundle_3 must be an instance of LineBundle.")
            

        homDims = _dimHom_Line_and_SingleTwist(line_bundle_1, line_bundle_2, line_bundle_3, self.degree)  

        bundle_vector = []
        dimension_vector = [] 
        shift_vector = []

        # create the necessary lists for the ChainComplex constructor 
        for i in range(len(homDims)):
            if homDims[i] == 0:
                continue
            dimension_vector.append(homDims[i])
            # The homDims tuple represents the degrees (-1, 0, 1, 2, 3). This is represented by cohomological
            # shifts [1], [0], [-1], [-2], and [-3] respectively.
            shift_vector.append(1 - i)
            bundle_vector.append(LineBundle(line_bundle_1.degree, self.catagory))

        object1 = ChainComplex(sheaf_vector=bundle_vector, shift_vector=shift_vector, dimension_vector=dimension_vector)

        object2 = SphericalTwistSum([(line_bundle_2, line_bundle_3)], dimension_vector=[1], shift_vector=[0], degree=self.degree)
        object3 = DerivedCategoryObject(string=f"Tw_{line_bundle_1.degree} Tw_{line_bundle_2.degree} O({line_bundle_3.degree})", catagory=self.catagory)

        return DistinguishedTriangle(object1, object2, object3)
    
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

        if self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")

            alpha = args[0]
            beta = args[1]
            d = args[2]

            ch = self.chernCharacter()
            
            return complex(2*d*alpha * ch[1] - ch[2] - ch[0] + (beta**2 - alpha**2)*d*ch[0], 
                           2*d*ch[1] - 2*d*alpha*beta*ch[0])
        else:
            raise NotImplementedError("Only P1, P2 and K3 catagories are implemented")

        
    
    def __str__(self):
        """
        Returns a string representation of the spherical twist by printing the defining triangle

        Returns:
        -------
        str
            A string representation of the spherical twist
        """
        return str(self.defining_triangle.object3)
    
    def secondary_canonical_triangle(self):
        """
        Method to compute the secondary canonical triangle of the spherical twist. The secondary canonical triangle
        is obtained by applying Tw_a to the defining triangle of Tw_b O(c); specifically, this gives

        Tw_a (Hom(O(b), O(c)) ⊗ O(b)) ----> Tw_a O(c) ----> Tw_a Tw_b O(c)

        Returns:
        -------
        DistinguishedTriangle
            The secondary canonical triangle of the spherical twist
        """
        # Get the object Tw_a O(c)
        object2 = SphericalTwistSum([(self.line_bundle_1, self.line_bundle_3)], dimension_vector=[1], shift_vector=[0], degree=self.degree)
        

        object3 = self.defining_triangle.object3

        single_twist_original = SphericalTwist(self.line_bundle_2, self.line_bundle_3, self.degree)
        hom_complex_single_twist = single_twist_original.defining_triangle.object1

        dimension_vec = []
        shift_vec = []
        line_bundle_pair_vec = []

        for i in range(len(hom_complex_single_twist.sheaf_vector)):
            if hom_complex_single_twist.dimension_vector[i] == 0:
                continue
            dimension_vec.append(hom_complex_single_twist.dimension_vector[i])
            shift_vec.append(hom_complex_single_twist.shift_vector[i])
            line_bundle_pair_vec.append((self.line_bundle_1, self.line_bundle_2))

        object1 = SphericalTwistSum(line_bundle_pair_vec, dimension_vec, shift_vec, self.degree)

        return DistinguishedTriangle(object1, object2, object3)
            
            

        

        

    
    def is_semistable(self, *args):
        
        
        if self.catagory == 'K3':

            if len(args) != 3:
                raise ValueError("Central charge of K3 requires three real number parameters: alpha, beta, and the degree")
            if not all(isinstance(x, (float, int)) for x in args):
                raise TypeError("K3 central charges should have three real number parameters: alpha, beta, and the degree")
            if not isinstance(args[2], int):
                raise TypeError("The degree of the K3 surface must be an integer")
        else:
            raise NotImplementedError("Only P1, P2 and K3 catagories are implemented")

        
        ###########
        # First check the defining triangle
        ###########

        # Write triangle as Tw_b O(c) -> Tw_a Tw b O(c) -> B + B[shift] + ...
        subobject = SphericalTwist(self.line_bundle_2, self.line_bundle_3, self.degree)

        modified_defining_triangle = self.defining_triangle.rotateLeft()
        quotient_complex = modified_defining_triangle.object3 

        print("First checking defining triangle:")
        print("-------------------\n")
        print(modified_defining_triangle)

        if subobject.is_semistable(*args):
            print(f"{subobject} is semistable")
            # get the phase of the single twist
            left_side_phase = cmath.phase(subobject.central_charge(*args)) / math.pi
            right_side_min_phase = quotient_complex.get_smallest_phase(*args)
            right_side_max_phase = quotient_complex.get_largest_phase(*args)

            if left_side_phase <= right_side_min_phase:
                return True
            elif left_side_phase > right_side_max_phase:
                return False
        else:
            print(f"{subobject} is not semistable")

        # nothing can be gathered from first triangle

        ###########
        # next check secondary canonical triangle
        ###########

        secondary_canonical_triangle = self.secondary_canonical_triangle().rotateLeft()
        print("\n Second checking secondary canonical triangle:")
        print("-------------------\n")
        print(secondary_canonical_triangle)

        
            









###################################################################
#                  Static Helper Methods                          #
###################################################################


def _dimHom_LineBundles(line_bundle_1, line_bundle_2, degree_K3 = 1):

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
    elif catagory == 'K3':
        return _dimHom_LineBundlesK3(line_bundle_1, line_bundle_2, degree_K3)

    else:
        raise NotImplementedError("Only P1, P2 and K3 catagories are implemented")


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
        
def _dimHom_LineBundlesK3(line_bundle_1, line_bundle_2, degree_K3):

    degree_dif = line_bundle_2.degree - line_bundle_1.degree

    if degree_dif == 0:
        return (1, 0, 1)
    elif degree_dif > 0:
        return (degree_K3 * degree_dif**2 + 2, 0, 0)
    else:
        return (0, 0, degree_K3 * degree_dif**2 + 2)

        

def _dimHom_Line_and_SingleTwist(line_bundle_1, line_bundle_2, line_bundle_3, degree_K3):
    if not isinstance(line_bundle_1, LineBundle):
        raise TypeError("line_bundle_1 must be an instance of LineBundle.")
    if not isinstance(line_bundle_2, LineBundle):
        raise TypeError("line_bundle_2 must be an instance of LineBundle.")
    if not isinstance(line_bundle_3, LineBundle):
        raise TypeError("line_bundle_3 must be an instance of LineBundle.")
    

    if line_bundle_1.catagory != line_bundle_2.catagory or line_bundle_1.catagory != line_bundle_3.catagory:
        raise ValueError("Line bundles must be defined on the same variety")
    

    if line_bundle_1.catagory == 'K3':
        return _dimHom_Line_and_SingleTwistK3(line_bundle_1, line_bundle_2, line_bundle_3, degree_K3)
    else:
        raise NotImplementedError("Double spherical twist is currently only implemented for K3 surfaces.")




# TODO: Implement the Hom spaces for the case of a single twist and a line bundle
def _dimHom_Line_and_SingleTwistK3(line_bundle_1, line_bundle_2, line_bundle_3, degree_K3):
    """
    Helper function which computes the dimensions of the derived hom spaces between a line bundle 
    and a spherical twist of line bundles; that is, this function computes the dimensions of 
    
               RHom(O(a), Tw_{O(b)} O(c))
    
    which is necessarily concentrated in cohomological degrees -1 through 3. As a consequence, this 
    function returns a tuple of five integers (hom-1, hom0, hom1, hom2, hom3) corresponding to the 
    dimensions of the higher-ext groups. 


    Parameters:
    ----------
    line_bundle_1 : LineBundle
        The last line bundle to twist around: i.e. O(a) where we are computing Tw_a Tw_b O(c)
    line_bundle_2 : LineBundle
        The first line bundle to twist around: i.e. O(b) where we are computing Tw_a Tw_b O(c)
    line_bundle_3 : LineBundle
        The line bundle which the spherical twists are being applied to: i.e. O(c) where we are computing Tw_a Tw_b O(c)

    """


    if not isinstance(degree_K3, int):
        raise TypeError("The degree of the K3 surface must be an integer")
    if not degree_K3 > 0:
        raise ValueError("The degree of the K3 surface must be a positive integer")
    if not isinstance(line_bundle_1, LineBundle):
        raise TypeError("line_bundle_1 must be an instance of LineBundleP2.")
    if not isinstance(line_bundle_2, LineBundle):
        raise TypeError("line_bundle_2 must be an instance of LineBundleP2.")
    if not isinstance(line_bundle_3, LineBundle):
        raise TypeError("line_bundle_3 must be an instance of LineBundleP2.")

    dif_12 = degree_K3 * (line_bundle_2.degree - line_bundle_1.degree)**2 + 1
    dif_23 = degree_K3 * (line_bundle_3.degree - line_bundle_2.degree)**2 + 1
    dif_13 = degree_K3 * (line_bundle_3.degree - line_bundle_1.degree)**2 + 1

    if line_bundle_2.degree < line_bundle_3.degree:

        if line_bundle_1.degree < line_bundle_2.degree:
            if dif_12 * dif_23 >= dif_13:
                return (dif_12 * dif_23 - dif_13, 0, 0, 0, 0)
            else:
                # There is an edge case when (b-c)(b-a) = -1; this reduces to Tw_{n+1} Tw_n O(n-1)
                # or Tw_{n-1} Tw_n O(n+1). In this case, the Hom space is concentrated in degree 0
                return (0, dif_13 - dif_12*dif_23, 0, 0, 0)
        elif line_bundle_1.degree == line_bundle_2.degree:
            return (0, 0, dif_23, 0, 0)
        elif line_bundle_1.degree > line_bundle_2.degree and line_bundle_1.degree < line_bundle_3.degree:
            return (0, dif_13, dif_12*dif_23, 0, 0)
        elif line_bundle_1.degree == line_bundle_3.degree:
            return (0, 1, dif_23 - 1, 0, 0)
        else:
            # line_bundle_1.deg > line_bundle_3.deg
            if dif_12 * dif_23 >= dif_13:
                return (0, 0, dif_23*dif_12-dif_13, 0, 0)
            else:
                # There is an edge case when (b-c)(b-a) = -1; this reduces to Tw_{n+1} Tw_n O(n-1)
                # or Tw_{n-1} Tw_n O(n+1). In this case, the Hom space is concentrated in degree 0
                # However, this theoretically should not occur here since we are considering the case
                # where a > c > b
                return (0, 0, 0, dif_13 - dif_12*dif_23, 0)
            
    elif line_bundle_2.degree == line_bundle_3.degree:

        # If b = c, then Tw_b O(c) = O(c)[-1]. Thus,
        #          Hom(O(a), Tw_b O(c)) = Hom(O(a), O(c)[-1]) = Hom(O(a), O(c))[-1] 

        (a,b,c) = _dimHom_LineBundlesK3(line_bundle_1, line_bundle_2, degree_K3)
        return (0, 0, a, b, c)
    
    else:
        # line_bundle_2.deg > line_bundle_3.deg
        if line_bundle_1.degree > line_bundle_2.degree:
            return (0, 0, 0, dif_13, dif_23*dif_13)
        elif line_bundle_1.degree == line_bundle_2.degree:
            return (0, 0, 0, 0, dif_23)
        elif line_bundle_1.degree < line_bundle_2.degree and line_bundle_1.degree > line_bundle_3.degree:
            if dif_12 * dif_23 >= dif_13:
                return (0, 0, dif_23*dif_12 - dif_13, 0, 0)
            else:
                return (0, 0, 0, dif_13 - dif_12*dif_23, 0)
        elif line_bundle_1.degree == line_bundle_3.degree:
            return (0, 1, dif_12*dif_23 - 1, 0, 0)
        else:
            # line_bundle_1.deg < line_bundle_3.deg
            return (0, dif_13, dif_12*dif_23, 0, 0)
        



            

        


###################################################################
#                        Main                                     #
###################################################################


if __name__ == "__main__":
    # lb1 = LineBundle(-1, catagory='P2')
    # lb2 = LineBundle(0, catagory='P2')
    # lb3 = LineBundle(2, catagory='P2')

    # sph = SphericalTwist(line_bundle_1=lb1, line_bundle_2=lb3)

    # print(sph)


    # lb4 = LineBundle(-5, catagory='P1')
    # lb5 = LineBundle(2, catagory='P1')

    # sph2 = SphericalTwist(line_bundle_1=lb4, line_bundle_2=lb5, degree=1)

    # print(sph2)


    lb6 = LineBundle(4, catagory='K3')
    lb7 = LineBundle(10, catagory='K3')
    lb8 = LineBundle(7, catagory='K3')
    lb9 = LineBundle(5, catagory='K3')

    # sph_sum = SphericalTwistSum([(lb6, lb7), (lb8, lb9), (lb6, lb9), (lb6, lb7)], [1, 5, 0, 10], [-1, 2, -5, -1], degree=1)
    # print(sph_sum)
    # print(sph_sum.chernCharacter())

    


    # sph3 = SphericalTwist(lb7, lb8, degree=2)




    # print(sph3.shift(3))
    # print("Is semistable: ", sph3.is_semistable(1, 2, 1))

    

    sph4 = DoubleSphericalTwist(lb7, lb8, lb6, degree=1)

    sph4.is_semistable(1, 2, 1)

    

    