import math


class ChainComplex():
    def __init__(self):
        self.complex = []

    def __init__(self, complex_vector, shift):

        if not all(isinstance(obj, CoherentSheaf) for obj in complex_vector):
            raise TypeError("All elements of complex_vector must be instances of CoherentSheaf.")
        self.complex = complex_vector
        self.shift = shift

    def __str__(self):
        """
        Prints a chain complex with perfectly aligned indices and chain elements.
        
        Args:
            vector (list): A list of integers representing the chain complex.
            shift (int): An integer representing the starting shift in the index.
        """
        # Prepare the indices as strings
        indices = [str(self.shift + i) for i in range(len(self.complex))]
        
        # Calculate the total width of each element block (element + arrow length)
        element_widths = [len(f"O({v})") for v in self.complex]
        
        # Build the index line with precise padding
        index_line = indices[0]  # Start with the first index
        for i in range(1, len(indices)):
            padding = " " * (element_widths[i-1] + 3 - len(indices[i-1]))  # Adjust padding based on the previous element
            index_line += padding + indices[i]
        
        # Build the chain complex line
        complex_line = " ---> ".join(f"{v}" for v in self.complex)
        
        # Print the result
        return index_line + "\n" + complex_line
    

    def chernClass(self):
        '''
        Compute alternating sum of Chern classes of the complex
        '''
        cherns = [obj.chernClass() for obj in self.complex]

        ch0 = 0
        ch1 = 0
        ch2 = 0

        for i in range(len(cherns)):
            ch0 += (-1)**(i + self.shift) * cherns[i].ch0
            ch1 += (-1)**(i + self.shift) * cherns[i].ch1
            ch2 += (-1)**(i + self.shift) * cherns[i].ch2

        return ChernClass(ch0, ch1, ch2)
    
    def shiftComplex(self, shift):
        return ChainComplex(self.complex, self.shift + shift)






class DistinguishedTriangle():
    def __init__(self, chain_complex1, chain_complex2, chain_complex3):
        self.complex1 = chain_complex1
        self.complex2 = chain_complex2
        self.complex3 = chain_complex3

    def __str__(self):
        ret_str = 'A ---> B ---> C'
        ret_str += '\n\nA:\n' + str(self.complex1)
        ret_str += '\n\nB:\n' + str(self.complex2)
        ret_str += '\n\nC:\n' + str(self.complex3)
        return ret_str
    
    def shiftLeft(self):
        return DistinguishedTriangle(self.complex3.shiftComplex(-1), self.complex1, self.complex2)
    def shiftRight(self):
        return DistinguishedTriangle(self.complex2, self.complex3, self.complex1.shiftComplex(1))
    

    



class ChernClass():
    def __init__(self, ch0, ch1, ch2): 
        '''
        ch0 = rank
        ch1 = degree
        ch2 = c2 - c_1^2 / 2
        '''
        self.ch0 = ch0
        self.ch1 = ch1
        self.ch2 = ch2

    def __str__(self):
        return f'<{self.ch0}, {self.ch1}, {self.ch2}>'
    

class CoherentSheaf():
    def __init__(self, rank, deg):
        self.rank = int(rank)
        self.deg = int(deg)
        self.c2 = 0
    def __init__(self, rank, deg, c2):
        self.rank = int(rank)
        self.deg = int(deg)
        self.c2 = float(c2)
    
    def chernClass(self):
        return ChernClass(self.rank, self.deg, self.c2)

    def __str__(self):
        return f'CoherentSheaf of rank {self.rank} and degree {self.deg}'
    





class VectorBundle(CoherentSheaf):

    def __init__(self, rank, deg, c2):
        self.rank = int(rank)
        self.deg = int(deg)
        self.c2 = c2

    def __str__(self):
        return f'VectorBundle of rank {self.rank} and degree {self.deg}'
    
    def chernClass(self):
        return ChernClass(self.rank, self.deg, self.c2)
    
    def isLineBundleSum(self):
        return self.c2 == self.deg**2 / 2



class LineBundle(VectorBundle):
    def __init__(self, deg):
        self.c0 = 1
        self.c1 = int(deg)
        self.c2 = float(deg**2 / 2)

    def __str__(self):
        return f'O({self.c1})'
    
    def chernClass(self):
        return ChernClass(self.c0, self.c1, self.c2)


    
class CotangentBundle(VectorBundle):
    def __init__(self, deg):
        self.c0 = 2
        self.c1 = 2*int(deg) - 3
        self.c2 = float(int(deg)**2 - 3 * int(deg) + float(3/2))

    def __str__(self):
        return f'\u03a9({int((self.c1 + 3) / 2)})'
    
    def chernClass(self):
        return ChernClass(self.c0, self.c1, self.c2)




class P2():

    def __init__(self):
        self.hyperplane = LineBundle(1)

    def __str__(self):
        return 'Projective Plane'
    
    def dimHom(self, line_bundle_1, line_bundle_2):
        '''
        use standard (Hom0, Hom1, Hom2, shift)
         '''
        
        degree_diff = line_bundle_2.deg - line_bundle_1.deg
        
        if degree_diff >= 0:
            hom0_dim = math.comb(degree_diff + 2, 2)
            return (hom0_dim, 0, 0, 0)
        elif degree_diff < 0 and degree_diff > -3:
            return (0, 0, 0, 0)
        else:
            hom2_dim = math.comb( line_bundle_1.deg - line_bundle_2.deg - 1, 2  )
            return (0, 0, hom2_dim, 0)



class LocalP2:

    def spherical_twist(self, vec_bundle_1, vec_bundle_2):
        """
        Spherical twist method that handles four different cases based on the types of vec_bundle_1 and vec_bundle_2.
        """
        if isinstance(vec_bundle_1, LineBundle) and isinstance(vec_bundle_2, LineBundle):
            # Case 1: Both are LineBundle
            print("Case 1: Both vec_bundle_1 and vec_bundle_2 are LineBundle")
            return self.__sph_twist_LineBundles(vec_bundle_1,vec_bundle_2)

        elif isinstance(vec_bundle_1, LineBundle) and isinstance(vec_bundle_2, CotangentBundle):
            # Case 2: vec_bundle_1 is LineBundle, vec_bundle_2 is CotangentBundle
            print("Case 2: vec_bundle_1 is LineBundle, vec_bundle_2 is CotangentBundle")
            # Add your logic for this case here

        elif isinstance(vec_bundle_1, CotangentBundle) and isinstance(vec_bundle_2, LineBundle):
            # Case 3: vec_bundle_1 is CotangentBundle, vec_bundle_2 is LineBundle
            print("Case 3: vec_bundle_1 is CotangentBundle, vec_bundle_2 is LineBundle")
            # Add your logic for this case here

        elif isinstance(vec_bundle_1, CotangentBundle) and isinstance(vec_bundle_2, CotangentBundle):
            # Case 4: Both are CotangentBundle
            print("Case 4: Both vec_bundle_1 and vec_bundle_2 are CotangentBundle")
            # Add your logic for this case here

        else:
            raise TypeError("vec_bundle_1 and vec_bundle_2 must be instances of LineBundle or CotangentBundle.")



    def __sph_twist_LineBundles(self, line_bundle_1, line_bundle_2):
        '''
        use triangle  i^* i_* E -> E -> E x O(-3)[2]
        '''

        degree_diff = line_bundle_2.c1 - line_bundle_1.c1

        if degree_diff >= 0:
            rk0 = math.comb(degree_diff + 2, 2)
            rk1 = math.comb(degree_diff + 5, 2)

            vec0 = VectorBundle(rk0, rk0 * line_bundle_1.c1, rk0 * line_bundle_1.c1**2 / 2)
            vec1 = VectorBundle(rk1, rk1 * line_bundle_1.c1, rk1 * line_bundle_1.c1**2 / 2)

            complex1 = ChainComplex([vec1, vec0], -1)

            vec2 = VectorBundle(1, line_bundle_2.c1, line_bundle_2.c1**2 / 2)
            complex2 = ChainComplex([vec2], 0)

            chern2 = complex2.chernClass()
            chern1 = complex1.chernClass()

            ret_chern = ChernClass(chern2.ch0 - chern1.ch0, chern2.ch1 - chern1.ch1, chern2.ch2 - chern1.ch2)

            return ret_chern

        
        elif degree_diff < 0 and degree_diff > -3:
            rk0 = math.comb(degree_diff + 5, 2)

            vec1 = VectorBundle(rk0, rk0 * line_bundle_1.c1, rk0 * line_bundle_1.c1**2 / 2)
            complex1 = ChainComplex([vec1], -1)

            vec2 = VectorBundle(1, line_bundle_2.c1, line_bundle_2.c2)

            complex2 = ChainComplex([vec2], 0)

            chern2 = complex2.chernClass()
            chern1 = complex1.chernClass()

            ret_chern = ChernClass(chern2.ch0 - chern1.ch0, chern2.ch1 - chern1.ch1, chern2.ch2 - chern1.ch2)

            return ret_chern

        

if __name__ == "__main__":
    linebundle1 = LineBundle(-3)
    linebundle2 = LineBundle(-2)
    linebundle3 = LineBundle(-1)
    linebundle4 = LineBundle(0)
    linebundle5 = LineBundle(1)

    chaincomplex = ChainComplex([linebundle1, linebundle2, linebundle3, linebundle4], -5)

    print(chaincomplex)

    cotangent5 = CotangentBundle(5)
    print(cotangent5)

    twist = LocalP2().spherical_twist(linebundle4, linebundle5)
    print(twist)


    twist2 = LocalP2().spherical_twist(linebundle4, linebundle3)    
    print(twist2)