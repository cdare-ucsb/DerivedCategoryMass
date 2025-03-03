import numpy as np
import math

import plotly.graph_objects as go
from plotly.graph_objs import *


class K3Surface:

    def __init__(self, degree):
        self.degree = degree


    def create_mukai_vectors(self, granularity=15):
        """
        Method to create the mukai vectors of a spherical vector bundle on a K3 surface of degree d.


        Parameters
        ----------
        granularity : int
            The granularity of the mukai vectors to be computed. This roughly indicates how close the mukai vectors should get to integer points
            on the real axis

        Returns
        -------
        mukai_vectors : list
            A list of tuples representing the mukai vectors of the spherical vector bundles
        """

        mukai_vectors = []

        for n in range(-granularity, granularity):

            numerator = self.degree*n**2 + 1
            divisors = self._npDivs(numerator)
            for div in divisors:
                r = div
                s = numerator // div

                mukai_vectors.append((r, n, s))

        return mukai_vectors

    def plot_alpha_beta_plane(self):

        x_vals = []
        y_vals = []

        mukai_vectors = self.create_mukai_vectors()

        for r, l, s in mukai_vectors:
            alpha, beta = self._convert_to_alpha_beta(r, l, s)
            x_vals.append(alpha)
            y_vals.append(beta)

            x_vals.append(alpha+1)
            y_vals.append(beta)

            x_vals.append(alpha-1)
            y_vals.append(beta)
        
        fig = go.Figure(data=[go.Scatter(x=x_vals, y=y_vals, mode='markers')])

        # Add vertical lines from x-axis to each point
        for x_val, y_val in zip(x_vals, y_vals):
            fig.add_shape(
                type="line",
                x0=x_val,
                y0=0,
                x1=x_val,
                y1=y_val,
                line=dict(
                    color="LightSeaGreen",
                    width=1
                )
            )

        fig.update_layout(
            title="Mukai Vectors",
            xaxis_title="α",
            yaxis_title="β",
            showlegend=False
        )
        fig.update_xaxes(range=[-1.6, 1.6])
        fig.update_yaxes(range=[-0.2, math.sqrt(self.degree)+0.2])

        # fig.show(config={'scrollZoom': True})
        fig.show()



    def _npDivs(self, N):
        """
        Helper method to compute the divisors of a number N.

        Parameters
        ----------
        N : int
            The number whose divisors are to be computed

        Returns
        -------
        divs : numpy.ndarray
            An array of divisors of N
        """

        divs = np.arange(1,int(N**0.5)+1)  # potential divisors up to √N
        divs = divs[N % divs==0]           # divisors
        comp = N//divs[::-1]               # complement quotients
        return np.concatenate((divs,comp[divs[-1]==comp[0]:])) # combined


    def _convert_to_alpha_beta(self, r, l, s):
        """
        Helper method to compute the coordinates of a hole in the (α,β) plane corresponding to a mukai vector.


        Suppose (r, lH, s) is the mukai vector of a spherical vector bundle on a degree d K3 surface (H^2 = 2d). Then
        # dl^2 - rs = -1. For a corresponding hole in the (α,β) plane, we must have the real and imaginary part of 
        # < exp(αΗ + ι βΗ) , (r, lH, s) > vanish. This gives

        0 = <(1, αH, (α^2 - β^2)d), (r, lH, s)> = 2dal - s - rd(α^2 - β^2)
        0 = <(0, βΗ, 2αβd), (r, lH, s)> = 2βdl - 2rdαβ = 2dβ(l - rα)
        
        The second equation gives β = 0 or l = rα. Since βΗ is necessarily ample, we must have β > 0.
        Therefore, l = rα. Substituting into the first equation simplifies to

        dl^2 - rs + dr^2 β^2 = 0

        Using dl^2 - rs = -1, we find that β = 1/r√d. 

        Parameters
        ----------
        r : int
            The first component of the mukai vector
        l : int
            The second component of the mukai vector
        s : int
            The third component of the mukai vector

        Returns
        -------
        alpha : float
            The x-coordinate of the hole in the (α,β) plane
        beta : float
            The y-coordinate of the hole in the (α,β) plane
        """


        alpha = float(l) / r
        beta = float(1 / r*math.sqrt(self.degree))

        return alpha, beta
 
    
if __name__ == "__main__":

    k3 = K3Surface(1)
    k3.plot_alpha_beta_plane()

    

    
    

        




