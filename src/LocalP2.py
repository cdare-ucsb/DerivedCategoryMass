from .ChernCharacter import ChernCharacter
import matplotlib.pyplot as plt
import math
import numpy as np
import plotly.graph_objects as go



# class StabilityCondition():

#     def __init__(self, s, q):
#         self.s = s
#         self.q = q


#     def plot_central_charge(self, objects, labels=None):
#         """
#         Plots the central charge on the complex plane and draws a line from the origin (0,0) to each complex number.

#         Args:
#             central_charges (list of complex): A list of complex numbers to plot.
#             labels (list of str, optional): Labels corresponding to each complex number.
#         """

#         central_charges = [obj.central_charge(self.s, self.q) for obj in objects]

#         plt.figure(figsize=(6, 6))
        
#         for i, z in enumerate(central_charges):

#             stab_color = 'blue' if objects[i].is_semistable(self.s, self.q) else 'red'
#             line_color = 'cyan' if objects[i].is_semistable(self.s, self.q) else 'red'

#             # Plot point
#             plt.scatter(z.real, z.imag, color=stab_color, marker='o')
            
#             # Draw line from (0,0) to the complex number
#             plt.plot([0, z.real], [0, z.imag], linestyle='-', color=line_color, alpha=0.7)

#             # Add labels if provided
#             if labels and i < len(labels):
#                 plt.text(z.real, z.imag, f" {labels[i]}", fontsize=12, verticalalignment='bottom')

#         # Draw x-axis and y-axis
#         plt.axhline(0, color='black', linewidth=1)
#         plt.axvline(0, color='black', linewidth=1)
#         plt.grid(True, linestyle="--", alpha=0.6)

#         plt.xlabel("Real")
#         plt.ylabel("Imaginary")
#         plt.title("Central Charge on Complex Plane")

#         plt.show()


class LePotier():

    def __init__(self, granularity=5, width=5):
        self.granularity = granularity
        self.width = width

        lower_bound = -1*self.width * 2**self.granularity
        upper_bound = self.width * 2**(self.granularity ) + 1

        self.boundary_points = []

        for i in range(lower_bound, upper_bound):

            self.boundary_points.append(self._e_plus(i, granularity))
            self.boundary_points.append(self._e_left(i, granularity))
            self.boundary_points.append(self._e_right(i, granularity))
        
        self.boundary_points = sorted(self.boundary_points, key=lambda x: x[0])
        


    def _get_dyadic_character(self, p, m):

        if not isinstance(p, int) or not isinstance(m, int) or not m >= 0:
            raise ValueError("Input data must be integers, with m >= 0")

        
        if p % int(2**m) == 0:
            d = p // int(2**m)
            return ChernCharacter(1, int(d), float( d**2 ) / 2)
        elif p % 4 == 3:
            return 3 * self._get_dyadic_character(p+1, m).ch0 * self._get_dyadic_character(p-1, m) - self._get_dyadic_character(p-3, m)
        elif p % 4 == 1:
            return 3 * self._get_dyadic_character(p-1, m).ch0 * self._get_dyadic_character(p+1, m) - self._get_dyadic_character(p+3, m)
        elif p % 4 == 2:
            new_p = p // 2
            return self._get_dyadic_character(new_p, m-1)
        else:
            new_p = p // 4
            return self._get_dyadic_character(new_p, m-2)
        
    def _e_reg(self, p, m):
        chern = self._get_dyadic_character(p, m)
        return ( float(chern.ch1) / chern.ch0, chern.ch2 / chern.ch0 )
        
    def _e_plus(self, p, m):
        chern = self._get_dyadic_character(p, m)
        return ( float(chern.ch1) / chern.ch0, chern.ch2 / chern.ch0 - float(1 / (chern.ch0)**2) )
    
    def _e_left(self, p, m):
        x_1, y_1 = self._e_plus(p, m)
        x_2, y_2 = self._e_reg(p-1, m)

        m = float(y_2 - y_1) / (x_2 - x_1)

        x_3 = m + math.sqrt( m**2 - 2*m*x_1 + 2*y_1 + 1 )
        if x_1 <= x_3 <= x_2 or x_2 <= x_3 <= x_1:
            y_3 = m * x_3 + (y_1 - m * x_1)
            return (x_3, y_3)
        else:
            x_3 = m - math.sqrt( m**2 - 2*m*x_1 + 2*y_1 + 1 )
            y_3 = m * x_3 + (y_1 - m * x_1)
            return (x_3, y_3)
    
    def _e_right(self, p, m):
        x_1, y_1 = self._e_plus(p, m)
        x_2, y_2 = self._e_reg(p+1, m)

        m = float(y_2 - y_1) / (x_2 - x_1)

        x_3 = m + math.sqrt( m**2 - 2*m*x_1 + 2*y_1 + 1 )
        if x_1 <= x_3 <= x_2 or x_2 <= x_3 <= x_1:
            y_3 = m * x_3 + (y_1 - m * x_1)
            return (x_3, y_3)
        else:
            x_3 = m - math.sqrt( m**2 - 2*m*x_1 + 2*y_1 + 1 )
            y_3 = m * x_3 + (y_1 - m * x_1)
            return (x_3, y_3)
        


    def is_above_curve(self, x, y):
        return y > self.curve_estimate(x)
    
    def curve_estimate(self, x):
        x_values = [p[0] for p in self.boundary_points]

        # Check if x is within the range of the curve
        if x < x_values[0] or x > x_values[-1]:
            raise ValueError("x is outside the range of the curve")

        # Find the segment containing x
        for i in range(len(self.boundary_points) - 1):
            x1, y1 = self.boundary_points[i]
            x2, y2 = self.boundary_points[i + 1]

            if x1 <= x <= x2:
                # Perform linear interpolation to find y-interp at x
                y_interp = y1 + (y2 - y1) * ((x - x1) / (x2 - x1))
                
                return y_interp

        raise ValueError("x is outside the range of the curve")


        
    def plot_drezet_le_potier(self):

        lower_bound = -1*self.width * 2**self.granularity
        upper_bound = self.width * 2**(self.granularity ) + 1

        # Create figure
        fig = go.Figure()   

        for i in range(lower_bound, upper_bound):

            x1, y1 = self._e_plus(i, self.granularity)
            x2, y2 = self._e_left(i, self.granularity)
            x3, y3 = self._e_right(i, self.granularity)

            fig.add_trace(go.Scatter(
                x=[x1, x2], 
                y=[y1, y2], 
                mode='lines+markers',
                line=dict(color='blue', width=4),
                marker=dict(size=6),
                showlegend=False  # Hide legend for individual line segments
            ))

            fig.add_trace(go.Scatter
            (
                x=[x1, x3], 
                y=[y1, y3], 
                mode='lines+markers',
                line=dict(color='blue', width=4),
                marker=dict(size=6),
                showlegend=False  # Hide legend for individual line segments
            ))
        # Ensure equal scaling of x and y axes
        fig.update_layout(
            xaxis=dict(scaleanchor="y"),  # Locks the aspect ratio
            yaxis=dict(scaleanchor="x")   # Ensures square grid
        )    

        config = dict({'scrollZoom': True})
        # Show plot
        fig.show(config=config)


    def plot_drezet_le_potier3d(self):

        lower_bound = -1*self.width * 2**self.granularity
        upper_bound = self.width * 2**(self.granularity ) + 1

        # Create figure
        fig = go.Figure()   

        for i in range(lower_bound, upper_bound):

            x1, y1 = self._e_plus(i, self.granularity)
            x2, y2 = self._e_left(i, self.granularity)
            x3, y3 = self._e_right(i, self.granularity)

            fig.add_trace(go.Scatter3d(
                x=[x1, x2], 
                y=[y1, y2], 
                z=[0,0],
                mode='lines+markers',
                line=dict(color='blue', width=5),
                marker=dict(size=6),
                showlegend=False  # Hide legend for individual line segments
            ))

            fig.add_trace(go.Scatter3d
            (
                x=[x1, x3], 
                y=[y1, y3], 
                z=[0,0],
                mode='lines+markers',
                line=dict(color='blue', width=5),
                marker=dict(size=6),
                showlegend=False  # Hide legend for individual line segments
            ))
        # Ensure equal scaling of x and y axes
        # fig.update_layout(
        #     xaxis=dict(scaleanchor="y"),  # Locks the aspect ratio
        #     yaxis=dict(scaleanchor="x")   # Ensures square grid
        # )    

        # Show plot
        fig.show()




         
        


if __name__ == "__main__":
    DLP = LePotier(granularity=3, width=5)
    DLP.plot_drezet_le_potier()

    print(DLP.curve_estimate(0.48950))
    


