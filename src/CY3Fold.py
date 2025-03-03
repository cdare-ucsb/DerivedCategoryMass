import numpy as np
import cmath
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
import cmath
import plotly.graph_objects as go
from plotly.graph_objs import *
import plotly.utils
import json



def calcZ1(x, y, k, n):
    return cmath.exp(1j*(2*cmath.pi*k/n)) * (cmath.cos(x+y*1j)**(2/n))

def calcZ2(x, y, k, n):
    return cmath.exp(1j*(2*cmath.pi*k/n)) * (cmath.sin(x+y*1j)**(2/n))

def calcZ1Real(x, y, k, n):
    return (calcZ1(x, y, k, n)).real

def calcZ2Real(x, y, k, n):
    return (calcZ2(x, y, k, n)).real

def calcZ(x, y, k1_, k2_, n, a_):
    z1 = calcZ1(x, y, k1_, n)
    z2 = calcZ2(x, y, k2_, n)
    return z1.imag * math.cos(a_) + z2.imag*math.sin(a_)

def calcZ_alt(x,y,k1_, k2_, n, a_):
    z1 = calcZ1(x, y, k1_, n)
    z2 = calcZ2(x, y, k2_, n)
    return z1.imag * math.sin(a_) - z2.imag*math.cos(a_)



def generate_matplot_animation_1(filename, row=30, col=30, nframes=100, t_interval=175):
    
    n = 5

    # set param range
    x = np.linspace(0, math.pi/2, col)
    y = np.linspace(-math.pi/2, math.pi/2, row)
    x, y = np.meshgrid(x, y)

    fig = plt.figure(dpi=100)

    ax = fig.add_subplot(projection='3d')

    def update(t):
        ax.cla()

        for k1 in range(n):
            for k2 in range(n):

                    X = np.frompyfunc(calcZ1Real, 4, 1)(x, y, k1, n).astype('float32')
                    Y = np.frompyfunc(calcZ2Real, 4, 1)(x, y, k2, n).astype('float32')
                    Z1 = np.frompyfunc(calcZ, 6, 1)(x, y, k1, k2, n, t/10).astype('float32')
                    Z2 = np.frompyfunc(calcZ_alt, 6, 1)(x, y, k1, k2, n, t/10).astype('float32')

                    ax.plot_surface(X, Y, Z2, alpha=0.8, cmap=cm.afmhot)
        

        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_zlim(-2, 2)

    ani = FuncAnimation(fig = fig, func = update, frames = nframes, interval = t_interval)
    
    
    with open(filename, "w") as f:
        print(ani.to_html5_video(), file=f)




######### WARNING #########
# The following function generates an animation that is too large for the browser to handle
# It is recommended to use the generate_matplot_animation_1 function instead
def generate_plotly_animation_1(row=30, col = 30):
    
    frames = []

    fig = go.Figure()

    x = np.linspace(0, math.pi/2, col)
    y = np.linspace(-math.pi/2, math.pi/2, row)
    x, y = np.meshgrid(x, y)

    X1 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
    Y1 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
    Z1 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 0, 3, 0).astype('float32')

    X2 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
    Y2 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
    Z2 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 1, 3, 0).astype('float32')

    X3 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
    Y3 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
    Z3 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 2, 3, 0).astype('float32')
    
    X4 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
    Y4 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
    Z4 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 0, 3, 0).astype('float32')

    X5 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
    Y5 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
    Z5 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 1, 3, 0).astype('float32')

    X6 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
    Y6 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
    Z6 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 2, 3, 0).astype('float32')
    
    X7 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
    Y7 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
    Z7 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 0, 3, 0).astype('float32')

    X8 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
    Y8 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
    Z8 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 1, 3, 0).astype('float32')

    X9 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
    Y9 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
    Z9 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 2, 3, 0).astype('float32')


    fig.add_trace(go.Surface(x=X1, y=Y1, z=Z1, showscale=False))
    fig.add_trace(go.Surface(x=X2, y=Y2, z=Z2, showscale=False))
    fig.add_trace(go.Surface(x=X3, y=Y3, z=Z3, showscale=False))
    fig.add_trace(go.Surface(x=X4, y=Y4, z=Z4, showscale=False))
    fig.add_trace(go.Surface(x=X5, y=Y5, z=Z5, showscale=False))
    fig.add_trace(go.Surface(x=X6, y=Y6, z=Z6, showscale=False))
    fig.add_trace(go.Surface(x=X7, y=Y7, z=Z7, showscale=False))
    fig.add_trace(go.Surface(x=X8, y=Y8, z=Z8, showscale=False))
    fig.add_trace(go.Surface(x=X9, y=Y9, z=Z9, showscale=False))


    for alpha in range(1, 50):
        X1 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
        Y1 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
        Z1 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 0, 3, alpha).astype('float32')

        X2 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
        Y2 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
        Z2 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 1, 3, alpha).astype('float32')

        X3 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
        Y3 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
        Z3 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 2, 3, alpha).astype('float32')
        
        X4 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
        Y4 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
        Z4 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 0, 3, alpha).astype('float32')

        X5 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
        Y5 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
        Z5 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 1, 3, alpha).astype('float32')

        X6 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
        Y6 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
        Z6 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 2, 3, alpha).astype('float32')
        
        X7 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
        Y7 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
        Z7 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 0, 3, alpha).astype('float32')

        X8 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
        Y8 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
        Z8 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 1, 3, alpha).astype('float32')

        X9 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
        Y9 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
        Z9 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 2, 3, alpha).astype('float32')

        

        frame = go.Frame(data =[go.Surface(x=X1, y=Y1, z=Z1, showscale=False),
                                go.Surface(x=X2, y=Y2, z=Z2, showscale=False),
                                go.Surface(x=X3, y=Y3, z=Z3, showscale=False),
                                go.Surface(x=X4, y=Y4, z=Z4, showscale=False),
                                go.Surface(x=X5, y=Y5, z=Z5, showscale=False),
                                go.Surface(x=X6, y=Y6, z=Z6, showscale=False),
                                go.Surface(x=X7, y=Y7, z=Z7, showscale=False),
                                go.Surface(x=X8, y=Y8, z=Z8, showscale=False),
                                go.Surface(x=X9, y=Y9, z=Z9, showscale=False)
                                ])
        frames.append(frame)


    fig.update(frames=frames)
    fig.update_layout(
        updatemenus=[dict(
            type="buttons",
            showactive=False,
            buttons=[dict(label="Play",
                          method="animate",
                          args=[None, dict(frame=dict(duration=20, redraw=True),
                                            fromcurrent=True)]),

                    dict(label="Pause",
                          method="animate",
                            args=[[None],
                                   dict(frame=dict(duration=0, redraw=False),
                                         mode="immediate",
                                         transition=dict(duration=0))])]
        )]
    )

    
    fig.update_layout(showlegend=False)

    fig.show()


def generate_plotly_graph_1(row=30, col=30, return_json=False):
    n=4

    x = np.linspace(0, math.pi/2, col)
    y = np.linspace(-math.pi/2, math.pi/2, row)
    x, y = np.meshgrid(x, y)

    fig = go.Figure()

    for k1 in range(n):
        for k2 in range(n):
            
            X = np.frompyfunc(calcZ1Real, 4, 1)(x, y, k1, n).astype('float32')
            Y = np.frompyfunc(calcZ2Real, 4, 1)(x, y, k2, n).astype('float32')
            Z1 = np.frompyfunc(calcZ, 6, 1)(x, y, k1, k2, n, 0).astype('float32')
            Z2 = np.frompyfunc(calcZ_alt, 6, 1)(x, y, k1, k2, n, 0).astype('float32')

            fig.add_trace(go.Surface(x=X, y=Y, z=Z1, showscale=False, colorscale='blues'))

    if not return_json:
        fig.show()
        fig.write_html('k3_surface.html', full_html=False)
    else:
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    



def generate_plotly_sliders_1(row=30, col=30):

    n = 4

    fig = go.Figure()

    x = np.linspace(0, math.pi/2, col)
    y = np.linspace(-math.pi/2, math.pi/2, row)
    x, y = np.meshgrid(x, y)


    for step in np.arange(15.1, 15.3, 0.1):

        X1 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
        Y1 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
        Z1 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 0, 3, step).astype('float32')

        X2 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
        Y2 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
        Z2 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 1, 3, step).astype('float32')

        X3 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 0, 3).astype('float32')
        Y3 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
        Z3 = np.frompyfunc(calcZ, 6, 1)(x, y, 0, 2, 3, step).astype('float32')
        
        X4 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
        Y4 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
        Z4 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 0, 3, step).astype('float32')

        X5 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
        Y5 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
        Z5 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 1, 3, step).astype('float32')

        X6 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 1, 3).astype('float32')
        Y6 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
        Z6 = np.frompyfunc(calcZ, 6, 1)(x, y, 1, 2, 3, step).astype('float32')
        
        X7 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
        Y7 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 0, 3).astype('float32')
        Z7 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 0, 3, step).astype('float32')

        X8 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
        Y8 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 1, 3).astype('float32')
        Z8 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 1, 3, step).astype('float32')

        X9 = np.frompyfunc(calcZ1Real, 4, 1)(x, y, 2, 3).astype('float32')
        Y9 = np.frompyfunc(calcZ2Real, 4, 1)(x, y, 2, 3).astype('float32')
        Z9 = np.frompyfunc(calcZ, 6, 1)(x, y, 2, 2, 3, step).astype('float32')

        fig.add_trace(go.Surface(x=X1, y=Y1, z=Z1, showscale=False))
        fig.add_trace(go.Surface(x=X2, y=Y2, z=Z2, showscale=False))
        fig.add_trace(go.Surface(x=X3, y=Y3, z=Z3, showscale=False))
        fig.add_trace(go.Surface(x=X4, y=Y4, z=Z4, showscale=False))
        fig.add_trace(go.Surface(x=X5, y=Y5, z=Z5, showscale=False))
        fig.add_trace(go.Surface(x=X6, y=Y6, z=Z6, showscale=False))
        fig.add_trace(go.Surface(x=X7, y=Y7, z=Z7, showscale=False))
        fig.add_trace(go.Surface(x=X8, y=Y8, z=Z8, showscale=False))
        fig.add_trace(go.Surface(x=X9, y=Y9, z=Z9, showscale=False))

    # Make 10th trace visible
    fig.data[10].visible = True
    fig.data[11].visible = True
    fig.data[12].visible = True
    fig.data[13].visible = True
    fig.data[14].visible = True
    fig.data[15].visible = True
    fig.data[16].visible = True
    fig.data[17].visible = True
    fig.data[18].visible = True
    fig.data[19].visible = True

        # Create and add slider
    steps = []


    for i in range(len(fig.data) // 10):
        step = dict(
            method="update",
            args=[{"visible": [False] * len(fig.data)},
                {"title": "Slider switched to step: " + str(i)}],  # layout attribute
        )
        step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+1] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+2] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+3] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+4] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+5] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+6] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+7] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+8] = True  # Toggle i'th trace to "visible"
        step["args"][0]["visible"][i+9] = True  # Toggle i'th trace to "visible"

        steps.append(step)

    sliders = [dict(
        active=1,
        currentvalue={"prefix": "Frequency: "},
        pad={"t": 20},
        steps=steps
    )]

    fig.update_layout(
        sliders=sliders
    )

    fig.show()











if __name__ == "__main__":
    
    generate_plotly_graph_1(row=30, col=30)
    # generate_plotly_sliders_1(row=30, col = 30)