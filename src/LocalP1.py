import dash
from dash import dcc, html
import plotly.graph_objects as go
from dash.dependencies import Input, Output
import numpy as np

# Initialize Dash app
dash_app = dash.Dash(
        __name__,
        server=flask_app,
        routes_pathname_prefix='/dashapp/'
    )


# Generate evenly spaced grid data
grid_size = 100  # Number of points along each axis
x_values = np.linspace(-5, 5, grid_size)
y_values = np.linspace(-5, 5, grid_size)
xx, yy = np.meshgrid(x_values, y_values)
grid_points = np.column_stack((xx.ravel(), yy.ravel()))

# Function to create the figure with arrows
def create_figure(target_x=1, target_y=1):
    fig = go.Figure()

    # Add invisible scatter points
    fig.add_trace(go.Scatter(
        x=grid_points[:, 0], y=grid_points[:, 1],
        mode="markers",
        marker=dict(size=8, color="blue", opacity=0),  # Invisible markers
        name="Grid Points"
    ))

    # Calculate the third line's end point
    third_line_x = -1 + target_x
    third_line_y = target_y

    fourth_line_x = -2 + target_x
    fourth_line_y = target_y

    fifth_line_x = -3 + target_x
    fifth_line_y = target_y

    # Define shapes (lines) to be added to the figure
    shapes = [
        # First arrow (fixed)
        dict(type="line", x0=0, y0=0, x1=-1, y1=0,
             line=dict(color="black", width=3)),
        # Second arrow (moves on click)
        dict(type="line", x0=0, y0=0, x1=target_x, y1=target_y,
             line=dict(color="red", width=3)),
        # Third line (dynamic based on click)
        dict(type="line", x0=0, y0=0, x1=third_line_x, y1=third_line_y,
             line=dict(color="green", width=3)),
        dict(type="line", x0=0, y0=0, x1=fourth_line_x, y1=fourth_line_y,
             line=dict(color="blue", width=3)),
        dict(type="line", x0=0, y0=0, x1=fifth_line_x, y1=fifth_line_y,
             line=dict(color="purple", width=3)),
    ]

    # Define annotations (arrowheads) to be added to the figure
    annotations = [
        # Arrowhead for first arrow
        dict(x=-1, y=0, ax=0, ay=0, text="", arrowcolor="red", arrowsize=2, arrowwidth=3, showarrow=True),
        # Arrowhead for second arrow
        dict(x=target_x, y=target_y, ax=0, ay=0, text="", arrowcolor="red", arrowsize=2, arrowwidth=3, showarrow=True),
        # Arrowhead for third line
        dict(x=third_line_x, y=third_line_y, ax=0, ay=0, text="", arrowcolor="green", arrowsize=2, arrowwidth=3, showarrow=True),
    ]

    # Update figure layout with shapes and annotations
    fig.update_layout(
        xaxis=dict(range=[-5, 5], zeroline=True),
        yaxis=dict(range=[-5, 5], zeroline=True),
        title="Interactive Arrows with Invisible Grid Points",
        shapes=shapes,
        annotations=annotations
    )

    return fig

# Dash layout
app.layout = html.Div([
    html.H3("Click on a grid point to move the second and third arrows"),
    dcc.Graph(id="interactive-graph", figure=create_figure()),
])

# Callback to update the arrows' positions
@app.callback(
    Output("interactive-graph", "figure"),
    Input("interactive-graph", "clickData"),
)
def update_arrows(click_data):
    if click_data:
        # Retrieve the clicked coordinates
        x_click = click_data["points"][0]["x"]
        y_click = click_data["points"][0]["y"]
    else:
        # Default position if no click has occurred
        x_click, y_click = 1, 1

    # Update the figure with the new arrow positions
    return create_figure(x_click, y_click)

# Run the Dash app
if __name__ == "__main__":
    app.run_server(debug=True)
