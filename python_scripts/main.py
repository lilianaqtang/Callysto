from architecture_script import construct_building
import plotly
import plotly.graph_objs as go
from ipywidgets import interact, interact_manual
from IPython.display import display
import ipywidgets as widgets
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ipywidgets import interact, interact_manual, widgets, Layout, VBox, HBox, Button
from IPython.display import display, Javascript, Markdown, HTML, clear_output
import pandas as pd
import plotly.express as px 
import plotly.graph_objects as go
from IPython.display import HTML
import codecs
import re


def choose_subject(option):
    intro_text=widgets.HTML(value=" ")
    if option == 'Architecture':
         display(Markdown("### <span style='color:blue'>ARCHITECTURE FOCUS SELECTED"))
         display(Markdown("<br> *Introductory blurb about architecture from architecture team. Use a graphic from Pulp team.*"))
         display(Markdown("#### <span style='color:blue'>Variables in architecture"))
         display(Markdown("Select building size:"))
         interact_manual(construct_building, 
         x=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Length: X(m)',style=style), 
         y=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Width: Y(m)',style=style),
         z=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Height: Z(m)',style=style));

    elif option == 'Biology':
         display(Markdown("### <span style='color:blue'>BIOLOGY FOCUS SELECTED"))

    elif option == 'Computer Science':
         display(Markdown("### <span style='color:blue'>COMPUTER SCIENCE FOCUS SELECTED"))
    elif option == 'Math':
         display(Markdown("### <span style='color:blue'>MATH FOCUS SELECTED"))
    elif option == 'Physics':
         display(Markdown("### <span style='color:blue'>PHYSICS FOCUS SELECTED"))

def construct_building(x, y, z):
    """
    This function calculates the floor area after the student inputs the desired length and width for the building

    Args:
        x(float): length of the building
        y(float): width of the building

    Returns:
        floor area = x * y
    """

    print("Floor area (footprint) = length x width = {0:.{1}f}m\u00b2".format(x * y, 0))

    print("Building volume = length x width x height = {0:.{1}f}m\u00b3".format(x * y * z, 0))

    # Configure Plotly to be rendered inline in the notebook.
    plotly.offline.init_notebook_mode()

    # Configure the trace.
    trace = go.Scatter3d(
        x=[0, x],
        y=[0, y],
        z=[0, z],

        mode='markers',
        marker={
            'size': 10,
            'opacity': 0.8,
            'color': 'olive'
        },
        text=['Your Floor Areas', 'Your Wall Areas', 'Roof']

    )

    # Configure the layout.
    layout = go.Layout(
        margin={'l': 0, 'r': 0, 'b': 0, 't': 0},
        plot_bgcolor='rgb(12,163,135)'
    )

    data = [trace]

    plot_figure = go.Figure(data=data, layout=layout)

    # Render the plot.
    plotly.offline.iplot(plot_figure)