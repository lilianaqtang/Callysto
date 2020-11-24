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

    fig = go.Figure(data=[
        go.Mesh3d(
           x=[0, 0, x, x, 0, 0, x, x],
           y=[0, y, y, 0, 0, y, y, 0],
           z=[0, 0, 0, 0, z, z, z, z],

        colorbar_title='z',
        colorscale=[[0, 'gold'],
                    [0.5, 'mediumturquoise'],
                    [1, 'magenta']],
        # Intensity of each vertex, which will be interpolated and color-coded
        intensity = np.linspace(0, 1, 12, endpoint=True),
        intensitymode='cell',
        # i, j and k give the vertices of triangles
        i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
        j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
        name='y',
        showscale=False
    )
])

    fig.show()
    display(Markdown("Building faces which direction? (aspect)"))
    display(Markdown("<br> *Display the whole building, add aspect labels – allow rotation*"))

    display(Markdown("Window on N, E, S, W wall?"))
    display(Markdown("<br> *toggle for windows on each N, E, S, W wall*"))

    display(Markdown("Change size of window"))
    display(Markdown("<br> *(all windows are set to same, and only 0 or 1 per wall)*"))

    display(Markdown("W2W Ratio (window area : wall area)   (interactive calculation) m2  : (interactive calculation) m2 "))
    
    display(Markdown("Insulation (thickness of wall)"))
    display(Markdown("<br> *(Display only one wall)*"))

    display(Markdown("Insulation (glass specification)"))
    display(Markdown("<br> *(Display only one wall)*"))
    
    display(Markdown("Insulation (roof specification)"))
    display(Markdown("<br> *(Display building with shaded roof)*"))

    display(Markdown("#### <span style='color:blue'>Better Building Calculations"))
    display(Markdown("Cost to build (interactive calculation) $"))
    display(Markdown("<br> *blurb from architecture team to explain calculation  and why it matters*"))

    display(Markdown("Heat Loss (interactive calculation) kW"))
    display(Markdown("<br> *blurb from architecture team to explain calculation  and why it matters*"))

    display(Markdown("Heating Intensity (Peak Thermal Load)  (interactive calculation) W/m2"))
    display(Markdown("<br> *blurb from architecture team to explain calculation  and why it matters*"))
    display(Markdown("<br> *blurb from architecture team about what desirable aim is for a sustainable building*"))

