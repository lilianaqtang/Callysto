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
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch, Rectangle
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import random
import copy
import matplotlib.pyplot as plt
import seaborn as sns
from ipywidgets import interact

def choose_subject(option):
    intro_text=widgets.HTML(value=" ")
    if option == 'Architecture':
         display(Markdown("### <span style='color:blue'>ARCHITECTURE FOCUS SELECTED"))
         display(Markdown("In this notebook we explore the principles of modelling building envelopes and calculating key energy factors. The calculations performed inform sustainable design by considering how changes to a building’s variables such as its dimensions, floor area and volume,  its window to wall surface ratio, etc., impact the buildings need for energy and ability to retain energy. "))
         display(Markdown("Experiment with changing the dimensions of a simple building and observe their impact on the building’s area and volume, energy use and the building’s cost. Notice even in this simple model there are a seemingly endless number of ways to design a building, with each design differing in its need for and ability to retain energy as well as its cost. "))
         display(Markdown("#### <span style='color:blue'><u>Variables in architecture</u>"))
         display(Markdown("<span style='color:red'><u>Select building dimensions:</u>"))
         interact_manual(construct_building,
         x=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Length: X(m)',style=style),
         y=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Width: Y(m)',style=style),
         z=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Height: Z(m)',style=style));
         display(Markdown("<span style='color:red'><u>Select building orientation by clicking & holding the building to the direction you wish:</u>"))
         display(Markdown("<span style='color:red'><u>Select glazing:</u>"))
         display(Markdown("<br> Change size of window *(all windows are set to same, and only 0 or 1 per wall)*"))
         interact_manual(choose_direction,
             direction=widgets.RadioButtons(
                 options=['North(N)', 'South(S)', 'West(W)', 'East(E)'],
                 value='West(W)',
                 description='Window on N, E, S, W wall?',
                 disabled=False,
             ),
            x=widgets.FloatSlider(min=1, max=20, step=1, value=1, description='Length: X(m)', style=style),
            y=widgets.FloatSlider(min=1, max=20, step=1, value=1, description='Width: Y(m)', style=style),
            z=widgets.FloatSlider(min=1, max=20, step=1, value=1, description='Height: Z(m)', style=style),
            x1=widgets.FloatSlider(min=1, max=10, step=1, value=1, description='Length of the window: X(m)', style=style),
            y1=widgets.FloatSlider(min=1, max=10, step=1, value=1, description='Height of the window: Y(m)', style=style),
            z1=widgets.FloatSlider(min=1, max=10, step=1, value=1, description='Thickness of the window: Z(cm)', style=style));

         display(Markdown("<span style='color:red'><u>Insulation (thickness of wall):</u>"))
         #display(Markdown("<br> *(Display only one wall) - Image from Pulp team*"))
         display(HTML("""<img src="./images/wall.jpg" width="40%" height="40%">"""))

         display(Markdown("<span style='color:red'><u>Insulation (glass specification):</u>"))
         #display(Markdown("<br> *(Display only one wall) - Image from Pulp team*"))
         display(HTML("""<img src="./images/glass_insulation.jpg" width="40%" height="40%">"""))

         display(Markdown("<span style='color:red'><u>Insulation (roof specification):</u>"))
         #display(Markdown("<br> *(Display building with shaded roof) - Image from Pulp team*"))
         display(HTML("""<img src="./images/roof.jpg" width="40%" height="40%">"""))

         display(Markdown("#### <span style='color:blue'><u>Better Building Calculations<u>"))
         display(Markdown("#### <span style='color:red'>Cost to build ($):"))
         display(Markdown("<br> The cost of a building is important to every client, and there are lots of different considerations that need be understood, so that they can be weighed up and decided upon. For this example, the investigation is looking at how a more insulated building could help to reduce the size of the heating system. Both of these have impact in capital cost and operating costs for heating energy bills (natural gas, electricity, wood, oil, etc)."))

         display(Markdown("#### <span style='color:red'>Heat Loss (kW): "))
         display(Markdown("<br> The heat loss of a building is calculated by measuring the surface area (or envelope) of the building walls, windows, roof and floor that is exposed to the outside environment. Each component of the building has insulative properties that will slow down the amount of heat that is transmitted (or lost) through to outside. Advanced calculations also include heat lost by air leakages (cold drafts) and heat required to warm up air for ventilation purposes to keep occupant healthy with good air quality."))

         display(Markdown("#### <span style='color:red'>Heating Intensity (Peak Thermal Load) (W/m\u00b2): "))
         display(Markdown("<br> The total heat loss for the building is calculated by adding up the heat losses from each building envelope component. This number vary considerably from project to project and is dependent upon size of building, WWR and level of insulation standards. In order to give a simple comparison of relative energy efficiency, a heating intensity calculation is performed by dividing the total heat loss by the habitable building floor area. This calculation takes the whole sum of heating requirements and averages it out over the floor area. The lower the peak thermal load the better the building, and the lower costs to operate."))
         display(Markdown("<b> The desirable aim is for a sustainable building:"))
         display(Markdown("The goal in designing a sustainable building should be that the envelope is capable of keeping occupants comfortable and healthy in the middle of winter without the need for extra heating. Interestingly, Canada has been involved in research into how these buildings should be designed and <a href='https://www.src.sk.ca/blog/closer-look-saskatchewan-conservation-house-and-four-others' target='_blank'>some built examples </a> show that this can be achieved even in a cold winter prairie climate."))
         display(Markdown("However, despite proof that this technology works, it is unfortunate that many buildings that are built in Canada today still require large amounts of energy to keep them comfortable and this results in greenhouse gas emissions, climate change as well as high costs for building owners."))
         display(Markdown("The goal of this learning tool is provide a first step in education of how simple steps can make significant results in lower the need for heating energy and improving the global environment and comfort for building occupants."))

    elif option == 'Computer Science':
         display(Markdown("### <span style='color:blue'>SELECTED: COMPUTER SCIENCE FOCUS ON DESIGN FOR FORM AND FUNCTION"))
         display(Markdown('A Genetic algorithm (GA) is a set of instructions (algorithm) to a computer, inspired by the evolutionary process of natural selection, that identify optimal solutions to fit a goal (fitness function).'))
         display(HTML("""A GA is used in optimization problems, when the possible solution(s) are complex and numerous and the best one needs to be found. E.g. designs for buildings; <a href="https://en.wikipedia.org/wiki/Evolved_antenna" target="_blank">NASAs ST5 radio antenna design</a>; models of human movement. The best solution (fittest) is defined by the problem being solved e.g. building, cost, desirability and sustainability to build and operate; antenna range of transmission; similarity to human movement."""))
         display(Markdown('The GA takes a set of options (each option has a chromosome  which is a phenotype descriptor),  for solutions to a problem (population) competes the options against the fitness function to find a winner(s). The winner(s) and new options, created through a process of change of their description (mutation and crossover), create a next set of options  (a generation) and the process is repeated. The process stops (terminates) after the fitness criterion is satisfied OR after a number of tries (iterations).'))
         display(Markdown(
             "#### <span style='color:red'>Challenge: Design a GA to find routes from top left corner to bottom right corner on the map, avoiding the blocks.  We want to minimize the number of steps in the route and the route’s distance."))
         display(Markdown("#### <span style='color:blue'>First let's Initialize the Map!"))
         display(Markdown('This map represents the space we are going to look for our route in. It’s a square space with blocks that the route must avoid. We can define how large the space is and how densely the blocks are arranged. Density is controlled by two factors, the distance between the blocks (closeness) and blocks to no-blocks ratio (spareness)'))
         interact_manual(initialize_map,
                         sparseness_of_map=widgets.FloatSlider(min=0.05, max=0.99, step=0.01, value=0.95, description='Spareness: ', style=style),
                         size_of_map=widgets.IntSlider(min=1, max=2000, step=1, value=1000, description='Size: ', style=style),
                         number_of_groups=widgets.IntSlider(min=1, max=9, step=1, value=1, description='Closeness of Blocks: ', style=style));
         #display(Markdown("#### <span style='color:blue'>2. Let's Put Everything Together to Design the Best Routes in the Map!"))
         display(Markdown("After initializing the map, we need to apply our genetic algorithm to solve the route finding problem. The steps to do this are:"))
         display(Markdown("I.  Build the initial population: The process begins with a set of individuals which is called a Population. Each individual is a solution to the problem we want to solve. "))
         display(Markdown("II. Apply the fitness function: The fitness function determines how fit any individual solution is to compete with other individuals. Applying the function to the individual gives  it a fitness score. The probability that an individual will be selected as a winner is based on its fitness score."))
         display(Markdown("III. Make the winner(s) selection: The idea of selection phase is to select the fittest individuals from a population and let them pass their genes to the next generation."))
         display(Markdown("IV. Change up the population:"))
         display(Markdown("* Crossover: Crossover is the most significant phase in a genetic algorithm. For each pair of parents to be mated, a crossover point is chosen at random from within the genes. Offspring are created by exchanging the genes of parents among themselves until the crossover point is reached."))
         display(Markdown("* Mutation: In certain new offspring formed, some of their genes can be subjected to a mutation with a low random probability. This implies that some of the bits in the bit string can be flipped."))
         display(Markdown("The algorithm terminates if the population has converged (does not produce offspring which are significantly different from the previous generation). Then it is said that the genetic algorithm has provided a set of solutions to our problem."))
         display(Markdown("#### <span style='color:blue'>2. Let's Put Everything Together to Design the Best Routes for the Map!"))
         interact_manual(solve_map,
                         sparseness_of_map=widgets.FloatSlider(min=0.05, max=0.99, step=0.01, value=0.95, description='Spareness: ', style=style),
                         size_of_map=widgets.IntSlider(min=1, max=2000, step=1, value=1000, description='Size: ', style=style),
                         number_of_groups=widgets.IntSlider(min=1, max=9, step=1, value=1, description='Closeness of Blocks: ', style=style),
                         population_size=widgets.IntSlider(min=1, max=50, step=1, value=30, description='Population Size: ', style=style),
                         number_of_iterations=widgets.IntSlider(min=100, max=5000, step=1, value=1000, description='# Iterations: ', style=style),
                         number_of_couples=widgets.IntSlider(min=1, max=10, step=1, value=9, description='Crossover factor: ', style=style),
                         number_of_winners_to_keep=widgets.IntSlider(min=1, max=5, step=1, value=2, description='# Winners: ', style=style),
                         mutation_probability=widgets.FloatSlider(min=0.01, max=1, step=0.01, value=0.05, description='Mutation Probability: ', style=style));
    elif option == 'Math':
         display(Markdown("### <span style='color:blue'>MATH FOCUS SELECTED"))
         display(Markdown("__The Evolution of Stick Ungulates__"))
         display(Markdown("In this Math Lab we will use mathematical modelling to explore the process of adaptation through Natural Selection.  Before completing this lab, students should be familiar with our previously developed model for Natural Selection acting at a single diploid gene locus with two alleles. We will use this model to consider how populations of stick ungulates evolve under different selective pressures, created by different environments.  "))
         display(Markdown("*What is a stick ungulate you ask?*"))
         display(HTML("""<a href="https://en.wikipedia.org/wiki/Ungulate" target="_blank">Ungulates</a> are large hooved mammals that come in both an even toed
and odd toed variety. In nature, odd toed ungulates include mammals such as horses and rhinoceroses
and tapirs, while even toed ungulates include cows, deer, pronghorns, okapi and giraffes, to name a few.
In the Jupyter Notebook environment we will study stick ungulates, which are pictured below. The nice
thing about stick ungulates is they possess a single heritable trait, height, which is determined by the
expression of a single gene locus at which there are only ever 2 alleles present in the population, allele
and allele . All other characteristics of a stick ungulate are fixed relative to height. Specifically, the
proportions of a stick ungulate are as follows:"""))
         display(Markdown("* Hind legs 0.4 X height"))
         display(Markdown("* Fore legs are 0.5 X height"))
         display(Markdown("* Body length 0.4 X height"))
         display(Markdown("* Neck length 0.5 X height"))
         display(Markdown("* Head length varies"))
         display(HTML("""<img src="./images/math_figure1.jpg">"""))
         display(Markdown("*Figure 1. A stick ungulate. Stick ungulates are found in Jupyter Notebooks.  Ranging in heights from 1m≤h≤3m, stick ungulates have fixed body proportions as shown above.*"))
         display(Markdown("Stick ungulates inhabit environments within Jupyter notebooks. The survival and reproduction of stick ungulates is dependent upon the resources a, b and c which can vary among 3 types of environments; impacting the fitness of stick ungulates of differing heights.  Specifically, the fitness of stick ungulates in the three environments is given by:"))
         display(HTML("""<img src="./images/table1.jpg">"""))
         display(Markdown("Your task is to observe how a general population of stick ungulates evolves over the course of time, adapting to the three environments as a result of natural selection. "))
         display(Markdown("#### <span style='color:blue'>Lab Instructions"))
         display(Markdown("In this lab, you will run a number of trials.  "))
         display(Markdown("At the beginning of “Trial one” you will define three environments by specifying the resource values of a, b and c.  Next, you will introduce a population of stick ungulates into each environment. The height of the stick ungulates will be determined by two alleles A and a.  By specifying the heights that are determined by each allele, and the initial frequency of each allele you will create a population of ungulates of varying heights.  You will then track the evolution of the stick ungulates in each environment for 10 generations. At the end of the trial, you will note how the population changed in response to its environment documenting the change in the frequency of each allele, the change in frequency of each allele the change in the average fitness of the population and the change in the average height of the population."))
         display(Markdown("In the second trial, a “mutation event” occurs in which you are allowed to respecify the heights that are determined by each allele. Note that the frequency of each allele remains unchanged from the end of Trial 1.  With new heights specified for each allele, you will again track the evolution of the stick ungulates in each environment for 10 generations. Note that no changes to your three environments occur, which means that the selective pressures observed in each environment are the same as in Trial 1. At the end of the trial, you will note how the population changed in response to its environment documenting the change in the frequency of each allele, the change in frequency of each allele the change in the average fitness of the population and the change in the average height of the population."))
         display(Markdown("In the third trial, a second “mutation event” occurs in which you are allowed to respecify the heights that are determined by each allele. Note that the frequency of each allele remains unchanged from the end of Trial 2.  With new heights specified for each allele, you will again track the evolution of the stick ungulates in each environment for 10 generations. Note that no changes to your three environments occur, which means that the selective pressures observed in each environment are the same as in Trial 1. At the end of the trial, you will note how the population changed in response to its environment documenting the change in the frequency of each allele, the change in frequency of each allele the change in the average fitness of the population and the change in the average height of the population."))
         display(Markdown("At the end of Trial 3, you will be asked to report how your population of stick ungulates changed over time."))

         display(Markdown("#### <span style='color:blue'>Instructions"))
         display(Markdown("Students may want to work in small groups. Each student should complete the lab on their own.  However, students may wish to coordinate themselves in small groups as this will allow them to run a broader range of experiments. If this approach is taken, then students should devise an experimental strategy at the beginning of each trial. "))
         display(Markdown("**Trial 1.**"))
         display(Markdown("**Modelling the Environment**"))
         display(Markdown("In this section, you will define 3 different environments by specifying the resource values of a, b and c using the sliders below. Once you set the slider for the resource values of a, b and c the corresponding fitness function for that environment will be plotted immediately below. Use the graph of the fitness function to refine your selection of the resource values of a, b and c for each environment. In defining the three environments try to create three very different environments as reflected by the shape of the fitness function.  For example, you may wish to study a case in which there is selection for an extreme height, a case where there is selection for an intermediate height and a case where there is selection against an intermediate height."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         interact_manual(fitness_function,
                         a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2081', style=style),
                         b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2081', style=style),
                         c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2081', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         interact_manual(fitness_function,
                    a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2082', style=style),
                    b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2082', style=style),
                    c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2082', style=style));
         display(Markdown("### <span style='color:blue'>Environment 3"))
         interact_manual(fitness_function,
                         a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2083', style=style),
                         b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2083', style=style),
                         c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2083', style=style));
         display(Markdown("The environments you specify will remain unchanged over the course of the three trials."))
         display(Markdown("If you are working with other students, then coordinate amongst yourselves so to ensure you explore as many variations of environments as possible."))
         display(Markdown("**Modelling the Population**"))         
         display(Markdown("In this section, you will define the characteristics of an initial population of stick ungulates that will inhabit each of your three environments. The height of the stick ungulates will be determined by two alleles and . Using the appropriate sliders below specify the heights that are determined by each phenotype. Next the appropriate sliders below specify the initial frequency of each allele. Once you have set the heights that are determined by each phenotype and the initial frequency of each allele a table of detailing your population will be outputted. The table will identify the three genotypes in your population, their frequency, corresponding heights and fitness in each of the three environments. Population averages will also be reported."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));

         display(Markdown("### <span style='color:blue'>Environment 3"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("**Modelling Natural Selection**"))
         display(Markdown("In this section you will track the evolution of your stick ungulates as they evolve in each environment over 10 generations. The code to do so has already been set up and all you have to do is run it.  For the first generation, the code will output the various calculations so that you can follow the steps of our model.  For generations 2 to 10, the code will simply output a table that reports the frequency of each allele, the frequency of each phenotype, the average fitness of the population and the average height of the population."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 3"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));

         display(Markdown("**Trial 2.**"))
         display(Markdown("**Modelling the Environment**"))
         display(Markdown("In this section, you will define 3 different environments by specifying the resource values of a, b and c using the sliders below. Once you set the slider for the resource values of a, b and c the corresponding fitness function for that environment will be plotted immediately below. Use the graph of the fitness function to refine your selection of the resource values of a, b and c for each environment. In defining the three environments try to create three very different environments as reflected by the shape of the fitness function.  For example, you may wish to study a case in which there is selection for an extreme height, a case where there is selection for an intermediate height and a case where there is selection against an intermediate height."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         interact_manual(fitness_function,
                         a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2081', style=style),
                         b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2081', style=style),
                         c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2081', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         interact_manual(fitness_function,
                    a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2082', style=style),
                    b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2082', style=style),
                    c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2082', style=style));
         display(Markdown("### <span style='color:blue'>Environment 3"))
         interact_manual(fitness_function,
                         a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2083', style=style),
                         b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2083', style=style),
                         c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2083', style=style));
         display(Markdown("The environments you specify will remain unchanged over the course of the three trials."))
         display(Markdown("If you are working with other students, then coordinate amongst yourselves so to ensure you explore as many variations of environments as possible."))
         display(Markdown("**Modelling the Population**"))         
         display(Markdown("In this section, you will define the characteristics of an initial population of stick ungulates that will inhabit each of your three environments. The height of the stick ungulates will be determined by two alleles and . Using the appropriate sliders below specify the heights that are determined by each phenotype. Next the appropriate sliders below specify the initial frequency of each allele. Once you have set the heights that are determined by each phenotype and the initial frequency of each allele a table of detailing your population will be outputted. The table will identify the three genotypes in your population, their frequency, corresponding heights and fitness in each of the three environments. Population averages will also be reported."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));

         display(Markdown("### <span style='color:blue'>Environment 3"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("**Modelling Natural Selection**"))
         display(Markdown("In this section you will track the evolution of your stick ungulates as they evolve in each environment over 10 generations. The code to do so has already been set up and all you have to do is run it.  For the first generation, the code will output the various calculations so that you can follow the steps of our model.  For generations 2 to 10, the code will simply output a table that reports the frequency of each allele, the frequency of each phenotype, the average fitness of the population and the average height of the population."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 3"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("**Trial 3.**"))
         display(Markdown("**Modelling the Environment**"))
         display(Markdown("In this section, you will define 3 different environments by specifying the resource values of a, b and c using the sliders below. Once you set the slider for the resource values of a, b and c the corresponding fitness function for that environment will be plotted immediately below. Use the graph of the fitness function to refine your selection of the resource values of a, b and c for each environment. In defining the three environments try to create three very different environments as reflected by the shape of the fitness function.  For example, you may wish to study a case in which there is selection for an extreme height, a case where there is selection for an intermediate height and a case where there is selection against an intermediate height."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         interact_manual(fitness_function,
                         a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2081', style=style),
                         b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2081', style=style),
                         c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2081', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         interact_manual(fitness_function,
                    a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2082', style=style),
                    b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2082', style=style),
                    c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2082', style=style));
         display(Markdown("### <span style='color:blue'>Environment 3"))
         interact_manual(fitness_function,
                         a=widgets.IntSlider(min=0, max=10, step=1, value=1, description='a\u2083', style=style),
                         b=widgets.IntSlider(min=0, max=10, step=1, value=1, description='b\u2083', style=style),
                         c=widgets.IntSlider(min=0, max=10, step=1, value=1, description='c\u2083', style=style));
         display(Markdown("The environments you specify will remain unchanged over the course of the three trials."))
         display(Markdown("If you are working with other students, then coordinate amongst yourselves so to ensure you explore as many variations of environments as possible."))
         display(Markdown("**Modelling the Population**"))         
         display(Markdown("In this section, you will define the characteristics of an initial population of stick ungulates that will inhabit each of your three environments. The height of the stick ungulates will be determined by two alleles and . Using the appropriate sliders below specify the heights that are determined by each phenotype. Next the appropriate sliders below specify the initial frequency of each allele. Once you have set the heights that are determined by each phenotype and the initial frequency of each allele a table of detailing your population will be outputted. The table will identify the three genotypes in your population, their frequency, corresponding heights and fitness in each of the three environments. Population averages will also be reported."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));

         display(Markdown("### <span style='color:blue'>Environment 3"))
         #display(Markdown("<span style='color:blue'>Phenotypes"))
         interact_manual(generation,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("**Modelling Natural Selection**"))
         display(Markdown("In this section you will track the evolution of your stick ungulates as they evolve in each environment over 10 generations. The code to do so has already been set up and all you have to do is run it.  For the first generation, the code will output the various calculations so that you can follow the steps of our model.  For generations 2 to 10, the code will simply output a table that reports the frequency of each allele, the frequency of each phenotype, the average fitness of the population and the average height of the population."))
         display(Markdown("### <span style='color:blue'>Environment 1"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 2"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));
         display(Markdown("### <span style='color:blue'>Environment 3"))
         interact_manual(natural_selection,
                         A=widgets.FloatSlider(min=0, max=1, step=0.1, value=0.1, description='Initial allele frequencies: A', style=style),
                         a1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: AA - Height', style=style),
                         b1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: Aa - Height', style=style),
                         c1=widgets.IntSlider(min=0, max=10, step=1, value=1, description='Phenotypes: aa - Height', style=style));


def natural_selection(A, a1, b1, c1):
    from IPython.display import HTML, display
    freq_AA = 'p\u00b2 = ' + '{0:.3g}'.format(A/10 * A/10)
    freq_Aa = 'p(1-p) = ' + '{0:.3g}'.format(A/10 * (1 - A)/10)
    freq_aa = '(1 - p)\u00b2 = ' + '{0:.3g}'.format((1-A)/10 * (1-A)/10)
    x_AA = A/10 * A/10
    x_Aa = A/10 * (1 - A)/10
    x_aa = (1-A)/10 * (1-A)/10
    fitness_AA = 'w\u2081' + '(xAA) = ' + '{0:.3g}'.format(x_AA * x_AA * a1 + b1 * x_AA + c1)
    fitness_Aa = 'w\u2081' + '(xAa) = ' + '{0:.3g}'.format(x_Aa * x_Aa * a1 + b1 * x_Aa + c1)
    fitness_aa = 'w\u2081' + '(xaa) = ' + '{0:.3g}'.format(x_aa * x_aa * a1 + b1 * x_aa + c1)
    wAA = x_AA * x_AA * a1 + b1 * x_AA + c1
    wAa = x_Aa * x_Aa * a1 + b1 * x_Aa + c1
    waa = x_aa * x_aa * a1 + b1 * x_aa + c1
    display(Markdown("In the first generation, we have the frequency and fitness of each phenotype:"))
    data_generation_1 = [['Phenotype', 'Frequency', 'Fitness'],
            ['AA', freq_AA, fitness_AA],
            ['Aa', freq_Aa, fitness_Aa],
            ['aa', freq_aa, fitness_aa],
            ]
    display(HTML(
        '<table><tr>{}</tr></table>'.format(
            '</tr><tr>'.join(
                '<td>{}</td>'.format('</td><td>'.join(str(_) for _ in row)) for row in data_generation_1)
        )
    ))
    display(Markdown("The average fitness of the population then is the fitness of each phonotype weighted by their frequency."))
    w1_calculation = A*A*wAA + 2*A*(1-A)*wAa + (1-A)*(1-A)*waa  
    w1 = 'W\u2081 = ' + 'p\u00b2' + 'w\u2081' + 'x(AA) + 2p(1-p)' + 'w\u2081' + 'x(Aa) + ' + '(1-p)\u00b2' + 'w\u2081(x(aa)) = ' + '{0:.3g}'.format(w1_calculation)
    print(w1)
    display(Markdown("The average fitness of allele A is:"))
    wA_calculation = A*A*wAA + A*(1-A)*wAa
    wA = 'W(A) = ' + '(' + 'p\u00b2' + 'w(AA) + p(1-p)w(Aa)) = ' + '{0:.3g}'.format(wA_calculation)
    print(wA)
    display(Markdown("The frequency of A in next generation is:"))
    p = "p' = W(A)/W = " + "((1 - p)\u00b2" + "w(aa) + " + "p(1-p)w(Aa))/W = " + '{0:.3g}'.format(wA_calculation/w1_calculation)
    print(p)
    
    freqA_2 = wA_calculation/w1_calculation
    w2_calculation = cal_fitness_function(wAA, wAa,waa,wA_calculation/w1_calculation)
    wA2_calculation = cal_A_fitness(wA_calculation/w1_calculation,wAA, wAa)

    freqA_3 = wA2_calculation/w2_calculation
    w3_calculation = cal_fitness_function(wAA,wAa, waa,wA2_calculation/w2_calculation)
    wA3_calculation = cal_A_fitness(wA2_calculation/w2_calculation, wAA, wAa)

    freqA_4 = wA3_calculation/w3_calculation
    w4_calculation = cal_fitness_function(wAA,wAa,waa,wA3_calculation/w3_calculation)
    wA4_calculation = cal_A_fitness(wA3_calculation/w3_calculation, wAA, wAa)

    freqA_5 = wA4_calculation/w4_calculation
    w5_calculation = cal_fitness_function(wAA, wAA,waa,wA4_calculation/w4_calculation)
    wA5_calculation = cal_A_fitness(wA4_calculation/w4_calculation, wAA, wAa)

    freqA_6 = wA5_calculation / w5_calculation
    w6_calculation = cal_fitness_function(wAA, wAa, waa, wA5_calculation / w5_calculation)
    wA6_calculation = cal_A_fitness(wA5_calculation / w5_calculation, wAA, wAa)

    freqA_7 = wA6_calculation / w6_calculation
    w7_calculation = cal_fitness_function(wAA, wAa, waa, wA6_calculation / w6_calculation)
    wA7_calculation = cal_A_fitness(wA6_calculation / w6_calculation, wAA, wAa)

    freqA_8 = wA7_calculation / w7_calculation
    w8_calculation = cal_fitness_function(wAA, wAA, waa, wA7_calculation / w7_calculation)
    wA8_calculation = cal_A_fitness(wA7_calculation / w7_calculation, wAA, wAa)

    freqA_9 = wA8_calculation / w8_calculation
    w9_calculation = cal_fitness_function(wAA, wAa, waa, wA8_calculation / w8_calculation)
    wA9_calculation = cal_A_fitness(wA8_calculation / w8_calculation, wAA, wAa)

    freqA_10 = wA9_calculation / w9_calculation
    w10_calculation = cal_fitness_function(wAA, wAA, waa, wA9_calculation / w9_calculation)
    wA10_calculation = cal_A_fitness(wA9_calculation / w9_calculation, wAA, wAa)

    data = [['Generation', 'Frequency (A)', 'Frequency (a)', 'Frequency (AA)', 'Frequency (Aa)', 'Frequency (aa)', 'Fitness'],
            ['1', '{0:.3g}'.format(a1/10), '{0:.3g}'.format(1 - a1/10), '{0:.3g}'.format(a1/10 * a1/10), '{0:.3g}'.format(b1/10 * (1 - b1/10)), '{0:.3g}'.format((1-c1/10) * (1-c1/10)), '{0:.3g}'.format(w1_calculation)],
            ['2', '{0:.3g}'.format(freqA_2), '{0:.3g}'.format(1-freqA_2), '{0:.3g}'.format(freqA_2*freqA_2), '{0:.3g}'.format(freqA_2*(1-freqA_2)), '{0:.3g}'.format((1-freqA_2)*(1-freqA_2)), '{0:.3g}'.format(w2_calculation)],
            ['3', '{0:.3g}'.format(freqA_3), '{0:.3g}'.format(1-freqA_3), '{0:.3g}'.format(freqA_3*freqA_3), '{0:.3g}'.format(freqA_3*(1-freqA_3)), '{0:.3g}'.format((1-freqA_3)*(1-freqA_3)), '{0:.3g}'.format(w3_calculation)],
            ['4', '{0:.3g}'.format(freqA_4), '{0:.3g}'.format(1-freqA_4), '{0:.3g}'.format(freqA_4*freqA_4), '{0:.3g}'.format(freqA_4*(1-freqA_4)), '{0:.3g}'.format((1-freqA_4)*(1-freqA_4)), '{0:.3g}'.format(w4_calculation)],
            ['5', '{0:.3g}'.format(freqA_5), '{0:.3g}'.format(1-freqA_5), '{0:.3g}'.format(freqA_5*freqA_5), '{0:.3g}'.format(freqA_5*(1-freqA_5)), '{0:.3g}'.format((1-freqA_5)*(1-freqA_5)), '{0:.3g}'.format(w5_calculation)],
            ['6', '{0:.3g}'.format(freqA_6), '{0:.3g}'.format(1-freqA_6), '{0:.3g}'.format(freqA_6*freqA_6), '{0:.3g}'.format(freqA_6*(1-freqA_6)), '{0:.3g}'.format((1-freqA_6)*(1-freqA_6)), '{0:.3g}'.format(w6_calculation)],
            ['7', '{0:.3g}'.format(freqA_7), '{0:.3g}'.format(1-freqA_7), '{0:.3g}'.format(freqA_7*freqA_7), '{0:.3g}'.format(freqA_7*(1-freqA_7)), '{0:.3g}'.format((1-freqA_7)*(1-freqA_7)), '{0:.3g}'.format(w7_calculation)],
            ['8','{0:.3g}'.format(freqA_8), '{0:.3g}'.format(1-freqA_8), '{0:.3g}'.format(freqA_8*freqA_8), '{0:.3g}'.format(freqA_8*(1-freqA_8)), '{0:.3g}'.format((1-freqA_8)*(1-freqA_8)), '{0:.3g}'.format(w8_calculation)],
            ['9', '{0:.3g}'.format(freqA_9), '{0:.3g}'.format(1-freqA_9), '{0:.3g}'.format(freqA_9*freqA_9), '{0:.3g}'.format(freqA_9*(1-freqA_9)), '{0:.3g}'.format((1-freqA_9)*(1-freqA_9)), '{0:.3g}'.format(w9_calculation)],
            ['10', '{0:.3g}'.format(freqA_10), '{0:.3g}'.format(1-freqA_10), '{0:.3g}'.format(freqA_10*freqA_10), '{0:.3g}'.format(freqA_10*(1-freqA_10)), '{0:.3g}'.format((1-freqA_10)*(1-freqA_10)), '{0:.3g}'.format(w10_calculation)],
            ]
    display(Markdown("#### <span style='color:blue'>Table"))
    display(HTML(
        '<table><tr>{}</tr></table>'.format(
            '</tr><tr>'.join(
                '<td>{}</td>'.format('</td><td>'.join(str(_) for _ in row)) for row in data)
        )
    ))

def cal_fitness_function(wAA, wAa, waa, p):
    return p*p*wAA + p*(1-p)*wAa + (1-p)*(1-p)*waa
def cal_A_fitness(p,wAA, wAa):
    return p*p*wAA + p*(1-p)*wAa

def generation(A, a1, b1, c1):
    from IPython.display import HTML, display
    freq_AA = 'p\u00b2 = ' + '{0:.3g}'.format(A/10 * A/10)
    freq_Aa = 'p(1-p) = ' + '{0:.3g}'.format(A/10 * (1 - A)/10)
    freq_aa = '(1 - p)\u00b2 = ' + '{0:.3g}'.format((1-A)/10 * (1-A)/10)
    x_AA = A/10 * A/10
    x_Aa = A/10 * (1 - A)/10
    x_aa = (1-A)/10 * (1-A)/10
    fitness_AA = 'w\u2081' + '(xAA) = ' + '{0:.3g}'.format(a1 * x_AA * x_AA + b1* x_AA + c1)
    fitness_Aa = 'w\u2081' + '(xAa) = ' + '{0:.3g}'.format(a1 * x_Aa * x_Aa + b1* x_Aa + c1)
    fitness_aa = 'w\u2081' + '(xaa) = ' + '{0:.3g}'.format(a1 * x_aa * x_aa + b1* x_aa + c1)
    data = [['Phenotype', 'Frequency', 'Fitness'],
            ['AA', freq_AA, fitness_AA],
            ['Aa', freq_Aa, fitness_Aa],
            ['aa', freq_aa, fitness_aa],
            ]
    display(Markdown("#### <span style='color:blue'>Table"))
    display(HTML(
        '<table><tr>{}</tr></table>'.format(
            '</tr><tr>'.join(
                '<td>{}</td>'.format('</td><td>'.join(str(_) for _ in row)) for row in data)
        )
    ))
    display(Markdown("#### <span style='color:blue'>Phenotypic Distribution"))
    genes = ['AA', 'Aa', 'aa']
    x_pos = np.arange(len(genes))
    freqs = [a1/10 * a1/10, b1/10 * (1 - b1/10), (1-c1/10) * (1-c1/10)]
    # Build the plot
    fig, ax = plt.subplots()

    ax.bar(x_pos, freqs, align='center', alpha=0.5)
    ax.set_ylabel('Frequency')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(genes)
    ax.set_xlabel('Height')
    ax.yaxis.grid(True)
    # Save the figure and show
    plt.tight_layout()
    plt.show()
    display(Markdown("If you are working with other students, then coordinate amongst yourselves so to ensure you explore as many types of initial populations as possible. "))

def initialize_map(sparseness_of_map, size_of_map, number_of_groups):
    the_map = np.zeros((size_of_map, size_of_map))    
    for i in range(0, size_of_map):
        for j in range(0, i):
            group_i = int(i/(size_of_map/number_of_groups))
            group_j = int(j/(size_of_map/number_of_groups))            
            if random.random() > sparseness_of_map and abs(group_i - group_j) <= 1:
                the_map[i][j] = random.random()
                the_map[j][i] = the_map[i][j]
    ax = sns.heatmap(the_map)
    plt.show()
    #return the_map
def solve_map(sparseness_of_map, size_of_map, population_size, number_of_iterations, number_of_couples, number_of_winners_to_keep, mutation_probability,number_of_groups):
    
    # parameters
    #sparseness_of_map = 0.95
    #size_of_map = 1000
    #population_size = 30
    #number_of_iterations = 1000
    #number_of_couples = 9
    #number_of_winners_to_keep = 2
    #mutation_probability = 0.05
    #number_of_groups = 1
    
    # initialize the map and save it
    the_map = initialize_complex_map(sparseness_of_map, size_of_map, number_of_groups)

    # create the starting population
    population = create_starting_population(population_size, the_map)

    last_distance = 10000
    # for a large number of iterations do:
        
    for i in range(0,number_of_iterations):
        new_population = []
        
        # evaluate the fitness of the current population
        scores = score_population(population, the_map)

        best = population[np.argmin(scores)]
        number_of_moves = len(best)
        distance = fitness(best, the_map)
        if i == number_of_iterations:
            print('This is the best solution for the stimulation. The challenge has been completed')  # the end
        if distance != last_distance:
            print('Iteration %i: Best so far is %i steps for a distance of %f' % (i, number_of_moves, distance))
            print('Looking for a better solution...')
            plot_best(the_map, best, i)
        
        # allow members of the population to breed based on their relative score; 
            # i.e., if their score is higher they're more likely to breed
        for j in range(0, number_of_couples):  
            new_1, new_2 = crossover(population[pick_mate(scores)], population[pick_mate(scores)])
            new_population = new_population + [new_1, new_2]
  
        # mutate
        for j in range(0, len(new_population)):
            new_population[j] = np.copy(mutate(new_population[j], 0.05, the_map))
            
        # keep members of previous generation
        new_population += [population[np.argmin(scores)]]
        for j in range(1, number_of_winners_to_keep):
            keeper = pick_mate(scores)            
            new_population += [population[keeper]]
            
        # add new random members
        while len(new_population) < population_size:
            new_population += [create_new_member(the_map)]
            
        #replace the old population with a real copy
        population = copy.deepcopy(new_population)
                
        last_distance = distance

    # plot the results
def plot_best(the_map, route, iteration_number):
    ax = sns.heatmap(the_map)

    x=[0.5] + [x + 0.5 for x in route[0:len(route)-1]] + [len(the_map) - 0.5]
    y=[0.5] + [x + 0.5 for x in route[1:len(route)]] + [len(the_map) - 0.5]
    
    plt.plot(x, y, marker = 'o', linewidth=4, markersize=12, linestyle = "-", color='white')
    plt.show()

def print_pop(population):
    for i in population:
        print(i)


# Let's make a more complicated map that has at least 10 stops that have to be made and see what happens.

def initialize_complex_map(sparseness_of_map, size_of_map, number_of_groups):

    the_map = np.zeros((size_of_map, size_of_map))    
    for i in range(0, size_of_map):
        for j in range(0, i):
            group_i = int(i/(size_of_map/number_of_groups))
            group_j = int(j/(size_of_map/number_of_groups))            
            if random.random() > sparseness_of_map and abs(group_i - group_j) <= 1:
                the_map[i][j] = random.random()
                the_map[j][i] = the_map[i][j]
    ax = sns.heatmap(the_map)
    plt.show()
    return the_map

def create_starting_population(size, the_map):
    
    #this just creates a population of different routes of a fixed size. 
    
    population = []
    
    for i in range(0,size):
        population.append(create_new_member(the_map))
        
    return population

def fitness(route, the_map):
    
    score = 0
    
    for i in range(1, len(route)):
        if (the_map[route[i-1]][route[i]] == 0) and i != len(the_map)-1:
            print("WARNING: INVALID ROUTE")
            print(route)
            print(the_map)
        score = score + the_map[route[i-1]][route[i]]

    return score

def crossover(a, b):
        
    # I initially made an error here by allowing routes to crossover at any point, which obviously won't work
    # you have to insure that when the two routes cross over that the resulting routes produce a valid route
    # which means that crossover points have to be at the same position value on the map
    
    common_elements = set(a) & set(b)
    
    if len(common_elements) == 2:
        return (a, b)
    else:
        common_elements.remove(0)
        common_elements.remove(max(a)) 
        value = random.sample(common_elements, 1)        
    
    cut_a = np.random.choice(np.where(np.isin(a, value))[0])
    cut_b = np.random.choice(np.where(np.isin(b, value))[0])
    
    new_a1 = copy.deepcopy(a[0:cut_a])
    new_a2 = copy.deepcopy(b[cut_b:])
    
    new_b1 = copy.deepcopy(b[0:cut_b])
    new_b2 = copy.deepcopy(a[cut_a:])
    
    new_a = np.append(new_a1, new_a2)
    new_b = np.append(new_b1, new_b2)
       
    return (new_a, new_b)


def mutate(route, probability, the_map):
    
    new_route = copy.deepcopy(route)
    
    for i in range(1, len(new_route)):
        if random.random() < probability:
            
            go = True

            while go:

                possible_values = np.nonzero(the_map[new_route[i-1]])
                proposed_value = random.randint(0,len(possible_values[0])-1)
                route = np.append(new_route, possible_values[0][proposed_value])

                if new_route[i] == len(the_map)-1:
                    go = False
                else:
                    i += 1
    
    return new_route


def create_new_member(the_map):
    # here we are going to create a new route
    # the new route can have any number of steps, so we'll select that randomly
    # the structure of the route will be a vector of integers where each value is the next step in the route
    # Everyone starts at 0, so the first value in the vector will indicate where to attempt to go next.
    # That is, if v_i = 4, then that would correspond to X_0,4 in the map that was created at initialization
    
    # N is the size of the map, so we need to make sure that 
    # we don't generate any values that exceed the size of the map

    N = len(the_map)
    
    route = np.zeros(1, dtype=int)

    go = True
    
    i = 1
    
    while go:
        
        possible_values = np.nonzero(the_map[route[i-1]])
        proposed_value = random.randint(0,len(possible_values[0])-1)
        route = np.append(route, possible_values[0][proposed_value])
                
        if route[i] == N-1:
            go = False
        else:
            i += 1
    
    return route

def score_population(population, the_map):
    
    scores = []
    
    for i in range(0, len(population)):
        scores += [fitness(population[i], the_map)]
        
    return scores


def pick_mate(scores):

    array = np.array(scores)
    temp = array.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(array))

    fitness = [len(ranks) - x for x in ranks]
    
    cum_scores = copy.deepcopy(fitness)
    
    for i in range(1,len(cum_scores)):
        cum_scores[i] = fitness[i] + cum_scores[i-1]
        
    probs = [x / cum_scores[-1] for x in cum_scores]
    
    rand = random.random()
    
    for i in range(0, len(probs)):
        if rand < probs[i]:
            
            return i

def fitness_function(a, b, c):
    fitness = "Fitness function: W\u2081 = "
    x_square = "x\u00b2"
    final_string = fitness + str(a) + x_square + " + " + str(b) + "x" + " + " + str(c)
    print(final_string)
    print("Graph of Fitness:")
    import numpy as np
    import matplotlib.pyplot as plt

    x = np.linspace(0, 10, num=100)
    fx = []
    for i in range(len(x)):
        fx.append(a*x[i] ** 2 - b * x[i] + c)

    plt.plot(x, fx)
    plt.grid()
    plt.axvline()
    plt.axhline()
    plt.show()


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



def choose_direction(direction, x, y, z, x1, y1, z1):
    if direction == 'West(W)':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Color the floor area
        p = Rectangle((0, 0, 0), x, y, fc='gray')
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=0, zdir="z")

        # Draw a rectangle on the x=0 'wall'
        p1 = Rectangle((x1, y1), z / 2, y / 2, fc='none', ec='g', lw=z1)
        ax.add_patch(p1)
        art3d.pathpatch_2d_to_3d(p1, z=0, zdir="x")

        ax.set_xlim(0, x)
        ax.set_ylim(0, y)
        ax.set_zlim(0, z)

        plt.show()
        window_area = x1 * y1
        wall_area = x * y
        print('Window to Wall Ratio (WWR) = window area : wall area = ', str(window_area/wall_area))
    elif direction == 'North(N)':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(60)

        # Color the floor area
        p = Rectangle((0, 0, 0), x, y, fc='gray')
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=0, zdir="z")

        # Draw a rectangle on the x=0 'wall'
        p1 = Rectangle((x1, y1), z / 2, y / 2, fc='none', ec='g', lw=z1)
        ax.add_patch(p1)
        art3d.pathpatch_2d_to_3d(p1, z=0, zdir="x")

        ax.set_xlim(0, x)
        ax.set_ylim(0, y)
        ax.set_zlim(0, z)

        plt.show()
        window_area = x1 * y1
        wall_area = x * y
        print('Window to Wall Ratio (WWR) = window area : wall area = ', str(window_area/wall_area))
        # Windows
    elif direction == 'East(E)':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Color the floor area
        p = Rectangle((0, 0, 0), x, y, fc='gray')
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=0, zdir="z")

        # Draw a rectangle on the x=1 'wall'
        p2 = Rectangle((x1, y1), x/2, z/2, fc='none', ec='g', lw=z1)
        ax.add_patch(p2)
        art3d.pathpatch_2d_to_3d(p2, z=4, zdir="y")

        ax.set_xlim(0, x)
        ax.set_ylim(0, y)
        ax.set_zlim(0, z)
        plt.show()
        window_area = x1 * y1
        wall_area = x * y
        print('Window to Wall Ratio (WWR) = window area : wall area = ', str(window_area/wall_area))
    elif direction == 'South(S)':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(400)

        # Color the floor area
        p = Rectangle((0, 0, 0), x, y, fc='gray')
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=0, zdir="z")

        # Draw a rectangle on the x=0 'wall'
        p1 = Rectangle((x1, y1), z / 2, y / 2, fc='none', ec='g', lw=z1)
        ax.add_patch(p1)
        art3d.pathpatch_2d_to_3d(p1, z=0, zdir="x")

        ax.set_xlim(0, x)
        ax.set_ylim(0, y)
        ax.set_zlim(0, z)

        plt.show()
        window_area = x1 * y1
        wall_area = x * y
        print('Window to Wall Ratio (WWR) = window area : wall area = ', str(window_area/wall_area))
	