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
         display(Markdown("<br> *Introductory blurb about architecture from architecture team. Use a graphic from Pulp team.*"))
         display(Markdown("#### <span style='color:blue'>Variables in architecture"))
         display(Markdown("Select building size:"))
         interact_manual(construct_building,
         x=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Length: X(m)',style=style),
         y=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Width: Y(m)',style=style),
         z=widgets.FloatSlider(min=1, max=20, step=1, value=1,description='Height: Z(m)',style=style));
         display(Markdown("Select building & window sizes:"))
         display(Markdown("<br> *(all windows are set to same, and only 0 or 1 per wall)*"))
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
            z1=widgets.FloatSlider(min=1, max=10, step=1, value=1, description='Thickness of the window: Z(m)', style=style));
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

    elif option == 'Biology':
         display(Markdown("### <span style='color:blue'>BIOLOGY FOCUS SELECTED"))

    elif option == 'Computer Science':
         display(Markdown("### <span style='color:blue'>COMPUTER SCIENCE FOCUS SELECTED"))
         display(Markdown("#### <span style='color:red'>Challenge: Design the routes in the map to minimize the number of steps and the distance we need to take from one point to another."))
         display(Markdown("#### <span style='color:blue'>1. Let's Initialize the Map!"))
         interact_manual(initialize_map,
                         sparseness_of_map=widgets.FloatSlider(min=0.05, max=0.99, step=0.01, value=0.99, description='Spareness: ', style=style),
                         size_of_map=widgets.IntSlider(min=1, max=20, step=1, value=20, description='Size: ', style=style),
                         number_of_groups=widgets.IntSlider(min=1, max=9, step=1, value=9, description='Number of Groups: ', style=style));
         #display(Markdown("#### <span style='color:blue'>2. Let's Put Everything Together to Design the Best Routes in the Map!"))
         display(Markdown("After initializing the map, to solve the challenge, we need to apply genetic algorithm which includes five phases:"))
         display(Markdown("* Initial Population: The process begins with a set of individuals which is called a Population. Each individual is a solution to the problem you want to solve. "))
         display(Markdown("* Fitness Funciton: The fitness function determines how fit an individual is (the ability of an individual to compete with other individuals). It gives a fitness score to each individual. The probability that an individual will be selected for reproduction is based on its fitness score."))
         display(Markdown("* Selection: The idea of selection phase is to select the fittest individuals and let them pass their genes to the next generation."))
         display(Markdown("* Crossover: Crossover is the most significant phase in a genetic algorithm. For each pair of parents to be mated, a crossover point is chosen at random from within the genes. Offspring are created by exchanging the genes of parents among themselves until the crossover point is reached."))
         display(Markdown("* Mutation: In certain new offspring formed, some of their genes can be subjected to a mutation with a low random probability. This implies that some of the bits in the bit string can be flipped."))
         display(Markdown("The algorithm terminates if the population has converged (does not produce offspring which are significantly different from the previous generation). Then it is said that the genetic algorithm has provided a set of solutions to our problem."))
         display(Markdown("#### <span style='color:blue'>2. Let's Put Everything Together to Design the Best Routes for the Map!"))
         interact_manual(solve_map,
                         sparseness_of_map=widgets.FloatSlider(min=0.05, max=0.99, step=0.01, value=0.99, description='Spareness: ', style=style),
                         size_of_map=widgets.IntSlider(min=1, max=20, step=1, value=20, description='Size: ', style=style),
                         number_of_groups=widgets.IntSlider(min=1, max=9, step=1, value=9, description='Number of Groups: ', style=style),
                         population_size=widgets.IntSlider(min=1, max=5, step=1, value=5, description='Population Size: ', style=style),
                         number_of_iterations=widgets.IntSlider(min=100, max=500, step=1, value=500, description='# Iterations: ', style=style),
                         number_of_couples=widgets.IntSlider(min=1, max=10, step=1, value=10, description='Size: ', style=style),
                         number_of_winners_to_keep=widgets.IntSlider(min=1, max=5, step=1, value=1, description='# Winners: ', style=style),
                         mutation_probability=widgets.FloatSlider(min=0.01, max=1, step=0.01, value=0.05, description='Mutation Probability: ', style=style));
    elif option == 'Math':
         display(Markdown("### <span style='color:blue'>MATH FOCUS SELECTED"))
         display(Markdown("In this Math Lab we will use mathematical modelling to explore the process of adaptation through Natural Selection.  Before completing this lab, students should be familiar with our previously developed model for Natural Selection acting at a single diploid gene locus with two alleles. We will use this model to consider how populations of stick ungulates evolve under different selective pressures, created by different environments.  "))
         display(Markdown("*What is a stick ungulate you ask?*"))
         display(HTML("""<a href="https://austaff-my.sharepoint.com/personal/jgreenwoodlee_athabascau_ca/Documents/Callysto/Ungulates" target="_blank">Ungulates</a> are large hooved mammals that come in both an even toed
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
                         a=widgets.IntSlider(min=0, max=20, step=1, value=1, description='a\u2081', style=style),
                         b=widgets.IntSlider(min=0, max=20, step=1, value=1, description='b\u2081', style=style),
                         c=widgets.IntSlider(min=0, max=20, step=1, value=1, description='c\u2081',
                                               style=style));


    elif option == 'Physics':
         display(Markdown("### <span style='color:blue'>PHYSICS FOCUS SELECTED"))

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

    last_distance = 1000000000
    # for a large number of iterations do:
        
    for i in range(0,number_of_iterations):
        new_population = []
        
        # evaluate the fitness of the current population
        scores = score_population(population, the_map)

        best = population[np.argmin(scores)]
        number_of_moves = len(best)
        distance = fitness(best, the_map)
        if distance != last_distance:
            print('Iteration %i: Best so far is %i steps for a distance of %f' % (i, number_of_moves, distance))
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
        plot_best(the_map, route, iteration_number)
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
	