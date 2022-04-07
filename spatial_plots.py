# plots module -- consider making subplots, etc

import numpy as np
import matplotlib.pyplot as plt
import spatial_model as model

# exploring the spatial model through examples 


# default parameters...

# these determine our "control" scenario 
CL = 30  # closure length in years
f = 0.25 # fishing effort 
m = 1    # number of closures 
n = 10   # number of patches 
frac_nomove = 0.95 #fraction of fish that stay put 
poaching = 0 #percentage of fishing effort not displaced by regulation 
IC = False   #initial condition of patches; False for all low, True for all high 

# then we vary every other degree of freedom separately:

# --- Time Series --- #
'''
# varying fishing 
model.graph_sol(CL, f, m, n, frac_nomove, IC, poaching)
model.graph_sol(CL, 0.35, m, n, frac_nomove, IC, poaching)
model.graph_sol(CL, 0.55, m, n, frac_nomove, IC, poaching)


# varying closure length 
model.graph_sol(5, f, m, n, frac_nomove, IC, poaching)
model.graph_sol(10, f, m, n, frac_nomove, IC, poaching)
model.graph_sol(40, f, m, n, frac_nomove, IC, poaching)

# varying m
model.graph_sol(CL, f, 2, n, frac_nomove, IC, poaching)
model.graph_sol(CL, f, 5, n, frac_nomove, IC, poaching)
model.graph_sol(CL, f, 7, n, frac_nomove, IC, poaching)

# varying n
model.graph_sol(CL, f, m, 2, frac_nomove, IC, poaching)
model.graph_sol(CL, f, m, 3, frac_nomove, IC, poaching)
model.graph_sol(CL, f, m, 5, frac_nomove, IC, poaching)

# varying dispersal
model.graph_sol(CL, f, m, n, 0, IC, poaching)
model.graph_sol(CL, f, m, n, 0.5, IC, poaching)
model.graph_sol(CL, f, m, n, 1, IC, poaching)

# varying poaching
model.graph_sol(CL, f, m, n, frac_nomove, IC, 0.1)
model.graph_sol(CL, f, m, n, frac_nomove, IC, 0.2)
model.graph_sol(CL, f, m, n, frac_nomove, IC, 0.5)
'''
# ------------------- #


# hysteresis region for fishing effort (C vs f for high, low coral)

model.hysteresis_zone()

# try to make a bifurcation diagram?

# C vs m plots for a given n (eg what's the best frac_close?)
# low fishing...
# medium fishing... 
# high fishing...
# experiment with showing oscillation amplitude variance...

# C vs n plots for a given m 

# heatmap of coral vs closure length and m/n for many patches

# for some poaching
# same for another value of f...

# for each of these, overlay the tau = t_thresh / p plot

# might need to calculate new thresholds....

# heatmap of coral vs fishing and closure length

# now we try making fishing a function of density...

# should we try adding in stochastic disturbances to illustrate the variance problem? 

# should we make a feedbacks flowchart in addition to the table of models?

