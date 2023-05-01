# SCRIPT DESCRIPTION
#
# 1. Generate grid
# 2. Setup initial conditions
# 3. Time_step
#       a. Set boundary conditions
#       b. Compute fluxes
#       c. Compute RK_4 time step

import numpy as np

import two_D_incomp_NS as NS_2D

#*****************************INPUT PARAMETERS*********************************

Nx = 11
Ny = 5

x_min = 0
x_max = 1

y_min = -1
y_max =  1

#****************************START OF CODE*************************************

xcoord, ycoord = NS_2D.two_D_grid_gen(Nx, Ny, \
                                x_min, x_max, \
                                y_min, y_max)



