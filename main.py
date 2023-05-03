# SCRIPT DESCRIPTION
#
# 1. Read grid data
# 2. Setup initial conditions
# 3. Time_step
#       a. Set boundary conditions
#       b. Compute fluxes
#       c. Compute RK_4 time step
# 4. Write output as HDF5 file

# u_t + ax ux + ay uy = 0
# v_t + ax vx + ay vy = 0

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import incomp_NS_2D_lib as iNS_2D
import incomp_NS_2D_user_mod as um

#*****************************INPUT PARAMETERS*********************************

Nx = 121
Ny = 51

x_min = 0
x_max = 12
    
y_min = -5                              
y_max =  5

dx = (x_max - x_min)/(Nx - 1)
dy = (y_max - y_min)/(Ny - 1)

D = 1
U_MAX = 1

CFL = 0.1/U_MAX

N_tstep = 10

#****************************START OF CODE*************************************

xcoord, ycoord = um.two_D_grid_gen(Nx, Ny, \
                             x_min, x_max, \
                             y_min, y_max)

dt = iNS_2D.comp_time_step(xcoord, \
                          ycoord, \
                             CFL)
        
U_sol = um.set_init_cond(xcoord, \
                         ycoord, \
                              D)

for t_idx in range(0, N_tstep):

    time = dt * t_idx    

    U_sol = um.set_boundary_cond(xcoord, ycoord, \
                                              D, \
                                          U_sol)

    U_sol = iNS_2D.compute_fluxes(dx, dy, \
                                   U_sol)

    #U_sol = iNS_2D.comp_RK4_time_step(U_sol)    

    #iNS_2D.write_output_to_file(U_sol)
    
u_vel = U_sol.u_vel

#*********************

fig, ax = plt.subplots(1, 1)
fig.set_dpi(300)

ax.pcolormesh(xcoord, ycoord, u_vel, cmap = 'jet')
ax.axis('equal')


plt.show()

    
    
print('incomp_iNS_2D: Executed successfully')
    
    






















