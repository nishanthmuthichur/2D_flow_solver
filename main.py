# SCRIPT DESCRIPTION
#
# 1. Read grid data
# 2. Setup initial conditions
# 3. Time_step
#       a. Set boundary conditions
#       b. Compute fluxes
#       c. Compute RK_4 time step
# 4. Write output as HDF5 file

# ut + ax ux + ay uy = 0
# vt + ax vx + ay vy = 0

import numpy as np

import flow_solver_2D_lib as fs_2D
import flow_solver_2D_user_mod as um

import fin_diff_lib as fdl

import linear_conv_2D_sol as lconv_2D


import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#*****************************INPUT PARAMETERS*********************************

Nx = 121
Ny = 101

x_min = 0
x_max = 5
    
y_min = -1                              
y_max =  1

dx = (x_max - x_min)/(Nx - 1)
dy = (y_max - y_min)/(Ny - 1)

D = 1
U_MAX = 1

CFL = 0.1/U_MAX

N_tstep = 1000

#****************************START OF CODE*************************************

xcoord, ycoord = um.two_D_grid_gen(Nx, Ny, \
                             x_min, x_max, \
                             y_min, y_max)

delta_t = fs_2D.comp_time_step(xcoord, \
                               ycoord, \
                                  CFL)
        
Flow_vec = um.set_init_cond(xcoord, \
                            ycoord, \
                                 D)

u_vel = np.zeros((Nx, Ny, int(N_tstep/10)))
    
compute_fluxes = lconv_2D.compute_fluxes    
    
for time_idx in range(0, N_tstep):

    time = delta_t * time_idx    

    Flow_vec.time_idx = time_idx
    Flow_vec.time     = time

    Flow_vec = um.set_boundary_cond(xcoord, ycoord, \
                                                 D, \
                                          Flow_vec)

    Flow_vec = fdl.compute_RK4_time_step(compute_fluxes, Flow_vec, \
                                                           dx, dy, \
                                                          delta_t)

    #iNS_2D.write_output_to_file(U_sol)
    
    if ((time_idx % 10) == 0):
    
        u_vel[:, :, int(time_idx/10)] = Flow_vec.U_sol[0]

#*********************

fig, ax = plt.subplots(1, 1)
fig.set_dpi(300)


pp = ax.pcolormesh(xcoord, ycoord, u_vel[:, :, 99], cmap = 'jet', vmin = 0, vmax = 1)
ax.axis('equal')
fig.colorbar(pp)

plt.show()
    
print('flow_solver_2D: Executed successfully')
    
    






















