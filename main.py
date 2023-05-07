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

#*****************************INPUT PARAMETERS*********************************

wkdir_grid_read = 'C:\\Users\\Nishanth\\Desktop\\nish_work\\python_proj\\flow_solver_2D\\'
ip_fname = 'ip_grid_2D_101_x_101.h5'

wkdir_write_data = 'C:\\Users\\Nishanth\\Desktop\\nish_work\\python_proj\\flow_solver_2D\\op_data\\'
op_gen_fname = 'lconv_2D_results'

D = 1
U_MAX = 1

CFL = 0.1/U_MAX

N_tstep = 100

op_freq_idx = 10
op_idx = 0
#****************************START OF CODE*************************************

# Read grid data from file
xcoord, ycoord = fs_2D.read_grid_data(wkdir_grid_read, ip_fname)

(Nx, Ny) = xcoord.shape

x_max = xcoord.max()
x_min = xcoord.min()

y_max = ycoord.max()
y_min = ycoord.min()

dx = (x_max - x_min)/(Nx - 1)
dy = (y_max - y_min)/(Ny - 1)

# Compute time_step from CFL
delta_t = fs_2D.comp_time_step(xcoord, \
                               ycoord, \
                                  CFL)
   
# Set initial conditions from usermodule        
Flow_vec = um.set_init_cond(xcoord, \
                            ycoord, \
                                 D)

compute_fluxes = lconv_2D.compute_fluxes    
    
for time_idx in range(0, N_tstep):

    time = delta_t * time_idx    

    Flow_vec.time_idx = time_idx
    Flow_vec.time     = time

    if (np.mod(Flow_vec.time_idx, op_freq_idx) == 0):

        fs_2D.write_data_to_hdf5_file(wkdir_write_data, op_gen_fname, \
                                                              op_idx, \
                                                            Flow_vec)
    
        op_idx = op_idx + 1


    Flow_vec = um.set_boundary_cond(xcoord, ycoord, \
                                                 D, \
                                          Flow_vec)

    Flow_vec = fdl.compute_RK4_time_step(compute_fluxes, Flow_vec, \
                                                           dx, dy, \
                                                          delta_t)

#*********************
    
print('flow_solver_2D: Executed successfully')
    
    






















