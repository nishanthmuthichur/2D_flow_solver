import numpy as np
import h5py as h5

import fin_diff_lib as fdl
#***************************MACRO DEFINITIONS**********************************

INIT_DUM_VAL = -100000

U_VEL = 0
V_VEL = 1

#************************DATA STRUCTURE DEFINITIONS****************************

class flow_blk:
    
    def __init__(self):

        self.U_sol = list()
        self.F_sol = list()
    
        self.time_idx = 0
        self.time     = 0

    def __str__(self):
        
        return f"time_idx = {self.time_idx}"    

#*************************FUNCTION DEFINITIONS*********************************

def read_grid_data(wkdir_grid_read, ip_fname):
    
    ip_fname_abs = wkdir_grid_read + ip_fname

    print(f'flow_solver_2D: Reading grid data from file:{ip_fname_abs}')

    file_dp = h5.File(ip_fname_abs, 'r')

    xcoord = file_dp.get('/xcoord')
    ycoord = file_dp.get('/ycoord')

    xcoord = np.array(xcoord)
    ycoord = np.array(ycoord)    

    file_dp.close()

    return xcoord, ycoord

#****************************
def comp_time_step(xcoord, \
                   ycoord, \
                      CFL):
        
    dx_min = INIT_DUM_VAL
    dy_min = INIT_DUM_VAL
    
    (Nx, Ny) = xcoord.shape
    
    d_xcoord = xcoord[1 : (Nx - 1), :] - \
               xcoord[0 : (Nx - 2), :]
               
    d_ycoord = ycoord[:, 1 : (Ny - 1)] - \
               ycoord[:, 0 : (Ny - 2)]
    
    dx_min = d_xcoord.min()
    dy_min = d_ycoord.min()    
    
    if (dx_min < dy_min):
        
        dt = CFL * dx_min
        
    else:
    
        dt = CFL * dy_min    
    
    print(f'flow_solver_2D: CFL = {CFL}; time step = {dt} s')    
    
    return dt
    
#******************************************************************************

def write_data_to_hdf5_file(wkdir_write_data, op_gen_fname, \
                                                    op_idx, \
                                                  Flow_vec):

    op_fname_abs = wkdir_write_data + \
                       op_gen_fname + \
                                '_' + \
                        str(op_idx) + \
                               '.h5'
                              
    file_dp = h5.File(op_fname_abs, 'w')

    file_dp.create_dataset('/time'     , data = Flow_vec.time)
    file_dp.create_dataset('/time_idx' , data = Flow_vec.time_idx)    
    file_dp.create_dataset('/u_vel'    , data = Flow_vec.U_sol[U_VEL])
    file_dp.create_dataset('/v_vel'    , data = Flow_vec.U_sol[V_VEL])

    file_dp.close()                                   

def filter_sol(Flow_vec):
    
    N_var = len(Flow_vec.U_sol)
    
    (Nx, Ny) = Flow_vec.U_sol[0].shape
    
    for var_idx in range(0, N_var):
        
        sol_field = Flow_vec.U_sol[var_idx]
        
        #Filtering along y-direction
        for x_idx in range(0, Nx):
            
            sol_field[x_idx, :] = fdl.compute_CD10_filter(sol_field[x_idx, :])    
            
        for y_idx in range(0, Ny):
            
            sol_field[:, y_idx] = fdl.compute_CD10_filter(sol_field[:, y_idx])
            
        Flow_vec.U_sol[var_idx] = sol_field    
    
    return Flow_vec

         