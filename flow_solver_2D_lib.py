import numpy as np

#***************************MACRO DEFINITIONS**********************************

INIT_DUM_VAL = -100000

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
    







           
          
          
          