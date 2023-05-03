import numpy as np

import cdiff_scheme_lib as cds

#***************************MACRO DEFINITIONS**********************************

INIT_DUM_VAL = -100000

#************************DATA STRUCTURE DEFINITIONS****************************

class flow_sol:
    
    def __init__(self, u_vel_init, \
                       v_vel_init):

        (Nx, Ny) = u_vel_init.shape

        self.u_vel = u_vel_init
        self.v_vel = v_vel_init
        self.Fu    = np.zeros((Nx, Ny))
        self.Fv    = np.zeros((Nx, Ny))    
    
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
    
    print(f'Incomp_NS_2D: CFL = {CFL}; time step = {dt} s')    
    
    return dt
    
def compute_fluxes(dx, dy, \
                    U_sol):
    
    
    U_sol = compute_inviscid_fluxes(dx, dy, \
                                     U_sol)
        
    U_sol = compute_viscous_fluxes(dx, dy, \
                                    U_sol)    

    return U_sol
    
    
def compute_inviscid_fluxes(dx, dy, \
                            U_sol): 
    
    ax = 1
    ay = 1

    [Nx, Ny] = U_sol.u_vel.shape    
    
    dudx = np.zeros((Nx, Ny))
    dvdx = np.zeros((Nx, Ny))    
    
    dudy = np.zeros((Nx, Ny))
    dvdy = np.zeros((Nx, Ny))    
    
    for y_idx in range(0, Ny):
        
        u = U_sol.u_vel[:, y_idx]
        v = U_sol.v_vel[:, y_idx]
        
        dudx[:, y_idx] = -ax * cds.comp_CD8_deriv(u) / dx
        dvdx[:, y_idx] = -ax * cds.comp_CD8_deriv(v) / dx    
    
    for x_idx in range(0, Nx):
        
        u = U_sol.u_vel[x_idx, :]
        v = U_sol.v_vel[x_idx, :]
        
        dudy[x_idx, :] = -ay * cds.comp_CD8_deriv(u) / dy
        dvdy[x_idx, :] = -ay * cds.comp_CD8_deriv(v) / dy            
    
    U_sol.Fu = dudx + dudy
    U_sol.Fv = dvdx + dvdy

    return U_sol
    


def compute_viscous_fluxes(dx, dy, \
                            U_sol): 


    return U_sol






           
          
          
          