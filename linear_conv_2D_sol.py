import numpy as np

import fin_diff_lib as fdl

#*****************************

U_VEL = 0
V_VEL = 1

#*****************************


def compute_fluxes(dx, dy, \
                 Flow_vec):
    
    Flow_vec = compute_inviscid_fluxes(dx, dy, \
                                     Flow_vec)
        
    Flow_vec = compute_viscous_fluxes(dx, dy, \
                                    Flow_vec)    

    return Flow_vec
    
    
def compute_inviscid_fluxes(dx, dy, \
                          Flow_vec): 
    
    ax = 1
    ay = 0

    [Nx, Ny] = Flow_vec.U_sol[U_VEL].shape
    
    dudx = np.zeros((Nx, Ny))
    dvdx = np.zeros((Nx, Ny))    
    
    dudy = np.zeros((Nx, Ny))
    dvdy = np.zeros((Nx, Ny))    
    
    for y_idx in range(0, Ny):
        
        u = Flow_vec.U_sol[U_VEL][:, y_idx]
        v = Flow_vec.U_sol[V_VEL][:, y_idx]
        
        dudx[:, y_idx] = -ax * fdl.compute_CD8_deriv(u) / dx
        dvdx[:, y_idx] = -ax * fdl.compute_CD8_deriv(v) / dx    
    
    for x_idx in range(0, Nx):
        
        u = Flow_vec.U_sol[U_VEL][x_idx, :]
        v = Flow_vec.U_sol[V_VEL][x_idx, :]
        
        dudy[x_idx, :] = -ay * fdl.compute_CD8_deriv(u) / dy
        dvdy[x_idx, :] = -ay * fdl.compute_CD8_deriv(v) / dy            
    
    Flow_vec.F_sol[U_VEL] = dudx + dudy
    Flow_vec.F_sol[V_VEL] = dvdx + dvdy
    
    return Flow_vec


def compute_viscous_fluxes(dx, dy, \
                         Flow_vec): 


    return Flow_vec

