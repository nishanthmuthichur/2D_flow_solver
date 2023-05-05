import numpy as np

import cdiff_scheme_lib as cds

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
    ay = 1

    [Nx, Ny] = Flow_vec[U_VEL].U_sol.shape
    
    dudx = np.zeros((Nx, Ny))
    dvdx = np.zeros((Nx, Ny))    
    
    dudy = np.zeros((Nx, Ny))
    dvdy = np.zeros((Nx, Ny))    
    
    for y_idx in range(0, Ny):
        
        u = Flow_vec[U_VEL].U_sol[:, y_idx]
        v = Flow_vec[V_VEL].U_sol[:, y_idx]        

        dudx[:, y_idx] = -ax * cds.comp_CD8_deriv(u) / dx
        dvdx[:, y_idx] = -ax * cds.comp_CD8_deriv(v) / dx    
    
    for x_idx in range(0, Nx):
        
        u = Flow_vec[U_VEL].U_sol[x_idx, :]
        v = Flow_vec[V_VEL].U_sol[x_idx, :]        
        
        dudy[x_idx, :] = -ay * cds.comp_CD8_deriv(u) / dy
        dvdy[x_idx, :] = -ay * cds.comp_CD8_deriv(v) / dy            
    
    Flow_vec[U_VEL].F = dudx + dudy
    Flow_vec[V_VEL].F = dvdx + dvdy
    
    return Flow_vec


def compute_viscous_fluxes(dx, dy, \
                         Flow_vec): 


    return Flow_vec

