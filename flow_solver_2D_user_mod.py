import numpy as np

import incomp_NS_2D_lib as NS_2D

#*****************************

def two_D_grid_gen(Nx, Ny, \
                   x_min, x_max, \
                   y_min, y_max):
    
    xcoord = np.zeros((Nx, Ny))        
    ycoord = np.zeros((Nx, Ny))            
    
    dx = (x_max - x_min) / (Nx - 1)
    dy = (y_max - y_min) / (Ny - 1)    
    
    for x_idx in range(0, Nx):
      for y_idx in range(0, Ny):  
    
          xcoord[x_idx, y_idx] = x_min + dx * x_idx        
          ycoord[x_idx, y_idx] = y_min + dy * y_idx

    print('Incomp_NS_2D: 2D uniform grid generated')
    
    return xcoord, ycoord

#*****************************

def set_init_cond(xcoord, ycoord, \
                               D):

    r0 = 0.5 * D;
    
    (Nx, Ny) = xcoord.shape
    
    u_vel_init = np.zeros((Nx, Ny))
    v_vel_init = np.zeros((Nx, Ny))
    
    for x_idx in range(0, Nx):
      for y_idx in range(0, Ny):  
        
          y_co = ycoord[x_idx, y_idx]
        
          b = 12
        
          if (y_co == 0): 
            
              u_vel = 1    
        
          else: 
            
              u_vel = 0.5 * (1 - np.tanh(b * ((abs(y_co)/r0) - (r0/abs(y_co)))))
        
          u_vel_init[x_idx, y_idx] = u_vel
          v_vel_init[x_idx, y_idx] = 0.0
        
    
    U_sol_init = NS_2D.flow_sol(u_vel_init, \
                                v_vel_init)
    
    print('Incomp_NS_2D: Initial conditions has been computed')        
        
    return U_sol_init    

def set_boundary_cond(xcoord, ycoord, \
                                   D, \
                               U_sol):

    r0 = 0.5 * D
    
    (Nx, Ny) = xcoord.shape
    
    for y_idx in range(0, Ny):
    
        y_co = ycoord[0, y_idx]
        
        b = 12
        
        if (y_co == 0): 
            
            u_vel = 1    
        
        else: 
            
            u_vel = 0.5 * (1 - np.tanh(b * ((abs(y_co)/r0) - (r0/abs(y_co)))))
            
        U_sol.u_vel[0, y_idx] = u_vel
        U_sol.v_vel[0, y_idx] = 0.0    
    
    print('Incomp_NS_2D: Boundary conditions has been computed')            
    
    return U_sol
