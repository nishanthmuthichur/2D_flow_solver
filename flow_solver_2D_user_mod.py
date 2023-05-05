import numpy as np

import flow_solver_2D_lib as fs_2D

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

    print('flow_solver_2D: 2D uniform grid generated')
    
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
        
    Flow_vec = list()
    
    Flow_vec.append(fs_2D.flow_blk(u_vel_init))
    Flow_vec.append(fs_2D.flow_blk(v_vel_init))

    print('usermod: Initial conditions has been computed')        
        
    return Flow_vec

def set_boundary_cond(xcoord, ycoord, \
                                   D, \
                            Flow_vec):

    r0 = 0.5 * D
    
    (Nx, Ny) = xcoord.shape
    
    for y_idx in range(0, Ny):
    
        y_co = ycoord[0, y_idx]
        
        b = 12
        
        if (y_co == 0): 
            
            u_vel_in = 1    
        
        else: 
            
            u_vel_in = 0.5 * (1 - np.tanh(b * ((abs(y_co)/r0) - (r0/abs(y_co)))))
            
        Flow_vec[0].U_sol[0, y_idx] = u_vel_in
        Flow_vec[1].U_sol[0, y_idx] = 0.0            
            
        
    print('usermod: Boundary conditions has been computed')            
    
    return Flow_vec
