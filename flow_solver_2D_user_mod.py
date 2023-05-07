import numpy as np

import flow_solver_2D_lib as fs_2D

#*****************************

U_VEL = 0
V_VEL = 1

#*****************************

def set_init_cond(xcoord, ycoord, \
                               D):

    r0 = 0.5 * D;
    
    (Nx, Ny) = xcoord.shape
    
    u_vel_init = np.zeros((Nx, Ny))
    v_vel_init = np.zeros((Nx, Ny))

    Fu_init = np.zeros((Nx, Ny))
    Fv_init = np.zeros((Nx, Ny))    


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
        
    Flow_vec = fs_2D.flow_blk()
    
    Flow_vec.U_sol.append(u_vel_init)
    Flow_vec.U_sol.append(v_vel_init)

    Flow_vec.F_sol.append(Fu_init)
    Flow_vec.F_sol.append(Fv_init)

    print('usermod: Initial conditions has been computed')        
        
    return Flow_vec

def set_boundary_cond(xcoord, ycoord, \
                                   D, \
                            Flow_vec):

    r0 = 0.5 * D + 0.05 * D * np.sin(2 * np.pi * Flow_vec.time)
    
    (Nx, Ny) = xcoord.shape
    
    for y_idx in range(0, Ny):
    
        y_co = ycoord[0, y_idx]
        
        b = 12 
        
        if (y_co == 0): 
            
            u_vel_in = 1    
        
        else: 
            
            u_vel_in = 0.5 * (1 - np.tanh(b * ((abs(y_co)/r0) - (r0/abs(y_co)))))
            
        Flow_vec.U_sol[U_VEL][0, y_idx] = u_vel_in
        Flow_vec.U_sol[V_VEL][0, y_idx] = 0.0
            
    print(f'usermod: time_idx = {Flow_vec.time_idx}. Boundary conditions has been computed')            
    
    return Flow_vec
