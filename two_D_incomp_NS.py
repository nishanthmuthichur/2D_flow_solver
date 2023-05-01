import numpy as np

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
    
    return xcoord, ycoord

