FLOW_SOLVER_2D PROJECT ARCHITECTURE



           ===========================         
           :                         :
           : user_module             :  
           :   => set_init_cond      : 
           :   => set_boundary_cond  :     
           :                         :  
           :                         :
           ===========================
                        |
                        |
                        V 
           ===========================          ==========================       
           :                         :          :                        :
           : Generic_2D_solver       :          :  fin_diff_lib          :
           :   => set_init_cond      :<=========:    => CD8 diff         :
           :   => set_boundary_cond  :          :    => RK4              :
           :   => compute_fluxes     :          :    => CD10 filter      :
           :                         :          :                        :
           ===========================          ==========================
                        |  
                        |               
                        V   
           =======================================
           :                                     :
           : Flow_solver_2D_lib                  :  
           :   => read_grid_file                 : 
           :   => compute_time_step              : 
           :   => Time loop for Rk4              :
           :   => Write output data as HDF5 file : 
           :                                     :
           =======================================


FLOW SOLVER STATUS

    1. The flow solver contains many for loops. Need to see if they could be replaced with list comprehensions
    2. An interfaces module need to be implemented. This interface module would allow different features in the 
       solver is to be easily replaceable. 
    3. The structure of the solver needs to be analysed and modified if needed. 
    4. There appears to be a bug in the filtering section of the solver. Need to resolved. 


FLOW SOLVER FUTURE UPDATES

    0. Filtering problem needs to be resolved. 
    1. The flow solver should be computed in the generalized coordinates. 
    2. Sponge layers need to be added. 
    3. Adaptive time_stepping feature needs to be implemented.
    4. Viscous stresses need to be included. 
    5. Non-linearity needs to be included
    6. Pressure + mass conservation needs to be included.
