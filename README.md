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
