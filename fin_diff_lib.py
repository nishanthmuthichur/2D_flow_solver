from operator import add

import numpy as np

# This function is used to compute the derivative 'dYdX' of 'Y' using 
# an explicit 8th order central difference stencil

def compute_CD8_deriv(Y):

    N_pts = len(Y);
    dYdX = np.zeros(N_pts);
  
    m = 4;
  
    A = np.array([ \
        [   0    ,    0     ,   0    ,    0   , (-25/12) ,  4    ,   -3    , (4/3)   ,  (-1/4)  ], \
        [   0    ,    0     ,   0    , (-1/4) ,  (-5/6)  , (3/2) , (-1/2)  , (1/12)  ,     0    ], \
        [   0    ,    0     , (1/12) , (-2/3) ,     0    , (2/3) , (-1/12) ,   0     ,     0    ], \
        [   0    , (-1/60)  , (3/20) , (-3/4) ,     0    , (3/4) , (-3/20) , (1/60)  ,     0    ], \
        [(1/280) , (-4/105) , (1/5)  , (-4/5) ,     0    , (4/5) , (-1/5)  , (4/105) , (-1/280) ], \
        [   0    , (-1/60)  , (3/20) , (-3/4) ,     0    , (3/4) , (-3/20) , (1/60)  ,     0    ], \
        [   0    ,    0     , (1/12) , (-2/3) ,     0    , (2/3) , (-1/12) ,   0     ,     0    ], \
        [   0    , (-1/12)  , (1/2)  , (-3/2) ,   (5/6)  , (1/4) ,    0    ,   0     ,     0    ], \
        [ (1/4)  , (-4/3)   ,   3    ,    -4  ,  (25/12) ,   0   ,    0    ,   0     ,     0    ], \
    ])
  
    for idx in range(0, N_pts):
        
        if (idx == 0):
                  
            dYdX[idx] =  A[0,    m   ] * Y[idx    ] + \
                         A[0, (m + 1)] * Y[idx + 1] + \
                         A[0, (m + 2)] * Y[idx + 2] + \
                         A[0, (m + 3)] * Y[idx + 3] + \
                         A[0, (m + 4)] * Y[idx + 4]          
            
        elif (idx == 1):
             
            dYdX[idx] =  A[1, (m - 1)] * Y[idx - 1] + \
                         A[1,    m   ] * Y[idx    ] + \
                         A[1, (m + 1)] * Y[idx + 1] + \
                         A[1, (m + 2)] * Y[idx + 2] + \
                         A[1, (m + 3)] * Y[idx + 3]
                        
        elif (idx == 2):

            dYdX[idx] =  A[2, (m - 2)] * Y[idx - 2] + \
                         A[2, (m - 1)] * Y[idx - 1] + \
                         A[2,    m   ] * Y[idx    ] + \
                         A[2, (m + 1)] * Y[idx + 1] + \
                         A[2, (m + 2)] * Y[idx + 2]
                        
        elif (idx == 3):

            dYdX[idx] =  A[3, (m - 3)] * Y[idx - 3] + \
                         A[3, (m - 2)] * Y[idx - 2] + \
                         A[3, (m - 1)] * Y[idx - 1] + \
                         A[3,    m   ] * Y[idx    ] + \
                         A[3, (m + 1)] * Y[idx + 1] + \
                         A[3, (m + 2)] * Y[idx + 2] + \
                         A[3, (m + 3)] * Y[idx + 3]

        elif (idx == (N_pts - 4)):                       

            dYdX[idx] = A[5, (m - 3)] * Y[idx - 3] + \
                        A[5, (m - 2)] * Y[idx - 2] + \
                        A[5, (m - 1)] * Y[idx - 1] + \
                        A[5,    m   ] * Y[idx    ] + \
                        A[5, (m + 1)] * Y[idx + 1] + \
                        A[5, (m + 2)] * Y[idx + 2] + \
                        A[5, (m + 3)] * Y[idx + 3]
                         
        elif (idx == (N_pts - 3)):
            
            dYdX[idx] = A[6, (m - 2)] * Y[idx - 2] + \
                        A[6, (m - 1)] * Y[idx - 1] + \
                        A[6,    m   ] * Y[idx    ] + \
                        A[6, (m + 1)] * Y[idx + 1] + \
                        A[6, (m + 2)] * Y[idx + 2]
        
        elif (idx == (N_pts - 2)):
            
            dYdX[idx] = A[7, (m - 3)] * Y[idx - 3] + \
                        A[7, (m - 2)] * Y[idx - 2] + \
                        A[7, (m - 1)] * Y[idx - 1] + \
                        A[7,    m   ] * Y[idx    ] + \
                        A[7, (m + 1)] * Y[idx + 1] 
                          
        elif (idx == (N_pts - 1)): 
            
            dYdX[idx] = A[8, (m - 4)] * Y[idx - 4] + \
                        A[8, (m - 3)] * Y[idx - 3] + \
                        A[8, (m - 2)] * Y[idx - 2] + \
                        A[8, (m - 1)] * Y[idx - 1] + \
                        A[8,    m   ] * Y[idx    ]
                          
        else:
            
            dYdX[idx] = A[4, (m - 4)] * Y[idx - 4] + \
                        A[4, (m - 3)] * Y[idx - 3] + \
                        A[4, (m - 2)] * Y[idx - 2] + \
                        A[4, (m - 1)] * Y[idx - 1] + \
                        A[4,    m   ] * Y[idx    ] + \
                        A[4, (m + 1)] * Y[idx + 1] + \
                        A[4, (m + 2)] * Y[idx + 2] + \
                        A[4, (m + 3)] * Y[idx + 3] + \
                        A[4, (m + 4)] * Y[idx + 4]            
            
    return dYdX

def compute_RK4_time_step(compute_fluxes, Flow_vec_0, \
                                              dx, dy, \
                                             delta_t):

    Flow_vec_up = Flow_vec_0    
    
    #S1
    Flow_vec_1 = Flow_vec_0     
    #K1
    Flow_vec_1 = compute_fluxes(dx, dy, \
                            Flow_vec_1)
    #RK_stage_1    
    Flow_vec_up.U_sol = list( map( add, \
                                   Flow_vec_up.U_sol, \
                                   [var_idx * (delta_t / 6) for var_idx in Flow_vec_1.F_sol]   ) )    
    #Flow_vec_up.U_sol = list( map( add, Flow_vec_up.U_sol, ((delta_t / 6) * Flow_vec_1.F_sol) ) )
    #Flow_vec_up.U_sol = Flow_vec_up.U_sol + delta_t * (Flow_vec_1.F_sol/6)        

        
    #S2    
    Flow_vec_1.U_sol = list( map(add, Flow_vec_0.U_sol, [var_idx * (delta_t * 0.5) for var_idx in Flow_vec_1.F_sol] ) )
    #Flow_vec_1.U_sol = Flow_vec_0.U_sol + ((delta_t * 0.5) * Flow_vec_1.F_sol)
    #K2
    Flow_vec_1 = compute_fluxes(dx, dy, \
                            Flow_vec_1)    
    #RK_stage_2    
    Flow_vec_up.U_sol = list( map( add, \
                                   Flow_vec_up.U_sol, \
                                   [var_idx * (delta_t / 3) for var_idx in Flow_vec_1.F_sol]   ) )        
    #Flow_vec_up.U_sol = list( map( add, Flow_vec_up.U_sol, ((delta_t / 3) * Flow_vec_1.F_sol) ) )
    #Flow_vec_up.U_sol = Flow_vec_up.U_sol + delta_t * (Flow_vec_1.F_sol/3)                
        
    
    #S3    
    Flow_vec_1.U_sol = list( map(add, Flow_vec_0.U_sol, [var_idx * (delta_t * 0.5) for var_idx in Flow_vec_1.F_sol] ) )    
    #Flow_vec_1.U_sol = Flow_vec_0.U_sol + ((delta_t * 0.5) * Flow_vec_1.F_sol)    
    #K3
    Flow_vec_1 = compute_fluxes(dx, dy, \
                            Flow_vec_1)        
        
    #RK_stage_3    
    Flow_vec_up.U_sol = list( map( add, \
                                   Flow_vec_up.U_sol, \
                                   [var_idx * (delta_t / 3) for var_idx in Flow_vec_1.F_sol]   ) )        
    #Flow_vec_up.U_sol = list( map( add, Flow_vec_up.U_sol, ((delta_t / 3) * Flow_vec_1.F_sol) ) )
    #Flow_vec_up.U_sol = Flow_vec_up.U_sol + delta_t * (Flow_vec_1.F_sol/3)            
        
        
    #S4    
    Flow_vec_1.U_sol = list( map(add, Flow_vec_0.U_sol, [var_idx * (delta_t) for var_idx in Flow_vec_1.F_sol] ) )    
    #Flow_vec_1.U_sol = Flow_vec_0.U_sol + (delta_t * Flow_vec_1.F_sol)        
    #K4
    Flow_vec_1 = compute_fluxes(dx, dy, \
                            Flow_vec_1)            
        
    #RK_stage_4    
    Flow_vec_up.U_sol = list( map( add, \
                                   Flow_vec_up.U_sol, \
                                   [var_idx * (delta_t / 6) for var_idx in Flow_vec_1.F_sol]   ) )        
    #Flow_vec_up.U_sol = list( map( add, Flow_vec_up.U_sol, ((delta_t / 6) * Flow_vec_1.F_sol) ) )        
    #Flow_vec_up.U_sol = Flow_vec_up.U_sol + delta_t * (Flow_vec_1.F_sol/6)            

    return Flow_vec_up

    