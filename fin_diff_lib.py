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

def compute_CD10_filter(U):
    
    N_pts = len(U);
    fil_U = np.zeros(N_pts);
  
    m = 5;
  
    A = np.array([ \
        [ 0,   0,   0,   0,    0,   1,   -5,  10, -10,  5, -1], \
        [ 0,   0,   0,   0,   -5,  26,  -55,  60, -35, 10, -1], \
        [ 0,   0,   0,  10,  -55, 126, -155, 110, -45, 10, -1], \
        [ 0,   0, -10,  60, -155, 226, -205, 120, -45, 10, -1], \
        [ 0,   5, -35, 110, -205, 251, -210, 120, -45, 10, -1], \
        [-1,  10, -45, 120, -210, 252, -210, 120, -45, 10, -1], \
        [-1,  10, -45, 120, -210, 251, -205, 110, -35,  5,  0], \
        [-1,  10, -45, 120, -205, 226, -155,  60, -10,  0,  0], \
        [-1,  10, -45, 110, -155, 126,  -55,  10,   0,  0,  0], \
        [-1,  10, -35,  60,  -55,  26,   -5,   0,   0,  0,  0], \
        [-1,   5, -10,  10,   -5,   1,    0,   0,   0,  0,  0]    
    ])
  
    for idx in range(0, N_pts):
        
        if ( (idx >= 5) and (idx <= (N_pts - 6)) ):
        
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[5, (m - 5) : (m + 5)] * U[idx - 5 : idx + 5])
        
        elif (idx == 0):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m    ) : (m + 5)] * U[idx     : idx + 5])
            
        elif (idx == 1):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 1) : (m + 5)] * U[idx - 1 : idx + 5])
            
        elif (idx == 2):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 2) : (m + 5)] * U[idx - 2 : idx + 5])            
            
        elif (idx == 3):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 3) : (m + 5)] * U[idx - 3 : idx + 5])            
            
        elif (idx == 4):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 4) : (m + 5)] * U[idx - 4 : idx + 5])            
            
        elif (idx == (N_pts - 5)):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 5) : (m + 4)] * U[idx - 5 : idx + 4])                        
            
        elif (idx == (N_pts - 4)):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 5) : (m + 3)] * U[idx - 5 : idx + 3])            
            
        elif (idx == (N_pts - 3)):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 5) : (m + 2)] * U[idx - 5 : idx + 2])                        
            
        elif (idx == (N_pts - 2)):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 5) : (m + 1)] * U[idx - 5 : idx + 1])                                    
            
        elif (idx == (N_pts - 1)):
            
            fil_U[idx] = U[idx] - (1 / 1024) * sum(A[0, (m - 5) : (m    )] * U[idx - 5 : idx    ])            

        else: 
            
            print('fin_diff_lib: CD10_filter: Error! Input variable not entering any of the if..else conditions.')

    return fil_U



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

    