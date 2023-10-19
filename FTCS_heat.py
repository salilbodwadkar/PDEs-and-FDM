import numpy as np
import matplotlib.pyplot as plt

def FTCS_heat(g0, g1, f, t_final, X, T):
    
    # g0 and g1 = *functions* at the boundary values at x = 0 and x = 1
    # f = *function* for initial condition
    # 1/X is space step size
    # 1/T is time step size
    
    Dx = 1/X
    Dt = t_final/T
    
    x = np.linspace(0, 1, X+1)
    t = np.linspace(0, t_final, T+1)
    
    u = np.zeros((X+1, T+1))
    u[:,0] = f(x)
    
    r = Dt/(Dx**2)
    
    for j in range(T):
        u[1:-1, j+1] = u[1:-1, j] + r*(u[2:, j] - 2*u[1:-1,j] + u[0:-2, j])
        u[0,j+1] = g0(t[j+1])
        u[X,j+1] = g1(t[j+1])
        
    return u

    
