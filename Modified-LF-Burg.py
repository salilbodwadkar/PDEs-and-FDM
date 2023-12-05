import numpy as np
import matplotlib.pyplot as plt

def burgers_LF(f, X, T, a, t_final):
    # f = *function* for initial condition
    # 1/X is space step size
    # 1/T is time step size

    Dx = 1/X
    Dt = t_final/T

    x = np.linspace(0, 1, X+1)  
    t = np.linspace(0, t_final, T+1)
    u = np.zeros((T, X))  

    u[:,0] = f(x)

    # Modified LF scheme for Burger's equation
    for i in range(T+1):
        for j in range(1,X+1):
            u_t = (u[i + 1, j] - 0.5 * (u[i, j + 1] + u[i, j - 1])) / Dt
            u_x = 0.5 * ((u[i, j + 1])**2 - (u[i, j - 1])**2) / (2 * Dx)
            u[i + 1, j] = u[i, j] - Dt * a * u_x + Dt * u_t

    return u
