import numpy as np
import matplotlib.pyplot as plt

def burgers_ftbs(f, X, T, a, t_final):
    # f = *function* for initial condition
    # 1/X is space step size
    # 1/T is time step size

    Dx = 1/X
    Dt = t_final/T

    x = np.linspace(0, 1, X+1)  
    t = np.linspace(0, t_final, T+1)
    u = np.zeros((T, X))  

    u[:,0] = f(x)

    # FTBS scheme for Burger's equation
    for i in range(1, T+1):
        for j in range(X+1):
            u[i, j] = u[i - 1, j] - (a * Dt / Dx) * u[i - 1, j] * (u[i - 1, j] - u[i - 1, j - 1])

    return u
