import numpy as np
import matplotlib.pyplot as plt

def burgers_lax_wendroff_scheme(initial_condition, nx, nt, a, L, T):
    dx = L / (nx - 1)  # Spatial step size
    dt = T / (nt - 1)  # Temporal step size

    x = np.linspace(0, L, nx)  # Spatial grid
    u = np.zeros((nt, nx))  # Array to store solution at each time step

    # Initial condition
    u[0] = initial_condition(x)

    for i in range(1, nt):
        u[i, 1:-1] = (u[i-1, 1:-1] - (a * dt / (2 * dx)) * (u[i-1, 2:] - u[i-1, :-2]) +
                      (a**2 * dt**2 / (2 * dx**2)) * (u[i-1, 2:] - 2 * u[i-1, 1:-1] + u[i-1, :-2]))

    return x, u
