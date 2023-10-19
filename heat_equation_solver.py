import numpy as np

def tridiag(v, d, w, N):
    e = np.ones(N)
    A = v*np.diag(e[1:],-1)+d*np.diag(e)+w*np.diag(e[1:],1)
    return A

def heat_equation_solver(g0, g1, f, t_final, X, T, method):
    
    # g0 and g1 = *functions* at the boundary values at x = 0 and x = 1
    # f = *function* for initial condition
    # 1/X is space step size
    # 1/T is time step size
    # method inputs are 'FTCS', 'FTBS, 'CN' (for Crank-Nicholson)
    
    Dx = 1/X
    Dt = t_final/T
    
    x = np.linspace(0, 1, X+1)
    t = np.linspace(0, t_final, T+1)
    
    u = np.zeros((X+1, T+1))
    u[:,0] = f(x)
    r = Dt/(Dx**2)
    ingp = u[1:-1,0]
    A = tridiag(1, -2, 1, X-1)
    
    if method == 'FTCS':
        for j in range(T):
            u[1:-1, j+1] = u[1:-1, j] + r*(u[2:, j] - 2*u[1:-1,j] + u[0:-2, j])
            u[0,j+1] = g0(t[j+1])
            u[X,j+1] = g1(t[j+1])
    elif method == 'FTBS':
        K = np.eye(X-1) - r*A
        for j in range(T):
            v = np.copy(ingp) 
            v[0] = v[0] + r*g0(t[j+1])
            v[-1] = v[-1] + r*g1(t[j+1])
            ingp = np.linalg.solve(K, v)
            u[1:-1,j+1] = ingp
            u[0, j+1] = g0(t[j+1]) 
            u[X, j+1] = g1(t[j+1])
    elif method == 'CN':
        K = np.eye(X-1) - 0.5*r*A
        for j in range(T):
            b = np.dot(np.eye(X-1)+0.5*r*A, ingp)
            b[0] = b[0] + 0.5*r*(g0(t[j])+g0(t[j+1]))
            b[-1] = b[-1] + 0.5*r*(g1(t[j])+g1(t[j+1]))
            ingp = np.linalg.solve(K, v)
            u[1:-1,j+1] = ingp
            u[0, j+1] = g0(t[j+1]) 
            u[X, j+1] = g1(t[j+1])
    return u

    
