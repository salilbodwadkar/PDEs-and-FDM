function Lax_Friedrich(u0, c, dx, dt, t_final)
    #u0 is initial condition matrix, c is velocity, dx and dt are space and time steps repsectively

    n = length(u0)
    nt = Int(t_final / dt)

    u = zeros(n, nt+1)
    u[:, 1] = u0
    
    for i in 1:nt
        u[:, i+1] = 0.5 * (u[mod1.(1:n-1, n), i] + u[mod1.(2:n, n), i]) - 0.5 * c * dt / dx * (u[mod1.(2:n, n), i] - u[mod1.(1:n-1, n), i])
    end
    
    return u
end
