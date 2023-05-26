function visc_FTCS(u0, dx, dt, v, t_final)
    #u0 is initial condition matrix, c is velocity, dx and dt are space and time steps repsectively, v is viscosity coefficient
    n = length(u0)
    
    nt = Int(t_final / dt)

    u = zeros(n, nt+1)
    u[:, 1] = u0

    for i in 1:nt
        u[2:n-1, i+1] = u[2:n-1, i] - 0.5 * dt / dx * (u[3:n, i].^2 - u[1:n-2, i].^2) + v * dt / dx^2 * (u[3:n, i] - 2 * u[2:n-1, i] + u[1:n-2, i])
    end

    return u
end
