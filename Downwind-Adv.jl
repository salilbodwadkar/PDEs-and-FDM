function downwind_scheme(u0, c, dx, dt, t_final)
    # u0 is initial condition matrix, c is velocity, dx and dt are space and time steps repsectively
    n = length(u0)

    nt = Int(t_final / dt)
  
    u = zeros(n, nt+1)
    u[:, 1] = u0

    for i in 1:nt
        u[:, i+1] = u[:, i] - c * dt / dx * (u[:, i] - u[mod1.(1:n-1, n), i])
    end
    return u
end
