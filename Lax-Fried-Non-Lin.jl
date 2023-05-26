function Lax_Fried_NL(u0, dx, dt, t_final, f)
   #u0 is initial condition matrix, c is velocity, dx and dt are space and time steps repsectively, f is flux function

    n = length(u0)
  
    nt = Int(t_final / dt)

    u = zeros(n, nt+1)
    u[:, 1] = u0

    for i in 1:nt
        u[:, i+1] = 0.5 * (u[mod1.(1:n-1, n), i] + u[mod1.(2:n, n), i]) - 0.5 * dt / dx * (f(u[mod1.(2:n, n), i]) - f(u[mod1.(1:n-1, n), i]))
    end
    return u
end
