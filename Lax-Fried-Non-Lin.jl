function lax_friedrich_scheme(u0, dx, dt, t_final, f)
  
    n = length(u0)
  
    nt = Int(t_final / dt)

    u = zeros(n, nt+1)
    u[:, 1] = u0

    for i in 1:nt
        u[:, i+1] = 0.5 * (u[mod1.(1:n-1, n), i] + u[mod1.(2:n, n), i]) - 0.5 * dt / dx * (f(u[mod1.(2:n, n), i]) - f(u[mod1.(1:n-1, n), i]))
    end
    return u
end
