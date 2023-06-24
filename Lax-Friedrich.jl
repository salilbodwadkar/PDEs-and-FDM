function Lax_Friedrich(u0, c, L, N, dt, t_final)
    #u0 is initial traffic density, c is advection velocity, L is road legnth, N is # of grid points, dt is time step
    using Plots
    
    dx = L/N
    x = range(0, L, length=N+1)

    #Lax-Friedrichs scheme
    for t in 0:dt:t_final
        u_next = zeros(N+1)
        for i in 2:N
            u_next[i] = (u[i+1] + u[i-1])/2 - (c*dt/(2*dx))*(u[i+1] - u[i-1])
        end
        u = copy(u_next)

        plot!(x, u)
        plot!(legend=false)
        display(plot!)
    end
end
