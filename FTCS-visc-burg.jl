function visc_FTCS(u0, L, N, dt, v, t_final)
    #u0 is initial traffic density, L is length of road, N is # of grid, dt is time steps, v is viscosity coefficient
    
    #note that velocity is not used here; to use it, one must modify the advection and diffusion terms
    
    using Plots

    dx = L / N
    x = range(0, L, length=N+1)

    # FTCS for viscous Burgers with advection and diffusion
    u = copy(u0)
    plot(x, u, xlim=(0, L), ylim=(0, 2), xlabel="Position", ylabel="Density", title="Traffic Density at Time t = 0")

    for t in 0:dt:t_final
        u_next = zeros(N+1)
        for i in 2:N
            #diffusion term
            diff = (v / dx^2) * (u[i+1] - 2*u[i] + u[i-1])
            #advection term
            adv= (dt / (2*dx)) * (u[i+1]^2 - u[i-1]^2)
            u_next[i] = u[i] - adv + diff
        end
        u = copy(u_next)
        
        plot!(x, u)
        plot!(legend=false)
        display(plot!)
    end
end
