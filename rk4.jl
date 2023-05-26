function rk4(f, a, b, u0, N)
    # Solve ODE y'(t) = f(t, y(t)) using Runge-Kutta 4.
    h = (b - a) / N
    t = collect(a:h:b)
    w = zeros(N+1, length(u0))
    w[1, :] = u0
    for i in 1:N
        k1 = h * f(t[i], w[i, :])
        k2 = h * f(t[i] + h/2, w[i, :] + k1/2)
        k3 = h * f(t[i] + h/2, w[i, :] + k2/2)
        k4 = h * f(t[i] + h, w[i, :] + k3)
        w[i+1, :] = w[i, :] + (k1 + 2*k2 + 2*k3 + k4) / 6
    end
    return t, w
end
