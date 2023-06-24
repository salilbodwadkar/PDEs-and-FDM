visc_FTCS <- function(u0, L, N, dt, v, t_final) {
  # u0 is initial traffic density, L is length of road, N is number of grid, dt is time steps, v is viscosity coefficient
  # note that velocity is not used here; to use it, one must modify the advection and diffusion terms
  
  library(ggplot2)
  
  dx <- L / N
  x <- seq(0, L, length.out = N + 1)
  
  # FTCS for viscous Burgers with advection and diffusion  
  for (t in seq(0, t_final, by = dt)) {
    u_next <- rep(0, N + 1)
    for (i in 2:N) {
      # diffusion term
      diff <- (v / dx^2) * (u[i + 1] - 2 * u[i] + u[i - 1])
      # advection term
      adv <- (dt / (2 * dx)) * (u[i + 1]^2 - u[i - 1]^2)
      u_next[i] <- u[i] - adv + diff
    }
    u <- u_next
    
    p <- p + geom_line(aes(x, u))
    p <- p + theme(legend.position = "none")
    print(p)
  }
}
