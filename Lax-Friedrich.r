library(ggplot2)

Lax_Friedrich <- function(u0, c, L, N, dt, t_final) {
  #u0 is initial traffic density, c is advection velocity, L is road legnth, N is # of grid points, dt is time step
  
  dx <- L / N
  x <- seq(0, L, length.out = N + 1)
  
  # Lax-Friedrichs scheme
  u <- u0
  p <- ggplot() + xlim(0, L) + ylim(0, 1.2) + xlab("Position") + ylab("Density") + ggtitle("Traffic Density at Time t = 0")
  p <- p + geom_line(aes(x, u))
  print(p)
  
  for (t in seq(0, t_final, by = dt)) {
    u_next <- rep(0, N + 1)
    for (i in 2:N) {
      u_next[i] <- (u[i + 1] + u[i - 1]) / 2 - (c * dt / (2 * dx)) * (u[i + 1] - u[i - 1])
    }
    u <- u_next
    
    p <- p + geom_line(aes(x, u))
    p <- p + theme(legend.position = "none")
    print(p)
  }
}
