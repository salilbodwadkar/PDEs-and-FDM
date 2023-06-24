heat_ftcs <- function(N, L, dt, t_final, alpha) {
  N <- 100  # N is of spatial discretization points, L is legnth of domain, dt is time step size
  
  dx <- L / (N - 1)  # Spatial step size

  x <- seq(0, L, length.out = N)  # Spatial grid

  T0 <- sin(pi * x)  # Initial condition
  T <- T0  # Solution array

  r <- alpha * dt / dx^2  

  t <- 0
  while (t < t_final) {
    T_new <- T + r * (c(T[-1], 0) - 2 * T + c(0, T[-N]))  # FTCS update
    T <- T_new
    t <- t + dt
  }

  return(T)
}
