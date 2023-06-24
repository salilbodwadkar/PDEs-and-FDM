rk4 <- function(f, a, b, u0, N) {
  # Solve ODE y'(t) = f(t, y(t)) using Runge-Kutta 4.
  h <- (b - a) / N
  t <- seq(a, b, by = h)
  w <- matrix(0, nrow = N + 1, ncol = length(u0))
  w[1, ] <- u0
  for (i in 1:N) {
    k1 <- h * f(t[i], w[i, ])
    k2 <- h * f(t[i] + h/2, w[i, ] + k1/2)
    k3 <- h * f(t[i] + h/2, w[i, ] + k2/2)
    k4 <- h * f(t[i] + h, w[i, ] + k3)
    w[i + 1, ] <- w[i, ] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  return(list(t = t, w = w))
}
