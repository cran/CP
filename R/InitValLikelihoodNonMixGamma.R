InitValLikelihoodNonMixGamma <- function(data) {

  ## Calculates initial values for the maximum likelihood estimation
  ## in the non-mixure model with Gamma type survival, i. e.
  ##   S(t) = c^Gamma^(0)(a, b * t), a > 0, b > 0, 0 < c < 1, t >= 0,
  ## via method of moments for a and b as the parameters
  ## in the simple Gamma model (taking into account only
  ## the patiens who died) and using the relative frequency of survival
  ## for the estimation of c.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns with the group
  ##         in the first (all of the same group),
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##
  ## Results:
  ##   Returns calculated initial values for the maximum likelihood estimation.

  # auxiliary variables
  t       <- data[, 3][data[, 2] == 1]
  n.d     <- length(x = t)
  m1      <- sum(t) / n.d
  m2      <- sum(t^2) / n.d
  n.alive <- sum(1 - data[, 2])
  n       <- length(x = data[, 1])

  # initial values...
  b.0 <- m1 / (m2 - m1^2)
  a.0 <- m1 * b.0
  c.0 <- n.alive / n

  # ...and if applicable projection into feasible region
  eps <- 1e-6
  if (a.0 < eps) {
    a.0 <- eps
  }
  if (b.0 < eps) {
    b.0 <- eps
  }
  if (c.0 < eps) {
    c.0 <- eps
  }
  if (c.0 > 1 - eps) {
    c.0 <- 1 - eps
  }

  return(c(a.0, b.0, c.0))
}