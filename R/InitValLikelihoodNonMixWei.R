InitValLikelihoodNonMixWei <- function(data) {

  ## Calculates initial values for the maximum likelihood estimation
  ## in the non-mixure model with Weibull type survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t^k)], lambda > 0, k > 0, 0 < c < 1, t >= 0
  ## via least square method for lambda and k as the parameters
  ## in the simple Weibull model (taking into account only
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
  t               <- data[, 3][data[, 2] == 1]
  sort.t          <- sort(x = t)
  log.sort.t      <- log(sort.t)
  n.d             <- length(x = t)
  i               <- 1 : n.d
  F               <- (i - 0.3) / (n.d + 0.4)
  log.minus.log.S <- log(- log(1 - F))
  # parameters from least square method
  b               <- ((sum(log.sort.t * log.minus.log.S) - sum(log.sort.t) * sum(log.minus.log.S) / n.d)
                        / (sum(log.sort.t^2) - sum(log.sort.t)^2 / n.d))
  a               <- sum(log.minus.log.S) / n.d - b * sum(log.sort.t) / n.d
  # further auxiliary variables
  n.alive         <- sum(1 - data[, 2])
  n               <- length(x = data[, 1])

  # initial values...
  lambda.0 <- exp(a)
  k.0      <- b
  c.0      <- n.alive / n

  # ...and if applicable projection into feasible region
  eps <- 1e-6
  if (lambda.0 < eps) {
    lambda.0 <- eps
  }
  if (k.0 < eps) {
    k.0 <- eps
  }
  if (c.0 < eps) {
    c.0 <- eps
  }
  if (c.0 > 1 - eps) {
    c.0 <- 1 - eps
  }

  return(c(lambda.0, k.0, c.0))
}