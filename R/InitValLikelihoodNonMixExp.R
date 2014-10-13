InitValLikelihoodNonMixExp <- function(data) {

  ## Calculates initial values for the maximum likelihood estimation
  ## in the non-mixure model with exponential survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t)], lambda > 0, 0 < c < 1, t >= 0,
  ## via maximum likelihood estimation of lambda as the parameter
  ## in the simple exponential model (taking into account only
  ## the patients who died) and using the relative frequency of survival
  ## for the estimation of c.
  ##
  ## Args:
  ##   df: Data frame which consists of at least three columns
  ##       with the group in the first,
  ##       status (1 = event, 0 = censored) in the second
  ##       and event time in the third column.
  ##
  ## Results:
  ##   Returns calculated initial values for the maximum likelihood estimation.

  # auxiliary variables
  d       <- sum(data[, 2])
  o.d     <- sum(data[, 3][data[, 2] == 1])
  n.alive <- sum(1 - data[, 2])
  n       <- length(x = data[, 1])

  # initial values...
  lambda.0 <- d / o.d
  c.0      <- n.alive / n

  # ...and if applicable projection into feasible region
  eps <- 1e-6
  if (lambda.0 < eps) {
    lambda.0 <- eps
  }
  if (c.0 < eps) {
    c.0 <- eps
  }
  if (c.0 > 1 - eps) {
    c.0 <- 1 - eps
  }

  return(c(lambda.0, c.0))
}