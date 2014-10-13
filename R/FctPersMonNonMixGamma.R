FctPersMonNonMixGamma <- function(data, a.hat, b.hat, group.name) {

  ## Calculates the values of some function of the person months
  ## in the non-mixture model with Gamma type survival, i. e.
  ##   S(t) = c^Gamma^(0)(a, b * t), a > 0, b > 0, 0 < c < 1, t >= 0,
  ## with Gamma^(0) being the regularized incomplete Gamma function
  ## of the upper bound.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns for the group
  ##         in the first (all of the same group),
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   a.hat: Estimator of the parameter lambda.
  ##   b.hat: Estimator of the parameter k.
  ##   group.name: Name of the group.
  ##
  ## Results:
  ##   Returns the value of some function of the person months.

  # function of person months, and their verification
  o.stroke <- sum(pgamma(q     = data[, 3],
                         shape = a.hat,
                         rate  = b.hat))
  if (o.stroke <= 0) {
    stop(paste("Number of person months in", group.name, "must be bigger than 0.",
         call. = FALSE))
  }

  return(o.stroke)
}