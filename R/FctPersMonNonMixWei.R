FctPersMonNonMixWei <- function(data, lambda.hat, k.hat, group.name) {

  ## Calculates the values of some function of the patient time
  ## in the non-mixture model with Weibull type survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t^k)], lambda > 0, k > 0, 0 < c < 1, t >= 0.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns for the group
  ##         in the first (all of the same group),
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   lambda.hat: Estimator of the parameter lambda.
  ##   k.hat: Estimator of the parameter k.
  ##   group.name: Name of the group.
  ##
  ## Results:
  ##   Returns the value of some function of the patient time.

  # function of patient time, and their verification
  o.stroke <- sum(1 - exp(- lambda.hat * data[, 3]^k.hat))
  if (o.stroke <= 0) {
    stop(paste("Number of person months in", group.name, "must be bigger than 0.",
         call. = FALSE))
  }

  return(o.stroke)
}