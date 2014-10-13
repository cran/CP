LikelihoodNonMixExp <- function(data1, data2, data,
                                lambda1.0, c1.0,
                                lambda2.0, c2.0,
                                lambda.0, c.0) {

  ## Calculates the maximum likelihood estimators
  ## of the parameters lambda and c in the non-mixture model
  ## with exponential survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t)], lambda > 0, 0 < c < 1, t >= 0.
  ##
  ## Args:
  ##   data1: Data frame which consists of at least three columns
  ##          with the group in the first (all of group 1),
  ##          status (1 = event, 0 = censored) in the second
  ##          and event time in the third column.
  ##   data2: Data frame which consists of at least three columns
  ##          with the group in the first (all of group 2),
  ##          status (1 = event, 0 = censored) in the second
  ##          and event time in the third column.
  ##   data: Data frame which consists of at least three columns
  ##         with the group in the first,
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   lambda1.0: Initial value for the estimate of the parameter lambda in group 1.
  ##   c1.0: Initial value for the estimate of the parameter c in group 1.
  ##   lambda2.0: Initial value for the estimate of the parameter lambda in group 2.
  ##   c2.0: Initial value for the estimate of the parameter c in group 2.
  ##   lambda.0: Initial value for the estimate of the parameter lambda for all data.
  ##   c.0: Initial value for the estimate of the parameter c for all data.
  ##
  ## Returns:
  ##   Returns the maximum likelihood estimates for the parameters
  ##   lambda and c, i. e. estimates within the ordinary model
  ##   and estimates under the proportional hazard assumption.

  # GROUP 1 (data1) -> lambda1.hat, c1.hat
  # log-likelihood function of parameters lambda and c for data of group 1
  l1 <- function(x) {
    (log(x[1]) * sum(data1[, 2])
      - x[1] * sum(data1[, 2] * data1[, 3])
      + log(- log(x[2])) * sum(data1[, 2])
      + log(x[2]) * sum(1 - exp(- x[1] * data1[, 3])))
  }
  # gradient of log-likelihood function
  l1Grad <- function(x) {
    c(1 / x[1] * sum(data1[, 2])
       - sum(data1[, 2] * data1[, 3])
       + log(x[2]) * sum(data1[, 3] * exp(- x[1] * data1[, 3])),
      1 / (x[2] * log(x[2])) * sum(data1[, 2])
       + 1 / x[2] * sum(1 - exp(- x[1] * data1[, 3])))
  }
  # constraints for parameters lambda and c for data of group 1
  eps <- 1e-7
  A1  <- matrix(data  = c(1, 0,
                          0, 1,
                          0, - 1),
                nrow  = 3,
                ncol  = 2,
                byrow = TRUE)
  b1  <- c(eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters lambda and c for data of group 1
  res1 <- constrOptim(theta   = c(lambda1.0, c1.0),
                      f       = l1,
                      grad    = l1Grad,
                      ui      = A1,
                      ci      = b1,
                      control = list(fnscale = -1))

  # GROUP 2 (data2) -> lambda2.hat, c2.hat
  # log-likelihood function of parameters lambda and c for data of group 2
  l2 <- function(x) {
    (log(x[1]) * sum(data2[, 2])
      - x[1] * sum(data2[, 2] * data2[, 3])
      + log(- log(x[2])) * sum(data2[, 2])
      + log(x[2]) * sum(1 - exp(- x[1] * data2[, 3])))
  }
  # gradient of log-likelihood function
  l2Grad <- function(x) {
    c(1 / x[1] * sum(data2[, 2])
       - sum(data2[, 2] * data2[, 3])
       + log(x[2]) * sum(data2[, 3] * exp(- x[1] * data2[, 3])),
      1 / (x[2] * log(x[2])) * sum(data2[, 2])
       + 1 / x[2] * sum(1 - exp(- x[1] * data2[, 3])))
  }
  # constraints for parameters lambda and c for data of group 2
  eps <- 1e-7
  A2  <- matrix(data  = c(1, 0,
                          0, 1,
                          0, - 1),
                nrow  = 3,
                ncol  = 2,
                byrow = TRUE)
  b2  <- c(eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters lambda and c for data of group 2
  res2 <- constrOptim(theta   = c(lambda2.0, c2.0),
                      f       = l2,
                      grad    = l2Grad,
                      ui      = A2,
                      ci      = b2,
                      control = list(fnscale = -1))

  # ALL DATA (data) -> lambda.hat
  # log-likelihood function of parameters lambda and c for all data
  lAll <- function(x) {
    (log(x[1]) * sum(data[, 2])
      - x[1] * sum(data[, 2] * data[, 3])
      + log(- log(x[2])) * sum(data[, 2])
      + log(x[2]) * sum(1 - exp(- x[1] * data[, 3])))
  }
  # gradient of log-likelihood function
  lAllGrad <- function(x) {
    c(1 / x[1] * sum(data[, 2])
       - sum(data[, 2] * data[, 3])
       + log(x[2]) * sum(data[, 3] * exp(- x[1] * data[, 3])),
      1 / (x[2] * log(x[2])) * sum(data[, 2])
       + 1 / x[2] * sum(1 - exp(- x[1] * data[, 3])))
  }
  # constraints for parameters lambda and c for all data
  eps  <- 1e-7
  AAll <- matrix(data  = c(1, 0,
                           0, 1,
                           0, - 1),
                 nrow  = 3,
                 ncol  = 2,
                 byrow = TRUE)
  bAll <- c(eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters lambda and c for all data
  resAll <- constrOptim(theta   = c(lambda.0, c.0),
                        f       = lAll,
                        grad    = lAllGrad,
                        ui      = AAll,
                        ci      = bAll,
                        control = list(fnscale = -1))
  lambda.hat <- resAll$par[1]

  # GROUP 1 (data1) with fixed lambda -> c1.cond.hat
  # log-likelihood function of parameter c for data of group 1
  l1.cond <- function(x) {
    (log(lambda.hat) * sum(data1[, 2])
      - lambda.hat * sum(data1[, 2] * data1[, 3])
      + log(- log(x)) * sum(data1[, 2])
      + log(x) * sum(1 - exp(- lambda.hat * data1[, 3])))
  }
  # gradient of log-likelihood function
  l1Grad.cond <- function(x) {
    (1 / (x * log(x)) * sum(data1[, 2])
      + 1 / x * sum(1 - exp(- lambda.hat * data1[, 3])))
  }
  # constraints for parameter c for data of group 1
  eps     <- 1e-7
  A1.cond <- matrix(data  = c(1,
                              - 1),
                    nrow  = 2,
                    ncol  = 1,
                    byrow = TRUE)
  b1.cond <- c(eps, - (1 - eps))
  # maximum likelihood estimator for parameter c for data of group 1
  res1.cond <- constrOptim(theta   = c1.0,
                           f       = l1.cond,
                           grad    = l1Grad.cond,
                           ui      = A1.cond,
                           ci      = b1.cond,
                           control = list(fnscale = -1))

  # GROUP 2 (data2) with fixed lambda -> c2.cond.hat
  # log-likelihood function of parameter c for data of group 2
  l2.cond <- function(x) {
    (log(lambda.hat) * sum(data2[, 2])
      - lambda.hat * sum(data2[, 2] * data2[, 3])
      + log(- log(x)) * sum(data2[, 2])
      + log(x) * sum(1 - exp(- lambda.hat * data2[, 3])))
  }
  # gradient of log-likelihood function
  l2Grad.cond <- function(x) {
    (1 / (x * log(x)) * sum(data2[, 2])
      + 1 / x * sum(1 - exp(- lambda.hat * data2[, 3])))
  }
  # constraints for parameter c for data of group 2
  eps     <- 1e-7
  A2.cond <- matrix(data  = c(1,
                              - 1),
                    nrow  = 2,
                    ncol  = 1,
                    byrow = TRUE)
  b2.cond <- c(eps, - (1 - eps))
  # maximum likelihood estimator for parameter c for data of group 2
  res2.cond <- constrOptim(theta   = c2.0,
                           f       = l2.cond,
                           grad    = l2Grad.cond,
                           ui      = A2.cond,
                           ci      = b2.cond,
                           control = list(fnscale = -1))

  return(c(res1$par[1], res1$par[2],
           res2$par[1], res2$par[2],
           resAll$par[1],
           res1.cond$par,
           res2.cond$par,
           res1.cond$value,
           res2.cond$value))
}