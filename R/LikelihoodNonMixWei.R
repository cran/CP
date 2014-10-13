LikelihoodNonMixWei <- function(data1, data2, data,
                                lambda1.0, k1.0, c1.0,
                                lambda2.0, k2.0, c2.0,
                                lambda.0, k.0, c.0) {

  ## Calculates the maximum likelihood estimators
  ## of the parameters lambda, k and c in the non-mixture model
  ## with Weibull type survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t^k)], lambda > 0, k > 0, 0 < c < 1, t >= 0.
  ##
  ## Args:
  ##   data1: Data frame which consists of at least three columns with the group
  ##          in the first (all of group 1),
  ##          status (1 = event, 0 = censored) in the second
  ##          and event time in the third column.
  ##   data2: Data frame which consists of at least three columns with the group
  ##          in the first (all of group 2),
  ##          status (1 = event, 0 = censored) in the second
  ##          and event time in the third column.
  ##   data: Data frame which consists of at least three columns with the group
  ##         in the first (all of the same group),
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   lambda1.0: Initial value for the estimate of the parameter lambda in group 1.
  ##   k1.0: Initial value for the estimate of the parameter k in group 1.
  ##   c1.0: Initial value for the estimate of the parameter c in group 1.
  ##   lambda2.0: Initial value for the estimate of the parameter lambda in group 2.
  ##   k2.0: Initial value for the estimate of the parameter k in group 2.
  ##   c2.0: Initial value for the estimate of the parameter c in group 2.
  ##   lambda.0: Initial value for the estimate of the parameter lambda for all data.
  ##   k.0: Initial value for the estimate of the parameter k for all data.
  ##   c.0: Initial value for the estimate of the parameter c for all data.
  ##
  ## Returns:
  ##   Returns the maximum likelihood estimates for the parameters
  ##   lambda, k and c, i. e. estimates within the ordinary model
  ##   and estimates under the proportional hazard assumption.

  # GROUP 1 (data1) -> lambda1.hat, k1.hat, c1.hat
  # log-likelihood function of parameters lambda, k and c for data of group 1
  l1 <- function(x) {
    (log(x[1]) * sum(data1[, 2])
      + log(x[2]) * sum(data1[, 2])
      + (x[2] - 1) * sum(data1[, 2] * log(data1[, 3]))
      - x[1] * sum(data1[, 2] * data1[, 3]^x[2])
      + log(- log(x[3])) * sum(data1[, 2])
      + log(x[3]) * sum(1 - exp(- x[1] * data1[, 3]^x[2])))
  }
  # gradient of log-likelihood function
  l1Grad <- function(x) {
    c(1 / x[1] * sum(data1[, 2])
       - sum(data1[, 2] * data1[, 3]^x[2])
       + log(x[3]) * sum(data1[, 3]^x[2] * exp(- x[1] * data1[, 3]^x[2])),
      1 / x[2] * sum(data1[, 2])
       + sum(data1[, 2] * log(data1[, 3]))
       - x[1] * sum(data1[, 2] * data1[, 3]^x[2] * log(data1[, 3]))
       + x[1] * log(x[3]) * sum(exp(- x[1] * data1[, 3]^x[2]) * data1[, 3]^x[2] * log(data1[, 3])),
      1 / (x[3] * log(x[3])) * sum(data1[, 2])
       + 1 / x[3] * sum(1 - exp(- x[1] * data1[, 3]^x[2])))
  }
  # constraints for parameters lambda, k and c for data of group 1
  eps <- 1e-7
  A1  <- matrix(data  = c(1, 0, 0,
                          0, 1, 0,
                          0, 0, 1,
                          0, 0, - 1),
                nrow  = 4,
                ncol  = 3,
                byrow = TRUE)
  b1  <- c(eps, eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters lambda, k and c for data of group 1
  res1 <- constrOptim(theta   = c(lambda1.0, k1.0, c1.0),
                      f       = l1,
                      grad    = l1Grad,
                      ui      = A1,
                      ci      = b1,
                      control = list(fnscale = -1))

  # GROUP 2 (data2) -> lambda2.hat, k2.hat, c2.hat
  # log-likelihood function of parameters lambda, k and c for data of group 2
  l2 <- function(x) {
    (log(x[1]) * sum(data2[, 2])
      + log(x[2]) * sum(data2[, 2])
      + (x[2] - 1) * sum(data2[, 2] * log(data2[, 3]))
      - x[1] * sum(data2[, 2] * data2[, 3]^x[2])
      + log(- log(x[3])) * sum(data2[, 2])
      + log(x[3]) * sum(1 - exp(- x[1] * data2[, 3]^x[2])))
  }
  # gradient of log-likelihood function
  l2Grad <- function(x) {
    c(1 / x[1] * sum(data2[, 2])
       - sum(data2[, 2] * data2[, 3]^x[2])
       + log(x[3]) * sum(data2[, 3]^x[2] * exp(- x[1] * data2[, 3]^x[2])),
      1 / x[2] * sum(data2[, 2])
       + sum(data2[, 2] * log(data2[, 3]))
       - x[1] * sum(data2[, 2] * data2[, 3]^x[2] * log(data2[, 3]))
       + x[1] * log(x[3]) * sum(exp(- x[1] * data2[, 3]^x[2]) * data2[, 3]^x[2] * log(data2[, 3])),
      1 / (x[3] * log(x[3])) * sum(data2[, 2])
       + 1 / x[3] * sum(1 - exp(- x[1] * data2[, 3]^x[2])))
  }
  # constraints for parameters lambda, k and c for data of group 2
  eps <- 1e-7
  A2  <- matrix(data  = c(1, 0, 0,
                          0, 1, 0,
                          0, 0, 1,
                          0, 0, - 1),
                nrow  = 4,
                ncol  = 3,
                byrow = TRUE)
  b2  <- c(eps, eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters lambda, k and c for data of group 2
  res2 <- constrOptim(theta   = c(lambda2.0, k2.0, c2.0),
                      f       = l2,
                      grad    = l2Grad,
                      ui      = A2,
                      ci      = b2,
                      control = list(fnscale = -1))

  # ALL DATA (data) -> lambda.hat, k.hat
  # log-likelihood function of parameters lambda, k and c for all data
  lAll <- function(x) {
    (log(x[1]) * sum(data[, 2])
      + log(x[2]) * sum(data[, 2])
      + (x[2] - 1) * sum(data[, 2] * log(data[, 3]))
      - x[1] * sum(data[, 2] * data[, 3]^x[2])
      + log(- log(x[3])) * sum(data[, 2])
      + log(x[3]) * sum(1 - exp(- x[1] * data[, 3]^x[2])))
  }
  # gradient of log-likelihood function
  lAllGrad <- function(x) {
    c(1 / x[1] * sum(data[, 2])
       - sum(data[, 2] * data[, 3]^x[2])
       + log(x[3]) * sum(data[, 3]^x[2] * exp(- x[1] * data[, 3]^x[2])),
      1 / x[2] * sum(data[, 2])
       + sum(data[, 2] * log(data[, 3]))
       - x[1] * sum(data[, 2] * data[, 3]^x[2] * log(data[, 3]))
       + x[1] * log(x[3]) * sum(exp(- x[1] * data[, 3]^x[2]) * data[, 3]^x[2] * log(data[, 3])),
      1 / (x[3] * log(x[3])) * sum(data[, 2])
       + 1 / x[3] * sum(1 - exp(- x[1] * data[, 3]^x[2])))
  }
  # constraints for parameters lambda, k and c for all data
  eps  <- 1e-7
  AAll <- matrix(data  = c(1, 0, 0,
                           0, 1, 0,
                           0, 0, 1,
                           0, 0, - 1),
                 nrow  = 4,
                 ncol  = 3,
                 byrow = TRUE)
  bAll <- c(eps, eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters lambda, k and c for all data
  resAll <- constrOptim(theta   = c(lambda.0, k.0, c.0),
                        f       = lAll,
                        grad    = lAllGrad,
                        ui      = AAll,
                        ci      = bAll,
                        control = list(fnscale = -1))
  lambda.hat <- resAll$par[1]
  k.hat      <- resAll$par[2]

  # GROUP 1 (data1) with fixed lambda and k -> c1.cond.hat
  # log-likelihood function of parameter c for data of group 1
  l1.cond <- function(x) {
    (log(lambda.hat) * sum(data1[, 2])
      + log(k.hat) * sum(data1[, 2])
      + (k.hat - 1) * sum(data1[, 2] * log(data1[, 3]))
      - lambda.hat * sum(data1[, 2] * data1[, 3]^k.hat)
      + log(- log(x)) * sum(data1[, 2])
      + log(x) * sum(1 - exp(- lambda.hat * data1[, 3]^k.hat)))
  }
  # gradient of log-likelihood function
  l1Grad.cond <- function(x) {
    (1 / (x * log(x)) * sum(data1[, 2])
      + 1 / x * sum(1 - exp(- lambda.hat * data1[, 3]^k.hat)))
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

  # GROUP 2 (data2) with fixed lambda and k -> c2.cond.hat
  # log-likelihood function of parameter c for data of group 2
  l2.cond <- function(x) {
    (log(lambda.hat) * sum(data2[, 2])
      + log(k.hat) * sum(data2[, 2])
      + (k.hat - 1) * sum(data2[, 2] * log(data2[, 3]))
      - lambda.hat * sum(data2[, 2] * data2[, 3]^k.hat)
      + log(- log(x)) * sum(data2[, 2])
      + log(x) * sum(1 - exp(- lambda.hat * data2[, 3]^k.hat)))
  }
  # gradient of log-likelihood function
  l2Grad.cond <- function(x) {
    (1 / (x * log(x)) * sum(data2[, 2])
      + 1 / x * sum(1 - exp(- lambda.hat * data2[, 3]^k.hat)))
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

  return(c(res1$par[1], res1$par[2], res1$par[3],
           res2$par[1], res2$par[2], res2$par[3],
           resAll$par[1], resAll$par[2],
           res1.cond$par,
           res2.cond$par,
           res1.cond$value,
           res2.cond$value))
}