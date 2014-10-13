LikelihoodNonMixGamma <- function(data1, data2, data,
                                  a1.0, b1.0, c1.0,
                                  a2.0, b2.0, c2.0,
                                  a.0, b.0, c.0) {

  ## Calculates the maximum likelihood estimators
  ## of the parameters a, b and c in the non-mixture model
  ## with Gamma type survival, i. e.
  ##   S(t) = c^Gamma^(0)(a, b * t), a > 0, b > 0, 0 < c < 1, t >= 0,
  ## with Gamma^(0) being the regularized incomplete Gamma function
  ## of the upper bound.
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
  ##   a1.0: Initial value for the estimate of the parameter a in group 1.
  ##   b1.0: Initial value for the estimate of the parameter b in group 1.
  ##   c1.0: Initial value for the estimate of the parameter c in group 1.
  ##   a2.0: Initial value for the estimate of the parameter a in group 2.
  ##   b2.0: Initial value for the estimate of the parameter b in group 2.
  ##   c2.0: Initial value for the estimate of the parameter c in group 2.
  ##   a.0: Initial value for the estimate of the parameter a for all data.
  ##   b.0: Initial value for the estimate of the parameter b for all data.
  ##   c.0: Initial value for the estimate of the parameter c for all data.
  ##
  ## Returns:
  ##   Returns the maximum likelihood estimates for the parameters
  ##   a, b and c, i. e. estimates within the ordinary model
  ##   and estimates under the proportional hazard assumption.

  # GROUP 1 (data1) -> a1.hat, b1.hat, c1.hat
  # log-likelihood function of parameters a, b and c for data of group 1
  l1 <- function(x) {
    (x[1] * log(x[2]) * sum(data1[, 2])
      - lgamma(x = x[1]) * sum(data1[, 2])
      + (x[1] - 1) * sum(data1[, 2] * log(data1[, 3]))
      - x[2] * sum(data1[, 2] * data1[, 3])
      + log(- log(x[3])) * sum(data1[, 2])
      + log(x[3]) * sum(pgamma(q     = data1[, 3],
                               shape = x[1],
                               rate  = x[2])))
  }
  # auxiliary function
  Integrand1 <- function(x, par) {
    log(x) * x^(par - 1) * exp(-x)
  }
  # gradient of log-likelihood function
  l1Grad <- function(x) {
    c(log(x[2]) * sum(data1[, 2])
       - digamma(x = x[1]) * sum(data1[, 2])
       + sum(data1[, 2] * log(data1[, 3]))
       + log(x[3]) * sum((gamma(x = x[1])
            * integrate(f     = Integrand1,
                        lower = 0,
                        upper = x[2] * data1[, 3],
                        par   = x[1])$value
           - gamma(x = x[1]) * pgamma(q     = data1[, 3],
                                      shape = x[1],
                                      rate  = x[2])
            * integrate(f     = Integrand1,
                        lower = 0,
                        upper = Inf,
                        par   = x[1])$value)
          / gamma(x = x[1])^2),
      x[1] / x[2] * sum(data1[, 2])
       - sum(data1[, 2] * data1[, 3])
       + log(x[3]) * sum((x[2] * data1[, 3])^(x[1] - 1) * exp(- x[2] * data1[, 3]) * data1[, 3]
          / gamma(x = x[1])),
      1 / (x[3] * log(x[3])) * sum(data1[, 2])
       + 1 / x[3] * sum(pgamma(q     = data1[, 3],
                               shape = x[1],
                               rate  = x[2])))
  }
  # constraints for parameters a, b and c for data of group 1
  eps <- 1e-7
  A1  <- matrix(data  = c(1, 0, 0,
                          0, 1, 0,
                          0, 0, 1,
                          0, 0, - 1),
                nrow  = 4,
                ncol  = 3,
                byrow = TRUE)
  b1  <- c(eps, eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters a, b and c for data of group 1
  res1 <- constrOptim(theta   = c(a1.0, b1.0, c1.0),
                      f       = l1,
                      grad    = l1Grad,
                      ui      = A1,
                      ci      = b1,
                      control = list(fnscale = -1))

  # GROUP 2 (data2) -> a2.hat, b2.hat, c2.hat
  # log-likelihood function of parameters a, b and c for data of group 2
  l2 <- function(x) {
    (x[1] * log(x[2]) * sum(data2[, 2])
      - lgamma(x = x[1]) * sum(data2[, 2])
      + (x[1] - 1) * sum(data2[, 2] * log(data2[, 3]))
      - x[2] * sum(data2[, 2] * data2[, 3])
      + log(- log(x[3])) * sum(data2[, 2])
      + log(x[3]) * sum(pgamma(q     = data2[, 3],
                               shape = x[1],
                               rate  = x[2])))
  }
  # auxiliary function
  Integrand2 <- function(x, par) {
    log(x) * x^(par - 1) * exp(-x)
  }
  # gradient of log-likelihood function
  l2Grad <- function(x) {
    c(log(x[2]) * sum(data2[, 2])
       - digamma(x = x[1]) * sum(data2[, 2])
       + sum(data2[, 2] * log(data2[, 3]))
       + log(x[3]) * sum((gamma(x = x[1])
            * integrate(f     = Integrand2,
                        lower = 0,
                        upper = x[2] * data2[, 3],
                        par   = x[1])$value
           - gamma(x = x[1]) * pgamma(q     = data2[, 3],
                                      shape = x[1],
                                      rate  = x[2])
            * integrate(f     = Integrand2,
                        lower = 0,
                        upper = Inf,
                        par   = x[1])$value)
          / gamma(x = x[1])^2),
      x[1] / x[2] * sum(data2[, 2])
       - sum(data2[, 2] * data2[, 3])
       + log(x[3]) * sum((x[2] * data2[, 3])^(x[1] - 1) * exp(- x[2] * data2[, 3]) * data2[, 3]
          / gamma(x = x[1])),
      1 / (x[3] * log(x[3])) * sum(data2[, 2])
       + 1 / x[3] * sum(pgamma(q     = data2[, 3],
                               shape = x[1],
                               rate  = x[2])))
  }
  # constraints for parameters a, b and c for data of group 2
  eps <- 1e-7
  A2  <- matrix(data  = c(1, 0, 0,
                          0, 1, 0,
                          0, 0, 1,
                          0, 0, - 1),
                nrow  = 4,
                ncol  = 3,
                byrow = TRUE)
  b2  <- c(eps, eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters a, b and c for data of group 2
  res2 <- constrOptim(theta   = c(a2.0, b2.0, c2.0),
                      f       = l2,
                      grad    = l2Grad,
                      ui      = A2,
                      ci      = b2,
                      control = list(fnscale = -1))

  # ALL DATA (data) -> a.hat, b.hat
  # log-likelihood function of parameters a, b and c for all data
  lAll <- function(x) {
    (x[1] * log(x[2]) * sum(data[, 2])
      - lgamma(x = x[1]) * sum(data[, 2])
      + (x[1] - 1) * sum(data[, 2] * log(data[, 3]))
      - x[2] * sum(data[, 2] * data[, 3])
      + log(- log(x[3])) * sum(data[, 2])
      + log(x[3]) * sum(pgamma(q     = data[, 3],
                               shape = x[1],
                               rate  = x[2])))
  }
  # auxiliary function
  IntegrandAll <- function(x, par) {
    log(x) * x^(par - 1) * exp(-x)
  }
  # gradient of log-likelihood function
  lAllGrad <- function(x) {
    c(log(x[2]) * sum(data[, 2])
       - digamma(x = x[1]) * sum(data[, 2])
       + sum(data[, 2] * log(data[, 3]))
       + log(x[3]) * sum((gamma(x = x[1])
            * integrate(f     = IntegrandAll,
                        lower = 0,
                        upper = x[2] * data[, 3],
                        par   = x[1])$value
           - gamma(x = x[1]) * pgamma(q     = data[, 3],
                                      shape = x[1],
                                      rate  = x[2])
            * integrate(f     = IntegrandAll,
                        lower = 0,
                        upper = Inf,
                        par   = x[1])$value)
          / gamma(x = x[1])^2),
      x[1] / x[2] * sum(data[, 2])
       - sum(data[, 2] * data[, 3])
       + log(x[3]) * sum((x[2] * data[, 3])^(x[1] - 1) * exp(- x[2] * data[, 3]) * data[, 3]
          / gamma(x = x[1])),
      1 / (x[3] * log(x[3])) * sum(data[, 2])
       + 1 / x[3] * sum(pgamma(q     = data[, 3],
                               shape = x[1],
                               rate  = x[2])))
  }
  # constraints for parameters a, b and c for all data
  eps  <- 1e-7
  AAll <- matrix(data  = c(1, 0, 0,
                           0, 1, 0,
                           0, 0, 1,
                           0, 0, - 1),
                 nrow  = 4,
                 ncol  = 3,
                 byrow = TRUE)
  bAll <- c(eps, eps, eps, - (1 - eps))
  # maximum likelihood estimators for parameters a, b and c for all data
  resAll <- constrOptim(theta   = c(a.0, b.0, c.0),
                        f       = lAll,
                        grad    = lAllGrad,
                        ui      = AAll,
                        ci      = bAll,
                        control = list(fnscale = -1))
  a.hat <- resAll$par[1]
  b.hat <- resAll$par[2]

  # GROUP 1 (data1) with fixed a and b -> c1.cond.hat
  # log-likelihood function of parameter c for data of group 1
  l1.cond <- function(x) {
    (a.hat * log(b.hat) * sum(data1[, 2])
      - lgamma(x = a.hat) * sum(data1[, 2])
      + (a.hat - 1) * sum(data1[, 2] * log(data1[, 3]))
      - b.hat * sum(data1[, 2] * data1[, 3])
      + log(- log(x)) * sum(data1[, 2])
      + log(x) * sum(pgamma(q     = data1[, 3],
                            shape = a.hat,
                            rate  = b.hat)))
  }
  # auxiliary function
  Integrand1.cond <- function(x, par) {
    log(x) * x^(par - 1) * exp(-x)
  }
  # gradient of log-likelihood function
  l1Grad.cond <- function(x) {
    (1 / (x * log(x)) * sum(data1[, 2])
      + 1 / x * sum(pgamma(q     = data1[, 3],
                    shape = a.hat,
                    rate  = b.hat)))
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

  # GROUP 2 (data2) with fixed a and b -> c2.cond.hat
  # log-likelihood function of parameter c for data of group 2
  l2.cond <- function(x) {
    (a.hat * log(b.hat) * sum(data2[, 2])
      - lgamma(x = a.hat) * sum(data2[, 2])
      + (a.hat - 1) * sum(data2[, 2] * log(data2[, 3]))
      - b.hat * sum(data2[, 2] * data2[, 3])
      + log(- log(x)) * sum(data2[, 2])
      + log(x) * sum(pgamma(q     = data2[, 3],
                            shape = a.hat,
                            rate  = b.hat)))
  }
  # auxiliary function
  Integrand2.cond <- function(x, par) {
    log(x) * x^(par - 1) * exp(-x)
  }
  # gradient of log-likelihood function
  l2Grad.cond <- function(x) {
    (1 / (x * log(x)) * sum(data2[, 2])
      + 1 / x * sum(pgamma(q     = data2[, 3],
                           shape = a.hat,
                           rate  = b.hat)))
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