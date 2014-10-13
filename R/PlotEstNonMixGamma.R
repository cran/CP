PlotEstNonMixGamma <- function(data1, data2,
                               a1.hat, b1.hat, c1.hat,
                               a2.hat, b2.hat, c2.hat,
                               group1.name, group2.name) {

  ## Plots the survival curves
  ## under the non-mixture model with Gamma type survival, i. e.
  ##   S(t) = c^Gamma^(0)(a, b * t), a > 0, b > 0, 0 < c < 1, t >= 0,
  ## with Gamma^(0) being the regularized incomplete Gamma function of the upper bound,
  ## with the estimated parameters.
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
  ##   a1.hat: Maximum likelihood estimator for a1.
  ##   b1.hat: Maximum likelihood estimator for b1.
  ##   c1.hat: Maximum likelihood estimator for c1.
  ##   a2.hat: Maximum likelihood estimator for a2.
  ##   b1.hat: Maximum likelihood estimator for b2.
  ##   c2.hat: Maximum likelihood estimator for c2.
  ##   group1.name: Expression for group 1.
  ##   group2.name: Expression for group 2.
  ##
  ## Results:
  ##   Plots the survival curves
  ##   under the non-mixture model with Gamma type survival
  ##   in an 1 x 2 plot array
  ##   in which the first component consists of the Kaplan-Meier curves
  ##   and the second component is still empty.

  # range of time for survival curves
  # of group 1 and group 2
  max1 <- max(data1[, 3])
  max2 <- max(data2[, 3])

  # choice of suitable stepwidth of time variables
  # of group 1 and group 2
  if (max1 / 0.01 <= 1000) {
    t1 <- seq(from = 0,
              to   = max1,
              by   = 0.01)
  }
  else {
    t1 <- seq(from       = 0,
              to         = max1,
              length.out = 1000)
  }
  if (max2 / 0.01 <= 1000) {
    t2 <- seq(from = 0,
              to   = max2,
              by   = 0.01)
  }
  else {
    t2 <- seq(from       = 0,
              to         = max2,
              length.out = 1000)
  }

  # survival curves of group 1 and group 2
  S1 <- c1.hat^pgamma(q     = t1,
                      shape = a1.hat,
                      rate  = b1.hat)
  S2 <- c2.hat^pgamma(q     = t2,
                      shape = a2.hat,
                      rate  = b2.hat)

  # plot of survival curve of group 1
  lines(x    = t1,
        y    = S1,
        col  = "blue")

  # plot of survival curve of group 2
  lines(x   = t2,
        y   = S2,
        col = "green")

  legend(x      = "topright",
         legend = c(group1.name, group2.name),
         col    = c("blue", "green"),
         lty    = c(1, 1),
         pch    = c(3, 3),
         bg     = "white")
}