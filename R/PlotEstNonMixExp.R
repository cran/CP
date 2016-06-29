PlotEstNonMixExp <- function(data1, data2,
                             lambda1.hat, c1.hat,
                             lambda2.hat, c2.hat,
                             group1.name, group2.name) {

  ## Plots the survival curves
  ## under the non-mixture model with exponential survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t)], lambda > 0, 0 < c < 1, t >= 0,
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
  ##   lambda1.hat: Maximum likelihood estimator for lambda1.
  ##   c1.hat: Maximum likelihood estimator for c1.
  ##   lambda2.hat: Maximum likelihood estimator for lambda2.
  ##   c2.hat: Maximum likelihood estimator for c2.
  ##   group1.name: Expression for group 1.
  ##   group2.name: Expression for group 2.
  ##
  ## Results:
  ##   Plots the survival curves
  ##   under the non-mixture model with exponential survival
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
  S1 <- c1.hat^(1 - exp(- lambda1.hat * t1))
  S2 <- c2.hat^(1 - exp(- lambda2.hat * t2))

  # plot of survival curve of group 1
  graphics::lines(x    = t1,
                  y    = S1,
                  col  = "blue")

  # plot of survival curve of group 2
  graphics::lines(x   = t2,
                  y   = S2,
                  col = "green")

  graphics::legend(x      = "topright",
                   legend = c(group1.name, group2.name),
                   col    = c("blue", "green"),
                   lty    = c(1, 1),
                   pch    = c(3, 3),
                   bg     = "white")
}