CalcConPwrExp <- function(theta.0,
                          d1, o1, O1.star, lambda1.hat,
                          d2, o2, O2.star,
                          n.star,
                          alpha) {

  ## Calculates the conditional power
  ## for the exponential model, i. e.
  ##  S(t) = exp(- lambda * t), lambda > 0, t >= 0.
  ##
  ## Args:
  ##   Parameters from exponential power calculations.
  ##
  ## Results:
  ##   Returns conditional power values.

  # range of theta for conditional power function
  min <- exp(min(log(theta.0), - log(theta.0)) - 1)
  max <- exp(max(log(theta.0), - log(theta.0)) + 1)

  # choice of suitable stepwidth for theta's
  if ((max - min) / 0.01 <= 1000) {
    theta <- seq(from = min,
                 to   = max,
                 by   = 0.01)
  }
  else {
    theta <- seq(from       = min,
                 to         = max,
                 length.out = 1000)
  }

  # expected value under the null hypothesis
  mu.theta.null     <- (log((d2 + lambda1.hat * O2.star) / (o2 + O2.star))
                         - log((d1 + lambda1.hat * O1.star) / (o1 + O1.star)))
  # standard deviation under the null hypothesis
  sigma.theta.null  <- sqrt(2 / n.star)
  # expected value under the alternative hypothesis
  mu.theta.alter    <- (log((d2 + theta * lambda1.hat * O2.star) / (o2 + O2.star))
                         - log((d1 + lambda1.hat * O1.star) / (o1 + O1.star)))
  # standard deviation under the alternative hypothesis
  sigma.theta.alter <- sqrt(2 / n.star)
  # (asymptotically) conditional power function
  gamma.theta       <- (stats::pnorm((stats::qnorm(alpha / 2) * sigma.theta.null + mu.theta.null
                                - mu.theta.alter) / sigma.theta.alter)
                         + 1 - stats::pnorm((stats::qnorm(1 - alpha / 2) * sigma.theta.null + mu.theta.null
                                       - mu.theta.alter) / sigma.theta.alter))
  # values at theta.0
  mu.theta.0        <- (log((d2 + theta.0 * lambda1.hat * O2.star) / (o2 + O2.star))
                         - log((d1 + lambda1.hat * O1.star) / (o1 + O1.star)))
  sigma.theta.0     <- sqrt(2 / n.star)
  gamma.theta.0     <- (stats::pnorm((stats::qnorm(alpha / 2) * sigma.theta.null + mu.theta.null
                                - mu.theta.0) / sigma.theta.0)
                         + 1 - stats::pnorm((stats::qnorm(1 - alpha / 2) * sigma.theta.null + mu.theta.null
                                       - mu.theta.0) / sigma.theta.0))

  return(list(theta, gamma.theta, gamma.theta.0))
}