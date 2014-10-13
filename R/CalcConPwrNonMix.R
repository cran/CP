CalcConPwrNonMix <- function(theta.0,
                             d1, o1.stroke, O1.stroke.star, c1.hat,
                             d2, o2.stroke, O2.stroke.star, O2.stroke.star.null,
                             n.star,
                             alpha) {

  ## Calculates the conditional power
  ## for the non-mixture model.
  ##
  ## Args:
  ##   Parameters from power calculations.
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
  mu.theta.null     <- (log((d2 - log(c1.hat) * O2.stroke.star.null) / (o2.stroke + O2.stroke.star.null))
                         - log((d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star)))
  # standard deviation under the null hypothesis
  sigma.theta.null  <- (sqrt(1 / (n.star
                                   * (1 - exp(- (d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star))
                                           * (1 + ((d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star))^2)))
                              + 1 / (n.star
                                      * (1 - exp(- (d2 - log(c1.hat) * O2.stroke.star.null) / (o2.stroke + O2.stroke.star.null))
                                              * (1 + ((d2 - log(c1.hat) * O2.stroke.star.null) / (o2.stroke + O2.stroke.star.null))^2)))))
  # expected value under the alternative hypothesis
  mu.theta.alter    <- (log((d2 - theta * log(c1.hat) * O2.stroke.star) / (o2.stroke + O2.stroke.star))
                         - log((d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star)))
  # standard deviation under the alternative hypothesis
  sigma.theta.alter <- (sqrt(1 / (n.star
                                   * (1 - exp(- (d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star))
                                           * (1 + ((d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star))^2)))
                              + 1 / (n.star
                                      * (1 - exp(- (d2 - theta * log(c1.hat) * O2.stroke.star) / (o2.stroke + O2.stroke.star))
                                              * (1 + ((d2 - theta * log(c1.hat) * O2.stroke.star) / (o2.stroke + O2.stroke.star))^2)))))
  # (asymptotically) conditional power function
  gamma.theta       <- (pnorm((qnorm(alpha / 2) * sigma.theta.null + mu.theta.null
                                - mu.theta.alter) / sigma.theta.alter)
                         + 1 - pnorm((qnorm(1 - alpha / 2) * sigma.theta.null + mu.theta.null
                                       - mu.theta.alter) / sigma.theta.alter))
  # values at theta.0
  mu.theta.0        <- (log((d2 - theta.0 * log(c1.hat) * O2.stroke.star) / (o2.stroke + O2.stroke.star))
                         - log((d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star)))
  sigma.theta.0     <- (sqrt(1 / (n.star
                                   * (1 - exp(- (d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star))
                                           * (1 + ((d1 - log(c1.hat) * O1.stroke.star) / (o1.stroke + O1.stroke.star))^2)))
                              + 1 / (n.star
                                      * (1 - exp(- (d2 - theta.0 * log(c1.hat) * O2.stroke.star) / (o2.stroke + O2.stroke.star))
                                              * (1 + ((d2 - theta.0 * log(c1.hat) * O2.stroke.star) / (o2.stroke + O2.stroke.star))^2)))))
  gamma.theta.0     <- (pnorm((qnorm(alpha / 2) * sigma.theta.null + mu.theta.null
                                - mu.theta.0) / sigma.theta.0)
                         + 1 - pnorm((qnorm(1 - alpha / 2) * sigma.theta.null + mu.theta.null
                                       - mu.theta.0) / sigma.theta.0))

  return(list(theta, gamma.theta, gamma.theta.0))
}