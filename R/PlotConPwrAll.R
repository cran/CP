PlotConPwrAll <- function(theta,
                          gamma.theta.exp,
                          gamma.theta.nonmix.exp,
                          gamma.theta.nonmix.wei,
                          gamma.theta.nonmix.gamma,
                          theta.0,
                          gamma.theta.0.exp,
                          gamma.theta.0.nonmix.exp,
                          gamma.theta.0.nonmix.wei,
                          gamma.theta.0.nonmix.gamma,
                          group1.name, group2.name) {

  ## Plots the conditional power curves.
  ##
  ## Args:
  ##   theta: Vector of hazard ratios for plotting.
  ##   gamma.theta.exp: Conditional power within the exponential model.
  ##   gamma.theta.nonmix.exp: Conditional power within the non-mixture model
  ##                           with exponential survival.
  ##   gamma.theta.nonmix.wei: Conditional power within the non-mixture model
  ##                           with Weibull type survival.
  ##   gamma.theta.nonmix.gamma: Conditional power within the non-mixture model
  ##                             with Gamma type survival.
  ##   theta.0: Originally postulated clinically relevant difference
  ##            (hazard ratio = hazard of group 2 / hazard of group 1).
  ##   gamma.theta.0.exp: Conditional power within the exponential model
  ##                      according to theta.0.
  ##   gamma.theta.0.nonmix.exp: Conditional power within the non-mixture model
  ##                             with exponential survival
  ##                             according to theta.0.
  ##   gamma.theta.0.nonmix.wei: Conditional power within the non-mixture model
  ##                             with Weibull type survival
  ##                             according to theta.0.
  ##   gamma.theta.0.nonmix.gamma: Conditional power within the non-mixture model
  ##                               with Gamma type survival
  ##                               according to theta.0.
  ##   group1.name: Name of group 1.
  ##   group2.name: Name of group 2.
  ##
  ## Returns:
  ##   Plot of the condtional power curves.

  graphics::plot(x    = log(theta),
                 y    = gamma.theta.exp,
                 type = "l",
                 xlab = paste("log(Hazard Ratio) = log(Hazard ",
                              group2.name,
                              " / Hazard ",
                              group1.name,
                              ")",
                              sep = ""),
                 ylab = "Conditional Power",
                 col  = "red",
                 ylim = c(0, 1))

  graphics::lines(x   = log(theta),
                  y   = gamma.theta.nonmix.exp,
                  col = "blue")

  graphics::lines(x   = log(theta),
                  y   = gamma.theta.nonmix.wei,
                  col = "green")

  graphics::lines(x   = log(theta),
                  y   = gamma.theta.nonmix.gamma,
                  col = "yellow")

  graphics::abline(v   = log(theta.0),
                   lty = 3)

  graphics::legend(x      = "topright",
                   legend = c("Exponential",
                              "Non-Mixture-Exponential",
                              "Non-Mixture-Weibull",
                              "Non-Mixture-Gamma"),
                   col    = c("red",
                              "blue",
                              "green",
                              "yellow"),
                   lty    = c(1,
                              1,
                              1,
                              1),
                   bg     = "white")
}