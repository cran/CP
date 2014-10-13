PlotConPwr <- function(theta, gamma.theta,
                       theta.0, gamma.theta.0,
                       group1.name, group2.name,
                       model.name) {

  ## Plots the conditional power curve.
  ##
  ## Args:
  ##   theta: Vector of hazard ratios for plotting.
  ##   gamma.theta: Conditional power according to theta.
  ##   theta.0: Originally postulated clinically relevant difference
  ##            (hazard ratio = hazard of group 2 / hazard of group 1).
  ##   gamma.theta.0: Conditional power according to theta.0.
  ##   group1.name: Name of group 1.
  ##   group2.name: Name of group 2.
  ##   model.name: Name of the used model for estimation.
  ##
  ## Returns:
  ##   Plot of the condtional power curve.

  plot(x    = log(theta),
       y    = gamma.theta,
       type = "l",
       main = model.name,
       #sub  = paste("log(",
       #             theta.0,
       #             ") = ",
       #             formatC(x     = log(theta.0),
       #                    digits = 4,
       #                    format = "f"),
       #             "          ",
       #             "CP(",
       #             theta.0,
       #             ") = ",
       #             formatC(x      = gamma.theta.0,
       #                     digits = 4,
       #                     format = "f"),
       #             sep = ""),
       xlab = paste("log(Hazard Ratio) = log(Hazard ",
                    group2.name,
                    " / Hazard ",
                    group1.name,
                    ")",
                    sep = ""),
       ylab = "Conditional Power",
       col  = "red",
       ylim = c(0, 1))

  abline(v   = log(theta.0),
         lty = 3)

  abline(h   = gamma.theta.0,
         lty = 3)
}