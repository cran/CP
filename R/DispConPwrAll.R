DispConPwrAll <- function(gamma.theta.0.exp,
                          gamma.theta.0.nonmix.exp,
                          gamma.theta.0.nonmix.wei,
                          gamma.theta.0.nonmix.gamma,
                          group1.name, group2.name) {

  ## Prints the calculated conditional power.
  ##
  ## Args:
  ##   gamma.theta.0.exp: Conditional power within the exponential model.
  ##   gamma.theta.0.nonmix.exp: Conditional power within the non-mixture model
  ##                             with exponential survival.
  ##   gamma.theta.0.nonmix.wei: Conditional power within the non-mixture model
  ##                             with Weibull type survival.
  ##   gamma.theta.0.nonmix.gamma: Conditional power within the non-mixture model
  ##                               with Gamma type survival.
  ##   group1.name: Name of group 1.
  ##   group2.name: Name of group 2.
  ##
  ## Results:
  ##   Returns the calculated conditional power.

  res               <- data.frame(c(formatC(x      = gamma.theta.0.exp,
                                            digits = 4,
                                            format = "f"),
                                    formatC(x      = gamma.theta.0.nonmix.exp,
                                            digits = 4,
                                            format = "f"),
                                    formatC(x      = gamma.theta.0.nonmix.wei,
                                            digits = 4,
                                            format = "f"),
                                    formatC(x      = gamma.theta.0.nonmix.gamma,
                                            digits = 4,
                                            format = "f")),
                                  row.names = c("Exponential",
                                                "Non-Mixture-Exponential",
                                                "Non-Mixture-Weibull",
                                                "Non-Mixture-Gamma"))
  colnames(x = res) <- "Conditional Power"

  cat("Conditional Power", "\n")
  cat("-----------------", "\n")
  print(x = res)
  cat("\n")
  
  message(paste("Note: Conditional power is calculated in view of the hazard ratio which is defined as the ratio of the hazard of group ",
                group2.name,
                " to the hazard of group ",
                group1.name,
                ".",
                sep = ""))
}