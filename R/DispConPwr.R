DispConPwr <- function(gamma.theta.0, group1.name, group2.name) {

  ## Prints the calculated conditional power.
  ##
  ## Args:
  ##   gamma.theta.0: Conditional power.
  ##   group1.name: Name of group 1.
  ##   group2.name: Name of group 2.
  ##
  ## Results:
  ##   Returns the calculated conditional power.

  cat("Conditional Power",
      formatC(x      = gamma.theta.0,
              digits = 4,
              format = "f"),
      "\n\n")

  message(paste("Note: Conditional power is calculated in view of the hazard ratio which is defined as the ratio of the hazard of group ",
                group2.name,
                " to the hazard of group ",
                group1.name,
                ".",
                sep = ""))
}