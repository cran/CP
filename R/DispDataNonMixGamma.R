DispDataNonMixGamma <- function(group1.name, n1, d1, n1.alive, o1, a1.hat, b1.hat, c1.hat, O1.star,
                                group2.name, n2, d2, n2.alive, o2, a2.hat, b2.hat, c2.hat, O2.star,
                                theta.0, theta.hat) {

  ## Displays the passed parameters.
  ##
  ## Args:
  ##   Parameters from non-mixture
  ##   Gamma type power calculations.
  ##
  ## Results:
  ##   Displays the passed parameters.

  res1           <- data.frame(cbind(c(n1,
                                       d1,
                                       n1.alive,
                                       floor(x = o1),
                                       formatC(x      = a1.hat,
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = b1.hat,
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = c1.hat,
                                               digits = 4,
                                               format = "f"),
                                       floor(x = O1.star)),
                                     c(n2,
                                       d2,
                                       n2.alive,
                                       floor(x = o2),
                                       formatC(x      = a2.hat,
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = b2.hat,
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = c2.hat,
                                               digits = 4,
                                               format = "f"),
                                       floor(x = O2.star))),
                               row.names = c("Patients",
                                             "Death Events",
                                             "Censored",
                                             "Person Months",
                                             "Estimated Shape Parameter",
                                             "Estimated Rate Parameter",
                                             "Estimated Survival Fraction",
                                             "Further Person Months"))
  colnames(res1) <- c(group1.name, group2.name)

  res2           <- data.frame(c(formatC(x      = theta.0,
                                         digits = 4,
                                         format = "f"),
                                 formatC(x      = theta.hat,
                                         digits = 4,
                                         format = "f")),
                               row.names = c("Postulated Hazard Ratio",
                                             "Estimated Hazard Ratio"))
  colnames(res2) <- ""

  print(x = res1)
  print(x = res2)
  cat("\n")
}