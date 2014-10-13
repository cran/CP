DispDataExp <- function(group1.name, n1, d1, n1.alive, o1, lambda1.hat, O1.star,
                        group2.name, n2, d2, n2.alive, o2, lambda2.hat, O2.star,
                        theta.0, theta.hat) {

  ## Displays the passed parameters.
  ##
  ## Args:
  ##   Parameters from exponential power calculations.
  ##
  ## Results:
  ##   Displays the passed parameters.

  res1           <- data.frame(cbind(c(n1,
                                       d1,
                                       n1.alive,
                                       floor(x = o1),
                                       formatC(x      = lambda1.hat,
                                               digits = 4,
                                               format = "f"),
                                       floor(x = O1.star)),
                                     c(n2,
                                       d2,
                                       n2.alive,
                                      floor(x = o2),
                                      formatC(x      = lambda2.hat,
                                              digits = 4,
                                              format = "f"),
                                      floor(x = O2.star))),
                               row.names = c("Patients",
                                             "Death Events",
                                             "Censored",
                                             "Person Months",
                                             "Estimated Hazard",
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