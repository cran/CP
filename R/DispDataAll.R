DispDataAll <- function(group1.name, group2.name,
                        interim.data1, interim.data2,
                        summary.exp,
                        summary.nonmix.exp,
                        summary.nonmix.wei,
                        summary.nonmix.gamma,
                        theta.0) {

  ## Displays the passed parameters.
  ##
  ## Args:
  ##   Parameters from power calculations.
  ##
  ## Results:
  ##   Displays the passed parameters.

  # interim analysis
  res1           <- data.frame(cbind(c(interim.data1[3],
                                       interim.data1[1],
                                       interim.data1[4],
                                       floor(x = interim.data1[2])),
                                     c(interim.data2[3],
                                       interim.data2[1],
                                       interim.data2[4],
                                       floor(x = interim.data2[2]))),
                               row.names = c("Patients",
                                             "Death Events",
                                             "Censored",
                                             "Person Months"))
  colnames(res1) <- c(group1.name, group2.name)

  # exponential model
  res2           <- data.frame(cbind(c(formatC(x      = summary.exp[1],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.exp[3],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.exp[5],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.exp[7])),
                                     c(formatC(x      = summary.exp[2],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.exp[4],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.exp[6],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.exp[8]))),
                               row.names = c("log(Likelihood)",
                                             "AIC",
                                             "lambda",
                                             "Further Person Months"))
  colnames(res2) <- c(group1.name, group2.name)

  # non-mixture model with exponential survival
  res3           <- data.frame(cbind(c(formatC(x      = summary.nonmix.exp[1],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.exp[3],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.exp[5],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.exp[6],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.nonmix.exp[9])),
                                     c(formatC(x      = summary.nonmix.exp[2],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.exp[4],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.exp[7],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.exp[8],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.nonmix.exp[10]))),
                               row.names = c("log(Likelihood)",
                                             "AIC",
                                             "lambda",
                                             "c",
                                             "Further Person Months"))
  colnames(res3) <- c(group1.name, group2.name)

  # non-mixture model with Weibull type survival
  res4           <- data.frame(cbind(c(formatC(x      = summary.nonmix.wei[1],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[3],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[5],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[6],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[7],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.nonmix.wei[11])),
                                     c(formatC(x      = summary.nonmix.wei[2],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[4],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[8],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[9],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.wei[10],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.nonmix.wei[12]))),
                               row.names = c("log(Likelihood)",
                                             "AIC",
                                             "lambda",
                                             "k",
                                             "c",
                                             "Further Person Months"))
  colnames(res4) <- c(group1.name, group2.name)

  # non-mixture model with Gamma type survival
  res5           <- data.frame(cbind(c(formatC(x      = summary.nonmix.gamma[1],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[3],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[5],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[6],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[7],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.nonmix.gamma[11])),
                                     c(formatC(x      = summary.nonmix.gamma[2],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[4],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[8],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[9],
                                               digits = 4,
                                               format = "f"),
                                       formatC(x      = summary.nonmix.gamma[10],
                                               digits = 4,
                                               format = "f"),
                                       floor(x = summary.nonmix.gamma[12]))),
                               row.names = c("log(Likelihood)",
                                             "AIC",
                                             "a",
                                             "b",
                                             "c",
                                             "Further Person Months"))
  colnames(res5) <- c(group1.name, group2.name)

  cat("Interim Analysis", "\n")
  cat("----------------", "\n")
  print(x = res1)
  cat("\n")
  cat("Postulated Hazard Ratio:", theta.0, "\n\n\n\n")

  cat("Exponential", "\n")
  cat("AIC =",
      formatC(x      = (summary.exp[3] + summary.exp[4]),
              digits = 4,
              format = "f"),
      "\n")
  cat("-----------", "\n")
  print(x = res2)
  cat("\n")
  cat("Estimated Hazard Ratio:",
      formatC(x      = summary.exp[9],
              digits = 4,
              format = "f"),
      "\n\n\n\n")

  cat("Non-Mixture-Exponential", "\n")
  cat("AIC =",
      formatC(x      = (summary.nonmix.exp[3] + summary.nonmix.exp[4]),
              digits = 4,
              format = "f"),
      "\n")
  cat("-----------------------", "\n")
  print(x = res3)
  cat("\n")
  cat("Estimated Hazard Ratio:",
      formatC(x      = summary.nonmix.exp[11],
              digits = 4,
              format = "f"),
      "\n\n\n\n")

  cat("Non-Mixture-Weibull", "\n")
  cat("AIC =",
      formatC(x      = (summary.nonmix.wei[3] + summary.nonmix.wei[4]),
              digits = 4,
              format = "f"),
      "\n")
  cat("-------------------", "\n")
  print(x = res4)
  cat("\n")
  cat("Estimated Hazard Ratio:",
      formatC(x      = summary.nonmix.wei[13],
              digits = 4,
              format = "f"),
      "\n\n\n\n")

  cat("Non-Mixture-Gamma", "\n")
  cat("AIC =",
      formatC(x      = (summary.nonmix.gamma[3] + summary.nonmix.gamma[4]),
              digits = 4,
              format = "f"),
      "\n")
  cat("-----------------", "\n")
  print(x = res5)
  cat("\n")
  cat("Estimated Hazard Ratio:",
      formatC(x      = summary.nonmix.gamma[13],
              digits = 4,
              format = "f"),
      "\n\n\n\n")
}