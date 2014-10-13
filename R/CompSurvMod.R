CompSurvMod <- function(data, cont.time,
                        new.pat = c(0, 0), theta.0 = 1, alpha = 0.05,
                        disp.data = FALSE, plot.km = FALSE) {

  ## Calculates the conditional power and plots the conditional power curve
  ## for the exponential model and the non-mixture models with
  ## exponential, Weibull type and Gamma type survival
  ## with respect to two different treatments and no dropouts.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns with the group
  ##         (two different expressions) in the first,
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   cont.time: Period of time of continuing the trial.
  ##   new.pat: 2-dimensional vector which consists of numbers of new patients
  ##            who will be recruited each time unit
  ##            (first component = group 1, second component = group 2)
  ##            with default at (0, 0).
  ##   theta.0: Originally postulated clinically relevant difference
  ##            (hazard ratio = hazard of group 2 / hazard of group 1)
  ##            with default at 1.
  ##   alpha: Significance level for conditional power calculations
  ##          with default at 0.05.
  ##   disp.data: Logical value indicating if all calculated data should be displayed
  ##              with default at FALSE.
  ##   plot.km: Logical value indicating if Kaplan-Meier curves
  ##            and estimated survival curves according to
  ##            the mentioned models should be plotted
  ##            with default at FALSE.
  ##
  ## Results:
  ##   Displays the calculated conditional power
  ##   and optionally an overview of the other calculated values,
  ##   and plots the conditional power curves
  ##   and optionally the Kaplan-Meier curves.
  ##   Returns the estimates of the parameters, the hazard ratios
  ##   and the conditional power.

  # check of passed parameters
  IsValid(data, cont.time, new.pat, theta.0, alpha, disp.data, plot.km)

  # split data frame into two data frames, each for one group,
  # and converting group expressions for internal calculations
  # into values 1 and 2
  split.data  <- SplitData(data)
  data1       <- split.data[[1]]
  group1.name <- split.data[[2]]
  data2       <- split.data[[3]]
  group2.name <- split.data[[4]]



  ########################
  # EXPONENTIAL SURVIVAL #
  ########################

  # maximum likelihood estimators for hazards of group 1 and group 2
  d1              <- sum(data1[, 2])
  o1              <- sum(data1[, 3])
  d2              <- sum(data2[, 2])
  o2              <- sum(data2[, 3])
  lambda1.hat.exp <- d1 / o1
  lambda2.hat.exp <- d2 / o2

  # maximum of log-likehood function without prefactor of censoring
  # and AIC = - 2 * log-likelihood + 2 * Parameter
  # of group 1 and group 2
  log.likelihood1.hat.exp <- sum(data1[, 2] * log(lambda1.hat.exp) - lambda1.hat.exp * data1[, 3])
  log.likelihood2.hat.exp <- sum(data2[, 2] * log(lambda2.hat.exp) - lambda2.hat.exp * data2[, 3])
  AIC1.exp                <- - 2 * log.likelihood1.hat.exp + 2 * 1
  AIC2.exp                <- - 2 * log.likelihood2.hat.exp + 2 * 1

  # estimator for hazard ratio theta = lambda2 / lambda1
  theta.hat.exp <- lambda2.hat.exp / lambda1.hat.exp

  # estimation of person months in group 1 and group 2
  n1.alive    <- sum(1 - data1[, 2])
  n2.alive    <- sum(1 - data2[, 2])
  O1.star.exp <- PersMonExp(d1, o1, n1.alive, new.pat[1], cont.time)
  O2.star.exp <- PersMonExp(d2, o2, n2.alive, new.pat[2], cont.time)

  # number of patients
  n.alive <- n1.alive + n2.alive
  rel     <- n.alive / length(x = data[, 1])
  n.star  <- floor(x = (n.alive + ((new.pat[1] + new.pat[2]) * cont.time * rel)))

  # conditional power calculations
  calc.conpwr.exp   <- CalcConPwrExp(theta.0,
                                     d1, o1, O1.star.exp, lambda1.hat.exp,
                                     d2, o2, O2.star.exp,
                                     n.star,
                                     alpha)
  theta             <- calc.conpwr.exp[[1]]
  gamma.theta.exp   <- calc.conpwr.exp[[2]]
  gamma.theta.0.exp <- calc.conpwr.exp[[3]]

  # summary vector
  summary.exp <- c(log.likelihood1.hat.exp, log.likelihood2.hat.exp,
                   AIC1.exp, AIC2.exp,
                   lambda1.hat.exp, lambda2.hat.exp,
                   O1.star.exp, O2.star.exp,
                   theta.hat.exp)



  #####################################
  # NON-MIXTURE: EXPONENTIAL SURVIVAL #
  #####################################

  # calculate initial values for maximum likelihood estimation
  # of parameters in group 1
  # and if applicable projection into feasible region
  init.val.data1.likelihood.nonmix.exp <- InitValLikelihoodNonMixExp(data1)
  # initial values for maximum likelihood estimation
  # of parameters in group 1
  lambda1.0                            <- init.val.data1.likelihood.nonmix.exp[1]
  c1.0                                 <- init.val.data1.likelihood.nonmix.exp[2]
  # calculate initial values for maximum likelihood estimation
  # of parameters in group 2
  # and if applicable projection into feasible region
  init.val.data2.likelihood.nonmix.exp <- InitValLikelihoodNonMixExp(data2)
  # initial values for maximum likelihood estimation
  # of parameters in group 2
  lambda2.0                            <- init.val.data2.likelihood.nonmix.exp[1]
  c2.0                                 <- init.val.data2.likelihood.nonmix.exp[2]
  # calculate initial values for maximum likelihood estimation
  # of parameters for all data
  # and if applicable projection into feasible region
  init.val.data.likelihood.nonmix.exp  <- InitValLikelihoodNonMixExp(data)
  # initial values for maximum likelihood estimation
  # of parameters for all data
  lambda.0                             <- init.val.data.likelihood.nonmix.exp[1]
  c.0                                  <- init.val.data.likelihood.nonmix.exp[2]
  # maximum likelihood estimation of parameters in group 1, group 2
  # and for all data
  likelihood.nonmix.exp                <- LikelihoodNonMixExp(data1, data2, data,
                                                              lambda1.0, c1.0,
                                                              lambda2.0, c2.0,
                                                              lambda.0, c.0)
  # maximum likelihood estimators of parameters,
  # maximum of log-likelihood function without prefactor of censoring
  # and AIC = - 2 * log-likelihood + 2 * Parameter
  # of group 1 and group 2
  lambda1.hat.nonmix.exp             <- likelihood.nonmix.exp[1]
  c1.hat.nonmix.exp                  <- likelihood.nonmix.exp[2]
  lambda2.hat.nonmix.exp             <- likelihood.nonmix.exp[3]
  c2.hat.nonmix.exp                  <- likelihood.nonmix.exp[4]
  lambda.hat.nonmix.exp              <- likelihood.nonmix.exp[5]
  c1.cond.hat.nonmix.exp             <- likelihood.nonmix.exp[6]
  c2.cond.hat.nonmix.exp             <- likelihood.nonmix.exp[7]
  log.likelihood1.hat.nonmix.exp     <- likelihood.nonmix.exp[8]
  AIC1.nonmix.exp                    <- - 2 * log.likelihood1.hat.nonmix.exp + 2 * 2
  log.likelihood2.hat.nonmix.exp     <- likelihood.nonmix.exp[9]
  AIC2.nonmix.exp                    <- - 2 * log.likelihood2.hat.nonmix.exp + 2 * 2

  # estimator for hazard ratio theta = log(c2) / log(c1)
  # under the assumption lambda1 = lambda2
  theta.hat.nonmix.exp <- log(c2.cond.hat.nonmix.exp) / log(c1.cond.hat.nonmix.exp)

  # estimation of person months in group 1 and group 2
  n1.alive           <- sum(1 - data1[, 2])
  n2.alive           <- sum(1 - data2[, 2])
  O1.star.nonmix.exp <- PersMonNonMixExp(lambda1.hat.nonmix.exp, c1.hat.nonmix.exp,
                                         n1.alive, new.pat[1], cont.time)
  O2.star.nonmix.exp <- PersMonNonMixExp(lambda2.hat.nonmix.exp, c2.hat.nonmix.exp,
                                         n2.alive, new.pat[2], cont.time)

  # functions of person months in group 1 , group 2
  # and in group 2 under the null hypothesis
  o1.stroke      <- FctPersMonNonMixExp(data1, lambda1.hat.nonmix.exp, group1.name)
  o2.stroke      <- FctPersMonNonMixExp(data2, lambda2.hat.nonmix.exp, group2.name)
  o2.stroke.null <- FctPersMonNonMixExp(data2, lambda1.hat.nonmix.exp, group2.name)

  # further functions of person months in group 1, group 2
  # and in group 2 under the null hypothesis
  n1                  <- length(x = data1[, 1])
  n2                  <- length(x = data2[, 1])
  O1.stroke.star      <- o1.stroke / n1 * (n1.alive + new.pat[1] * (cont.time + 1) / 2)
  O2.stroke.star      <- o2.stroke / n2 * (n2.alive + new.pat[2] * (cont.time + 1) / 2)
  O2.stroke.star.null <- o2.stroke.null / n2 * (n2.alive + new.pat[2] * (cont.time + 1) / 2)

  # number of patients
  n.alive <- n1.alive + n2.alive
  rel     <- n.alive / (n1 + n2)
  n.star  <- floor(x = (n.alive + ((new.pat[1] + new.pat[2]) * cont.time * rel)))

  # conditional power calculations
  d1                       <- sum(data1[, 2])
  d2                       <- sum(data2[, 2])
  calc.conpwr.nonmix       <- CalcConPwrNonMix(theta.0,
                                               d1, o1.stroke, O1.stroke.star, c1.cond.hat.nonmix.exp,
                                               d2, o2.stroke, O2.stroke.star, O2.stroke.star.null,
                                               n.star,
                                               alpha)
  theta                    <- calc.conpwr.nonmix[[1]]
  gamma.theta.nonmix.exp   <- calc.conpwr.nonmix[[2]]
  gamma.theta.0.nonmix.exp <- calc.conpwr.nonmix[[3]]

  # summary vector
  summary.nonmix.exp <- c(log.likelihood1.hat.nonmix.exp, log.likelihood2.hat.nonmix.exp,
                          AIC1.nonmix.exp, AIC2.nonmix.exp,
                          lambda.hat.nonmix.exp, c1.cond.hat.nonmix.exp,
                          lambda.hat.nonmix.exp, c2.cond.hat.nonmix.exp,
                          O1.star.nonmix.exp, O2.star.nonmix.exp,
                          theta.hat.nonmix.exp)



  ######################################
  # NON-MIXTURE: WEIBULL TYPE SURVIVAL #
  ######################################

  # calculate initial values for maximum likelihood estimation
  # of parameters in group 1
  # and if applicable projection into feasible region
  init.val.data1.likelihood.nonmix.wei <- InitValLikelihoodNonMixWei(data1)
  # initial values for maximum likelihood estimation
  # of parameters in group 1
  lambda1.0                            <- init.val.data1.likelihood.nonmix.wei[1]
  k1.0                                 <- init.val.data1.likelihood.nonmix.wei[2]
  c1.0                                 <- init.val.data1.likelihood.nonmix.wei[3]
  # calculate initial values for maximum likelihood estimation
  # of parameters in group 2
  # and if applicable projection into feasible region
  init.val.data2.likelihood.nonmix.wei <- InitValLikelihoodNonMixWei(data2)
  # initial values for maximum likelihood estimation
  # of parameters in group 2
  lambda2.0                            <- init.val.data2.likelihood.nonmix.wei[1]
  k2.0                                 <- init.val.data2.likelihood.nonmix.wei[2]
  c2.0                                 <- init.val.data2.likelihood.nonmix.wei[3]
  # calculate initial values for maximum likelihood estimation
  # of parameters for all data
  # and if applicable projection into feasible region
  init.val.data.likelihood.nonmix.wei  <- InitValLikelihoodNonMixWei(data)
  # initial values for maximum likelihood estimation
  # of parameters for all data
  lambda.0                             <- init.val.data.likelihood.nonmix.wei[1]
  k.0                                  <- init.val.data.likelihood.nonmix.wei[2]
  c.0                                  <- init.val.data.likelihood.nonmix.wei[3]
  # maximum likelihood estimation of parameters in group 1, group 2
  # and for all data
  likelihood.nonmix.wei                <- LikelihoodNonMixWei(data1, data2, data,
                                                              lambda1.0, k1.0, c1.0,
                                                              lambda2.0, k2.0, c2.0,
                                                              lambda.0, k.0, c.0)
  # maximum likelihood estimators of parameters,
  # maximum of log-likelihood function without prefactor of censoring
  # and AIC = - 2 * log-likelihood + 2 * Parameter
  # of group 1 and group 2
  lambda1.hat.nonmix.wei              <- likelihood.nonmix.wei[1]
  k1.hat.nonmix.wei                   <- likelihood.nonmix.wei[2]
  c1.hat.nonmix.wei                   <- likelihood.nonmix.wei[3]
  lambda2.hat.nonmix.wei              <- likelihood.nonmix.wei[4]
  k2.hat.nonmix.wei                   <- likelihood.nonmix.wei[5]
  c2.hat.nonmix.wei                   <- likelihood.nonmix.wei[6]
  lambda.hat.nonmix.wei               <- likelihood.nonmix.wei[7]
  k.hat.nonmix.wei                    <- likelihood.nonmix.wei[8]
  c1.cond.hat.nonmix.wei              <- likelihood.nonmix.wei[9]
  c2.cond.hat.nonmix.wei              <- likelihood.nonmix.wei[10]
  log.likelihood1.hat.nonmix.wei      <- likelihood.nonmix.wei[11]
  AIC1.nonmix.wei                     <- - 2 * log.likelihood1.hat.nonmix.wei + 2 * 2
  log.likelihood2.hat.nonmix.wei      <- likelihood.nonmix.wei[12]
  AIC2.nonmix.wei                     <- - 2 * log.likelihood2.hat.nonmix.wei + 2 * 2

  # estimator for hazard ratio theta = log(c2) / log(c1)
  # under the assumption lambda1 = lambda2 and k1 = k2
  theta.hat.nonmix.wei <- log(c2.cond.hat.nonmix.wei) / log(c1.cond.hat.nonmix.wei)

  # estimation of person months in group 1 and group 2
  n1.alive           <- sum(1 - data1[, 2])
  n2.alive           <- sum(1 - data2[, 2])
  O1.star.nonmix.wei <- PersMonNonMixWei(lambda1.hat.nonmix.wei, k1.hat.nonmix.wei, c1.hat.nonmix.wei,
                                         n1.alive, new.pat[1], cont.time)
  O2.star.nonmix.wei <- PersMonNonMixWei(lambda2.hat.nonmix.wei, k2.hat.nonmix.wei, c2.hat.nonmix.wei,
                                         n2.alive, new.pat[2], cont.time)

  # functions of person months in group 1 , group 2
  # and in group 2 under the null hypothesis
  o1.stroke      <- FctPersMonNonMixWei(data1, lambda1.hat.nonmix.wei, k1.hat.nonmix.wei, group1.name)
  o2.stroke      <- FctPersMonNonMixWei(data2, lambda2.hat.nonmix.wei, k2.hat.nonmix.wei, group2.name)
  o2.stroke.null <- FctPersMonNonMixWei(data2, lambda1.hat.nonmix.wei, k1.hat.nonmix.wei, group2.name)

  # further functions of person months in group 1, group 2
  # and in group 2 under the null hypothesis
  n1                  <- length(x = data1[, 1])
  n2                  <- length(x = data2[, 1])
  O1.stroke.star      <- o1.stroke / n1 * (n1.alive + new.pat[1] * (cont.time + 1) / 2)
  O2.stroke.star      <- o2.stroke / n2 * (n2.alive + new.pat[2] * (cont.time + 1) / 2)
  O2.stroke.star.null <- o2.stroke.null / n2 * (n2.alive + new.pat[2] * (cont.time + 1) / 2)

  # number of patients
  n.alive <- n1.alive + n2.alive
  rel     <- n.alive / (n1 + n2)
  n.star  <- floor(x = (n.alive + ((new.pat[1] + new.pat[2]) * cont.time * rel)))

  # conditional power calculations
  d1                       <- sum(data1[, 2])
  d2                       <- sum(data2[, 2])
  calc.conpwr.nonmix       <- CalcConPwrNonMix(theta.0,
                                               d1, o1.stroke, O1.stroke.star, c1.cond.hat.nonmix.wei,
                                               d2, o2.stroke, O2.stroke.star, O2.stroke.star.null,
                                               n.star,
                                               alpha)
  theta                    <- calc.conpwr.nonmix[[1]]
  gamma.theta.nonmix.wei   <- calc.conpwr.nonmix[[2]]
  gamma.theta.0.nonmix.wei <- calc.conpwr.nonmix[[3]]

  # summary vector
  summary.nonmix.wei <- c(log.likelihood1.hat.nonmix.wei, log.likelihood2.hat.nonmix.wei,
                          AIC1.nonmix.wei, AIC2.nonmix.wei,
                          lambda.hat.nonmix.wei, k.hat.nonmix.wei, c1.cond.hat.nonmix.wei,
                          lambda.hat.nonmix.wei, k.hat.nonmix.wei, c2.cond.hat.nonmix.wei,
                          O1.star.nonmix.wei, O2.star.nonmix.wei,
                          theta.hat.nonmix.wei)



  ####################################
  # NON-MIXTURE: GAMMA TYPE SURVIVAL #
  ####################################

  # calculate initial values for maximum likelihood estimation
  # of parameters in group 1
  # and if applicable projection into feasible region
  init.val.data1.likelihood.nonmix.gamma <- InitValLikelihoodNonMixGamma(data1)
  # initial values for maximum likelihood estimation
  # of parameters in group 1
  a1.0                                   <- init.val.data1.likelihood.nonmix.gamma[1]
  b1.0                                   <- init.val.data1.likelihood.nonmix.gamma[2]
  c1.0                                   <- init.val.data1.likelihood.nonmix.gamma[3]
  # calculate initial values for maximum likelihood estimation
  # of parameters in group 2
  # and if applicable projection into feasible region
  init.val.data2.likelihood.nonmix.gamma <- InitValLikelihoodNonMixGamma(data2)
  # initial values for maximum likelihood estimation
  # of parameters in group 2
  a2.0                                   <- init.val.data2.likelihood.nonmix.gamma[1]
  b2.0                                   <- init.val.data2.likelihood.nonmix.gamma[2]
  c2.0                                   <- init.val.data2.likelihood.nonmix.gamma[3]
  # calculate initial values for maximum likelihood estimation
  # of parameters for all data
  # and if applicable projection into feasible region
  init.val.data.likelihood.nonmix.gamma  <- InitValLikelihoodNonMixGamma(data)
  # initial values for maximum likelihood estimation
  # of parameters for all data
  a.0                                    <- init.val.data.likelihood.nonmix.gamma[1]
  b.0                                    <- init.val.data.likelihood.nonmix.gamma[2]
  c.0                                    <- init.val.data.likelihood.nonmix.gamma[3]
  # maximum likelihood estimation of parameters in group 1, group 2
  # and for all data
  likelihood.nonmix.gamma                <- LikelihoodNonMixGamma(data1, data2, data,
                                                                  a1.0, b1.0, c1.0,
                                                                  a2.0, b2.0, c2.0,
                                                                  a.0, b.0, c.0)
  # maximum likelihood estimators of parameters,
  # maximum of log-likelihood function without prefactor of censoring
  # and AIC = - 2 * log-likelihood + 2 * Parameter
  # of group 1 and group 2
  a1.hat.nonmix.gamma                  <- likelihood.nonmix.gamma[1]
  b1.hat.nonmix.gamma                  <- likelihood.nonmix.gamma[2]
  c1.hat.nonmix.gamma                  <- likelihood.nonmix.gamma[3]
  a2.hat.nonmix.gamma                  <- likelihood.nonmix.gamma[4]
  b2.hat.nonmix.gamma                  <- likelihood.nonmix.gamma[5]
  c2.hat.nonmix.gamma                  <- likelihood.nonmix.gamma[6]
  a.hat.nonmix.gamma                   <- likelihood.nonmix.gamma[7]
  b.hat.nonmix.gamma                   <- likelihood.nonmix.gamma[8]
  c1.cond.hat.nonmix.gamma             <- likelihood.nonmix.gamma[9]
  c2.cond.hat.nonmix.gamma             <- likelihood.nonmix.gamma[10]
  log.likelihood1.hat.nonmix.gamma     <- likelihood.nonmix.gamma[11]
  AIC1.nonmix.gamma                    <- - 2 * log.likelihood1.hat.nonmix.gamma + 2 * 2
  log.likelihood2.hat.nonmix.gamma     <- likelihood.nonmix.gamma[12]
  AIC2.nonmix.gamma                    <- - 2 * log.likelihood2.hat.nonmix.gamma + 2 * 2

  # estimator for hazard ratio theta = log(c2) / log(c1)
  # under the assumption a1 = a2 and b1 = b2
  theta.hat.nonmix.gamma <- log(c2.cond.hat.nonmix.gamma) / log(c1.cond.hat.nonmix.gamma)

  # estimation of person months in group 1 and group 2
  n1.alive             <- sum(1 - data1[, 2])
  n2.alive             <- sum(1 - data2[, 2])
  O1.star.nonmix.gamma <- PersMonNonMixGamma(a1.hat.nonmix.gamma, b1.hat.nonmix.gamma, c1.hat.nonmix.gamma,
                                             n1.alive, new.pat[1], cont.time)
  O2.star.nonmix.gamma <- PersMonNonMixGamma(a2.hat.nonmix.gamma, b2.hat.nonmix.gamma, c2.hat.nonmix.gamma,
                                             n2.alive, new.pat[2], cont.time)

  # functions of person months in group 1 , group 2
  # and in group 2 under the null hypothesis
  o1.stroke      <- FctPersMonNonMixGamma(data1, a1.hat.nonmix.gamma, b1.hat.nonmix.gamma, group1.name)
  o2.stroke      <- FctPersMonNonMixGamma(data2, a2.hat.nonmix.gamma, b2.hat.nonmix.gamma, group2.name)
  o2.stroke.null <- FctPersMonNonMixGamma(data2, a1.hat.nonmix.gamma, b1.hat.nonmix.gamma, group2.name)

  # further functions of person months in group 1, group 2
  # and in group 2 under the null hypothesis
  n1                  <- length(x = data1[, 1])
  n2                  <- length(x = data2[, 1])
  O1.stroke.star      <- o1.stroke / n1 * (n1.alive + new.pat[1] * (cont.time + 1) / 2)
  O2.stroke.star      <- o2.stroke / n2 * (n2.alive + new.pat[2] * (cont.time + 1) / 2)
  O2.stroke.star.null <- o2.stroke.null / n2 * (n2.alive + new.pat[2] * (cont.time + 1) / 2)

  # number of patients
  n.alive <- n1.alive + n2.alive
  rel     <- n.alive / (n1 + n2)
  n.star  <- floor(x = (n.alive + ((new.pat[1] + new.pat[2]) * cont.time * rel)))

  # conditional power calculations
  d1                         <- sum(data1[, 2])
  d2                         <- sum(data2[, 2])
  calc.conpwr.nonmix         <- CalcConPwrNonMix(theta.0,
                                                 d1, o1.stroke, O1.stroke.star, c1.cond.hat.nonmix.gamma,
                                                 d2, o2.stroke, O2.stroke.star, O2.stroke.star.null,
                                                 n.star,
                                                 alpha)
  theta                      <- calc.conpwr.nonmix[[1]]
  gamma.theta.nonmix.gamma   <- calc.conpwr.nonmix[[2]]
  gamma.theta.0.nonmix.gamma <- calc.conpwr.nonmix[[3]]

  # summary vector
  summary.nonmix.gamma <- c(log.likelihood1.hat.nonmix.gamma, log.likelihood2.hat.nonmix.gamma,
                            AIC1.nonmix.gamma, AIC2.nonmix.gamma,
                            a.hat.nonmix.gamma, b.hat.nonmix.gamma, c1.cond.hat.nonmix.gamma,
                            a.hat.nonmix.gamma, b.hat.nonmix.gamma, c2.cond.hat.nonmix.gamma,
                            O1.star.nonmix.gamma, O2.star.nonmix.gamma,
                            theta.hat.nonmix.gamma)



  # results
  # additional data (optional)
  if (disp.data == TRUE) {
    # calculate number of death events, person months, number of patients
    # and number of patients still alive of group1 and group 2
    interim.data1 <- InterimData(data1, group1.name)
    interim.data2 <- InterimData(data2, group2.name)

    DispDataAll(group1.name, group2.name,      # 2 x  1 elements
                interim.data1, interim.data2,  # 2 x  4 elements
                summary.exp,                   #      9 elements
                summary.nonmix.exp,            #     11 elements
                summary.nonmix.wei,            #     13 elements
                summary.nonmix.gamma,          #     13 elements
                theta.0)                       #      1 elements
  }
  # conditional power
  DispConPwrAll(gamma.theta.0.exp,
                gamma.theta.0.nonmix.exp,
                gamma.theta.0.nonmix.wei,
                gamma.theta.0.nonmix.gamma,
                group1.name, group2.name)
  # standardization of plot window
  par(las   = 1,
      mfrow = c(1, 1))
  # plots of Kaplan-Meier curves
  # and estimated survival curves
  # according to the four mentioned models (optional)
  if (plot.km == TRUE) {
    layout(mat = matrix(data  = c(1, 2, 3, 4, 5, 5),
                        nrow  = 3,
                        ncol  = 2,
                        byrow = TRUE))

    # exponential model
    PlotKM(data, "Exponential Model")
    PlotEstExp(data1, data2,
               lambda1.hat.exp, lambda2.hat.exp,
               group1.name, group2.name)

    # non-mixture model with exponential survival
    PlotKM(data, "Non-Mixture Model with Exponential Survival")
    PlotEstNonMixExp(data1, data2,
                     lambda1.hat.nonmix.exp, c1.hat.nonmix.exp,
                     lambda2.hat.nonmix.exp, c2.hat.nonmix.exp,
                     group1.name, group2.name)

    # non-mixture model with Weibull type survival
    PlotKM(data, "Non-Mixture Model with Weibull type Survival")
    PlotEstNonMixWei(data1, data2,
                     lambda1.hat.nonmix.wei, k1.hat.nonmix.wei, c1.hat.nonmix.wei,
                     lambda2.hat.nonmix.wei, k2.hat.nonmix.wei, c2.hat.nonmix.wei,
                     group1.name, group2.name)

    # non-mixture model with Gamma type survival
    PlotKM(data, "Non-Mixture Model with Gamma type Survival")
    PlotEstNonMixGamma(data1, data2,
                       a1.hat.nonmix.gamma, b1.hat.nonmix.gamma, c1.hat.nonmix.gamma,
                       a2.hat.nonmix.gamma, b2.hat.nonmix.gamma, c2.hat.nonmix.gamma,
                       group1.name, group2.name)
  }
  # plot of conditional power curves
  PlotConPwrAll(theta,
                gamma.theta.exp,
                gamma.theta.nonmix.exp,
                gamma.theta.nonmix.wei,
                gamma.theta.nonmix.gamma,
                theta.0,
                gamma.theta.0.exp,
                gamma.theta.0.nonmix.exp,
                gamma.theta.0.nonmix.wei,
                gamma.theta.0.nonmix.gamma,
                group1.name, group2.name)

  # return values
  return(value = invisible(x = list(lambda1.hat.exp         = lambda1.hat.exp,
                                    lambda2.hat.exp         = lambda2.hat.exp,
                                    theta.hat.exp           = theta.hat.exp,
                                    gamma.theta.0.exp       = gamma.theta.0.exp,
                                    lambda.hat.nm.exp       = lambda.hat.nonmix.exp,
                                    c1.hat.nm.exp           = c1.cond.hat.nonmix.exp,
                                    c2.hat.nm.exp           = c2.cond.hat.nonmix.exp,
                                    theta.hat.nm.exp        = theta.hat.nonmix.exp,
                                    gamma.theta.0.nm.exp    = gamma.theta.0.nonmix.exp,
                                    lambda.hat.nm.wei       = lambda.hat.nonmix.wei,
                                    k.hat.nm.wei            = k.hat.nonmix.wei,
                                    c1.hat.nm.wei           = c1.cond.hat.nonmix.wei,
                                    c2.hat.nm.wei           = c2.cond.hat.nonmix.wei,
                                    theta.hat.nm.wei        = theta.hat.nonmix.wei,
                                    gamma.theta.0.nm.wei    = gamma.theta.0.nonmix.wei,
                                    a.hat.nm.gamma          = a.hat.nonmix.gamma,
                                    b.hat.nm.gamma          = b.hat.nonmix.gamma,
                                    c1.hat.nm.gamma         = c1.cond.hat.nonmix.gamma,
                                    c2.hat.nm.gamma         = c2.cond.hat.nonmix.gamma,
                                    theta.hat.nm.gamma      = theta.hat.nonmix.gamma,
                                    gamma.theta.0.nm.gamma  = gamma.theta.0.nonmix.gamma)))
}