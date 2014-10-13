ConPwrNonMixWei <- function(data, cont.time,
                            new.pat = c(0, 0), theta.0 = 1, alpha = 0.05,
                            disp.data = FALSE, plot.km = FALSE) {

  ## Calculates the conditional power and plots the conditional power curve
  ## for the non-mixture model with Weibull type survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t^k)], lambda > 0, k > 0, 0 < c < 1, t >= 0,
  ## with respect to two different treatments and no drop outs.
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
  ##            the non-mixture model with Weibull type survival should be plotted
  ##            with default at FALSE.
  ##
  ## Returns:
  ##   Displays the calculated conditional power
  ##   and optionally an overview of the other calculated values,
  ##   and plots the conditional power curve
  ##   and optionally the Kaplan-Meier curves
  ##   plus the estimated survival curves.
  ##   Returns the estimates of the parameters, the hazard ratio
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
  # maximum likelihood estimators of parameters in group 1, group 2
  lambda1.hat                       <- likelihood.nonmix.wei[1]
  k1.hat                            <- likelihood.nonmix.wei[2]
  c1.hat                            <- likelihood.nonmix.wei[3]
  lambda2.hat                       <- likelihood.nonmix.wei[4]
  k2.hat                            <- likelihood.nonmix.wei[5]
  c2.hat                            <- likelihood.nonmix.wei[6]
  lambda.hat                        <- likelihood.nonmix.wei[7]
  k.hat                             <- likelihood.nonmix.wei[8]
  c1.cond.hat                       <- likelihood.nonmix.wei[9]
  c2.cond.hat                       <- likelihood.nonmix.wei[10]

  # estimator for hazard ratio theta = log(c2) / log(c1)
  # under the assumption lambda1 = lambda2 and k1 = k2
  theta.hat <- log(c2.cond.hat) / log(c1.cond.hat)

  # estimation of person months in group 1 and group 2
  n1.alive <- sum(1 - data1[, 2])
  n2.alive <- sum(1 - data2[, 2])
  O1.star  <- PersMonNonMixWei(lambda1.hat, k1.hat, c1.hat, n1.alive, new.pat[1], cont.time)
  O2.star  <- PersMonNonMixWei(lambda2.hat, k2.hat, c2.hat, n2.alive, new.pat[2], cont.time)

  # functions of person months in group 1 , group 2
  # and in group 2 under the null hypothesis
  o1.stroke      <- FctPersMonNonMixWei(data1, lambda1.hat, k1.hat, group1.name)
  o2.stroke      <- FctPersMonNonMixWei(data2, lambda2.hat, k2.hat, group2.name)
  o2.stroke.null <- FctPersMonNonMixWei(data2, lambda1.hat, k1.hat, group2.name)

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
  d1                 <- sum(data1[, 2])
  d2                 <- sum(data2[, 2])
  calc.conpwr.nonmix <- CalcConPwrNonMix(theta.0,
                                         d1, o1.stroke, O1.stroke.star, c1.cond.hat,
                                         d2, o2.stroke, O2.stroke.star, O2.stroke.star.null,
                                         n.star,
                                         alpha)
  theta              <- calc.conpwr.nonmix[[1]]
  gamma.theta        <- calc.conpwr.nonmix[[2]]
  gamma.theta.0      <- calc.conpwr.nonmix[[3]]

  # results
  # additional data (optional)
  if (disp.data == TRUE) {
    # calculate number of death events, person months, number of patients
    # and number of patients still alive of group1 and group 2
    interim.data1 <- InterimData(data1, group1.name)
    interim.data2 <- InterimData(data2, group2.name)
    d1            <- interim.data1[1]
    o1            <- interim.data1[2]
    n1            <- interim.data1[3]
    n1.alive      <- interim.data1[4]
    d2            <- interim.data2[1]
    o2            <- interim.data2[2]
    n2            <- interim.data2[3]
    n2.alive      <- interim.data2[4]

    DispDataNonMixWei(group1.name, n1, d1, n1.alive, o1, lambda.hat, k.hat, c1.cond.hat, O1.star,
                      group2.name, n2, d2, n2.alive, o2, lambda.hat, k.hat, c2.cond.hat, O2.star,
                      theta.0, theta.hat)
  }
  # conditional power
  DispConPwr(gamma.theta.0, group1.name, group2.name)
  # standardization of plot window
  par(las   = 1,
      mfrow = c(1, 1))
  # plots of Kaplan-Meier curves (optional)
  if (plot.km == TRUE) {
    par(mfrow = c(1, 2))

    PlotKM(data, "Non-Mixture Model with Weibull type Survival")
    PlotEstNonMixWei(data1, data2,
                     lambda1.hat, k1.hat, c1.hat,
                     lambda2.hat, k2.hat, c2.hat,
                     group1.name, group2.name)
  }
  # plot of conditional power curve
  PlotConPwr(theta, gamma.theta,
             theta.0, gamma.theta.0,
             group1.name, group2.name,
             "Non-Mixture Model with Weibull type Survival")

  # return values
  return(value = invisible(x = list(lambda.hat    = lambda.hat,  k.hat  = k.hat,
                                    c1.hat        = c1.cond.hat, c2.hat = c2.cond.hat,
                                    theta.hat     = theta.hat,
                                    gamma.theta.0 = gamma.theta.0)))
}