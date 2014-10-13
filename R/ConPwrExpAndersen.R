ConPwrExpAndersen <- function(data, cont.time,
                              new.pat = c(0, 0), theta.0 = 1, alpha = 0.05,
                              disp.data = FALSE, plot.km = FALSE) {

  ## Calculates the conditional power and plots the conditional power curve
  ## for the exponential model with constant hazards by Per Kragh Andersen, i. e.
  ##  S(t) = exp(- lambda * t), lambda > 0, t >= 0,
  ## with respect to two different treatments and no drop outs.
  ## The original formulae of the Andersen paper are used.
  ## (Andersen, P. K. (1987). Conditional power calculations as an aid
  ##  in the decision whether to continue a clinical trial.
  ##  Controlled Clinical Trials 8, 67-74.)
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
  ##            the exponential model should be plotted
  ##            with default at FALSE.
  ##
  ## Returns:
  ##   Displays the calculated conditional power
  ##   and optionally an overview of the other calculated values,
  ##   and plots the conditional power curve
  ##   and optionally the Kaplan-Meier curves
  ##   plus the estimated survival curves.
  ##   Returns the estimates of the hazards, the hazard ratio
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

  # maximum likelihood estimators for hazards of group 1 and group 2
  d1          <- sum(data1[, 2])
  o1          <- sum(data1[, 3])
  d2          <- sum(data2[, 2])
  o2          <- sum(data2[, 3])
  lambda1.hat <- d1 / o1
  lambda2.hat <- d2 / o2

  # estimator for hazard ratio theta = lambda2 / lambda1
  theta.hat <- lambda2.hat / lambda1.hat

  # estimation of patient times in group 1 and group 2
  n1.alive <- sum(1 - data1[, 2])
  n2.alive <- sum(1 - data2[, 2])
  O1.star  <- PersMonExp(d1, o1, n1.alive, new.pat[1], cont.time)
  O2.star  <- PersMonExp(d2, o2, n2.alive, new.pat[2], cont.time)

  # conditional power calculations
  # with the original formulae of the Andersen paper
  calc.conpwr.exp.andersen <- CalcConPwrExpAndersen(theta.0,
                                                    d1, o1, O1.star, lambda1.hat,
                                                    d2, o2, O2.star,
                                                    alpha)
  theta                    <- calc.conpwr.exp.andersen[[1]]
  gamma.theta              <- calc.conpwr.exp.andersen[[2]]
  gamma.theta.0            <- calc.conpwr.exp.andersen[[3]]



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
     
    DispDataExp(group1.name, n1, d1, n1.alive, o1, lambda1.hat, O1.star,
                group2.name, n2, d2, n2.alive, o2, lambda2.hat, O2.star,
                theta.0, theta.hat)
  }
  # conditional power
  DispConPwr(gamma.theta.0, group1.name, group2.name)
  # standardization of plot window
  par(las   = 1,
      mfrow = c(1, 1))
  # plot of Kaplan-Meier curves and estimated survival (optional)
  if (plot.km == TRUE) {
    par(mfrow = c(1, 2))

    PlotKM(data, "Exponential Model (Andersen)")
    PlotEstExp(data1, data2,
               lambda1.hat, lambda2.hat,
               group1.name, group2.name)
  }
  # plot of conditional power curve
  PlotConPwr(theta, gamma.theta,
             theta.0, gamma.theta.0,
             group1.name, group2.name,
             "Exponential Model (Andersen)")

  # return values
  return(value = invisible(x = list(lambda1.hat   = lambda1.hat,
                                    lambda2.hat   = lambda2.hat,
                                    theta.hat     = theta.hat,
                                    gamma.theta.0 = gamma.theta.0)))
}