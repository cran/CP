GenerateDataFrame <- function() {

  ## Generates a data frame for testing conditional power calculations
  ## where the data is created by random.
  ## The data frame consists of three columns:
  ##   First = group (two different expressions: 'A' and 'B')
  ##   Second = status (1 = event, 0 = censored)
  ##   Third = event time
  ##
  ## Args:
  ##   None
  ##
  ## Results:
  ##   Returns data frame for testing conditional power calculations.

  # number of patients (circa 200)
  n <- rpois(n      = 1,
             lambda = 200)

  # initialization of group, status and event time vectors
  group <- vector(mode   = "character",
                  length = 0)
  stat  <- vector(mode   = "numeric",
                  length = 0)
  time  <- vector(mode   = "numeric",
                  length = 0)

  # probability of being censored for each patient of group A
  # respectively group B between 0.4 and 0.6
  # to avoid extreme (implausible) values
  cens.prob.A <- runif(n   = 1,
                       min = 0.4,
                       max = 0.6)
  cens.prob.B <- runif(n   = 1,
                       min = 0.4,
                       max = 0.6)

  # generating data rowwise
  # by choice of group with probability 0.5
  for (i in 1 : n) {
      if (rbinom(n    = 1,
                 size = 1,
                 prob = 0.5) == 1) {
        group <- c(group, "A")
        stat  <- c(stat, rbinom(n    = 1,
                                size = 1,
                                prob = (1 - cens.prob.A)))
        time  <- c(time, (0.01 + round(x      = rexp(n    = 1,
                                                     rate = (1 - cens.prob.A)),
                                       digits = 2)))
      }
      else {
        group <- c(group, "B")
        stat  <- c(stat, rbinom(n    = 1,
                                size = 1,
                                prob = (1 - cens.prob.B)))
        time  <- c(time, (0.01 + round(x      = rexp(n    = 1,
                                                     rate = (1 - cens.prob.B)),
                                       digits = 2)))
      }
  }

  # data frame
  data <- data.frame(group, stat, time)

  return(data)
}