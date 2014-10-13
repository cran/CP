\name{CompSurvMod}

\alias{CompSurvMod}

\title{
 Conditional Power (Comparison)
}

\description{
 Calculates the conditional power within the exponential model
 and the non-mixture models with exponential, Weibull type
 and Gamma type survival.
}

\usage{
 CompSurvMod(data, cont.time, new.pat = c(0, 0),
             theta.0 = 1, alpha = 0.05,
             disp.data = FALSE, plot.km = FALSE)
}

\arguments{
 \item{data}{
  Data frame which consists of at least three columns with the group
  (two different expressions) in the first,
  status (1 = event, 0 = censored) in the second
  and event time in the third column.
 }
 \item{cont.time}{
  Period of time of continuing the trial.
 }
 \item{new.pat}{
  2-dimensional vector which consists of numbers of new patients
  who will be recruited each time unit
  (first component = group 1, second component = group 2)
  with default at (0, 0).
 }
 \item{theta.0}{
  Originally postulated clinically relevant difference
  (hazard ratio = hazard of group 2 / hazard of group 1)
  with default at 1.
 }
 \item{alpha}{
  Significance level for conditional power calculations
  with default at 0.05.
 }
 \item{disp.data}{
  Logical value indicating if all calculated data should be displayed
  with default at FALSE.
 }
 \item{plot.km}{
  Logical value indicating if Kaplan-Meier curves
  and estimated survival curves according to
  the four mentioned models should be plotted
  with default at FALSE.
 }
}

\details{
 This function calculates the conditional power within the exponential model
 and the non-mixture models with exponential, Weibull type
 and Gamma type survival and plots the conditional power curves.
 
 Optionally, further data will be displayed.
 This includes data from interim analysis, log-likelihoods, AICs,
 calculated estimators and further patient times.
 
 Moreover, it is possible to plot the Kaplan-Meier curves
 and the estimated survival curves according to
 the four mentioned models.
}

\value{
 See Details.
 
 Returns a list which consists of the following components:
 \item{lambda1.hat.exp}{
  estimated hazard of group 1
  within the exponential model
 }
 \item{lambda2.hat.exp}{
  estimated hazard of group 2
  within the exponential model
 }
 \item{theta.hat.exp}{
  estimated hazard ratio
   = estimated hazard of group 2 / estimated hazard of group 1
  within the exponential model
 }
 \item{gamma.theta.0.exp}{
  conditional power
  within the exponential model
 }
 \item{lambda1.hat.nm.exp}{
  estimated rate parameter of group 1
  within the non-mixture model with exponential survival
 }
 \item{c1.hat.nm.exp}{
  estimated survival fraction of group 1
  within the non-mixture model with exponential survival
 }
 \item{lambda2.hat.nm.exp}{
  estimated rate parameter of group 2
  within the non-mixture model with exponential survival
 }
 \item{c2.hat.nm.exp}{
  estimated survival fraction of group 2
  within the non-mixture model with exponential survival
 }
 \item{theta.hat.nm.exp}{
  estimated hazard ratio
   = \eqn{log(}estimated survival fraction of group 2\eqn{)}
      / \eqn{log(}estimated survival fraction of group 1\eqn{)}
  within the non-mixture model with exponential survival
 }
 \item{gamma.theta.0.nm.exp}{
  conditional power
  within the non-mixture model with exponential survival
 }
 \item{lambda1.hat.nm.wei}{
  estimated scale parameter of group 1
  within the non-mixture model with Weibull type survival
 }
 \item{k1.hat.nm.wei}{
  estimated shape parameter of group 1
  within the non-mixture model with Weibull type survival
 }
 \item{c1.hat.nm.wei}{
  estimated survival fraction of group 1
  within the non-mixture model with Weibull type survival
 }
 \item{lambda2.hat.nm.wei}{
  estimated scale parameter of group 2
  within the non-mixture model with Weibull type survival
 }
 \item{k2.hat.nm.wei}{
  estimated shape parameter of group 2
  within the non-mixture model with Weibull type survival
 }
 \item{c2.hat.nm.wei}{
  estimated survival fraction of group 2
  within the non-mixture model with Weibull type survival
 }
 \item{theta.hat.nm.wei}{
  estimated hazard ratio
   = \eqn{log(}estimated survival fraction of group 2\eqn{)}
      / \eqn{log(}estimated survival fraction of group 1\eqn{)}
  within the non-mixture model with Weibull type survival
 }
 \item{gamma.theta.0.nm.wei}{
  conditional power
  within the non-mixture model with Weibull type survival
 }
 \item{a1.hat.nm.gamma}{
  estimated shape parameter of group 1
  within the non-mixture model with Gamma type survival
 }
 \item{b1.hat.nm.gamma}{
  estimated rate parameter of group 1
  within the non-mixture model with Gamma type survival
 }
 \item{c1.hat.nm.gamma}{
  estimated survival fraction of group 1
  within the non-mixture model with Gamma type survival
 }
 \item{a2.hat.nm.gamma}{
  estimated shape parameter of group 2
  within the non-mixture model with Gamma type survival
 }
 \item{b2.hat.nm.gamma}{
  estimated rate parameter of group 2
  within the non-mixture model with Gamma type survival
 }
 \item{c2.hat.nm.gamma}{
  estimated survival fraction of group 2
  within the non-mixture model with Gamma type survival
 }
 \item{theta.hat.nm.gamma}{
  estimated hazard ratio
   = \eqn{log(}estimated survival fraction of group 2\eqn{)}
      / \eqn{log(}estimated survival fraction of group 1\eqn{)}
  within the non-mixture model with Gamma type survival
 }
 \item{gamma.theta.0.nm.gamma}{
  conditional power
  within the non-mixture model with Gamma type survival
 }
}

\note{
 There are several mechanisms to ensure that no illegal operations will be done
 and maximum likelihood calculations will be executed stable.
 That is why there should not be too less data in the data frame,
 for example one patient of each group and both being censored.
}

\author{
 Andreas Kuehnapfel
}

\references{
 Kuehnapfel, A. (2013).
 Die bedingte Power in der Ueberlebenszeitanalyse.
}

\seealso{
 \code{\link{CP}} \cr
 \code{\link{ConPwrExp}} \cr
 \code{\link{ConPwrNonMixExp}} \cr
 \code{\link{ConPwrNonMixWei}} \cr
 \code{\link{ConPwrNonMixGamma}} \cr
 \code{\link{ConPwrExpAndersen}} \cr
 \code{\link{GenerateDataFrame}} \cr
 \code{\link{test}}
}

\examples{
 # data frame 'test' generated by 'GenerateDataFrame'
 
 # conditional power calculations
 # within the four mentioned models
 CompSurvMod(data = test, cont.time = 12, new.pat = c(2.5, 2.5),
             theta.0 = 0.75, alpha = 0.05,
             disp.data = TRUE, plot.km = TRUE)
}