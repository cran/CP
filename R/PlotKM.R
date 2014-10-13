PlotKM <- function(data, model.name) {

  ## Plots the Kaplan-Meier curves for the passed data frame.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns with the group 
  ##         (values 1 and 2) in the first,
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   model.name: Name of the used model for estimation.
  ##
  ## Returns:
  ##   Plots the Kaplan-Meier curves
  ##   being the first component in an 1 x 2 array
  ##   in which the second component is still empty.

  # package 'survival' for survival analysis
  km <- survival::survfit(formula = survival::Surv(time  = data[, 3],
                                                   event = data[, 2]) ~ data[, 1])
  # plot of Kaplan-Meier curves
  plot(x    = km,
       main = model.name,
       xlab = "Time",
       ylab = "Survival",
       col  = c("blue", "green"),
       ylim = c(0, 1))
}