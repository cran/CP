InterimData <- function(data, group.name) {

  ## Calculates number of death events, person months, number of patients
  ## and number of patients still alive out of the passed data frame.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns for the group
  ##         in the first (all of the same group),
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   group.name: Name of the group.
  ##
  ## Results:
  ##   Returns number of death events, person months, number of patients
  ##   and number of patients still alive.

  # number of death events, and its verification
  d <- sum(data[, 2])
  if (d < 0) {
    stop(paste("Number of death events in", group.name, "must be bigger than or equal to 0.",
         call. = FALSE))
  }

  # person months, and its verification
  o <- sum(data[, 3])
  if (o <= 0) {
    stop(paste("Number of person months in", group.name, "must be bigger than 0.",
         call. = FALSE))
  }

  # number of patients at the beginning
  n <- length(data[, 1])

  # number of patients still alive
  n.alive <- sum(1 - data[, 2])

  return(c(d, o, n, n.alive))
}