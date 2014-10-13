SplitData <- function(data) {

  ## Splits data frame into two data frames, each for one group,
  ## and converts group expressions for internal calculations
  ## into values 1 and 2.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns with the group
  ##         (two different expressions) in the first,
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##
  ## Results:
  ##   Returns two data frames, each for one group,
  ##   and the original names of the two groups.

  data0         <- split(x = data,
                         f = data[, 1])
  data1         <- data0[[1]]
  data2         <- data0[[2]]
  group1.name <- as.character(x = data1[1, 1])
  group2.name <- as.character(x = data2[1, 1])
  data1[, 1]    <- as.numeric(x = 1)
  data2[, 1]    <- as.numeric(x = 2)

  return(list(data1, group1.name, data2, group2.name))
}