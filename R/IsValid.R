IsValid <- function(data, cont.time, new.pat, theta.0, alpha, disp.data, plot.km) {

  ## Checks the passed parameters.
  ##
  ## Args:
  ##   data: Data frame which consists of at least three columns with the group
  ##         (two different expressions) in the first,
  ##         status (1 = event, 0 = censored) in the second
  ##         and event time in the third column.
  ##   cont.time: Period of time of continuing the trial.
  ##   new.pat: 2-dimensional vector which consists of numbers of new patients
  ##            who will be recruited each time unit
  ##            (first component = group 1, second component = group 2).
  ##   theta.0: Originally postulated clinically relevant difference (hazard ratio).
  ##   alpha: Significance level for conditional power calculations.
  ##   disp.data: Logical value indicating if all calculated data should be displayed.
  ##   plot.km: Logical value indicating if Kaplan-Meier curves should be plotted.
  ##
  ## Returns:
  ##   Stops calculations if parameters are invalid, otherwise no abort.

  # check of data
  if (length(x = data) < 3) {
    stop("data must have three or more columns.",
         call. = FALSE)
  }
  if (length(x = data[, 1]) <= 0) {
    stop("Columns of data must have dimension bigger than 0.",
         call. = FALSE)
  }
  if (length(x = levels(x = as.factor(x = data[, 1]))) != 2) {
    stop("First column of data must consist of exactly two different groups.",
         call. = FALSE)
  }
  if (is.numeric(x = data[, 2]) == FALSE) {
    stop("Entries in the second column must be numerical.",
         call. = FALSE)
  }
  length <- length(x = unique(x = data[, 2]))
  if (length != 2 &&
      length != 1) {
    stop("Second column of data must consist of one or two different values.",
         call. = FALSE)
  }
  unique <- unique(x = data[, 2])
  if (length == 2) {
    if (sort(x = unique)[1] != 0 ||
        sort(x = unique)[2] != 1) {
      stop("Entries in the second column must be 0 and 1.",
           call. = FALSE)
    }
  }
  else {  # length == 1
    if (unique != 1) {
      stop("Entries in the second column must be 1 (and 0).",
           call. = FALSE)
    }
  }
  data0 <- split(x = data,
                 f = data[, 1])
  if (sum(data0[[1]][, 2]) <= 0 ||
      sum(data0[[2]][, 2]) <= 0) {
    stop("There must be minimum one event in each group.",
         call. = FALSE)
  }
  if (is.numeric(x = data[, 3]) == FALSE) {
    stop("Entries in the third column of data must be numerical.",
         call. = FALSE)
  }
  if (min(sort(x = data[, 3])) <= 0) {
    stop("Entries in the third column must be bigger than 0.",
         call. = FALSE)
  }

  # check of cont.time
  if (length(x = cont.time) != 1) {
    stop("cont.time must be a single value.",
         call. = FALSE)
  }
  if (is.numeric(x = cont.time) == FALSE) {
    stop("cont.time must be numerical.",
         call. = FALSE)
  }
  if (cont.time < 0) {
    stop("cont.time must be bigger than or equal to 0.",
         call. = FALSE)
  }

  # check of new.pat
  if (length(x = new.pat) != 2) {
    stop("new.pat must be a two-dimensional vector.",
         call. = FALSE)
  }
  if (is.numeric(x = new.pat) == FALSE) {
    stop("Entries in new.pat must be numerical.",
         call. = FALSE)
  }
  if (new.pat[1] < 0 ||
      new.pat[2] < 0) {
    stop("Entries in new.pat must be bigger than or equal to 0.",
         call. = FALSE)
  }

  # check of theta.0
  if (length(x = theta.0) != 1) {
    stop("theta.0 must be a single value.",
         call. = FALSE)
  }
  if (is.numeric(x = theta.0) == FALSE) {
    stop("theta.0 must be numerical.",
         call. = FALSE)
  }
  if (theta.0 <= 0) {
    stop("theta.0 must be bigger than 0.",
         call. = FALSE)
  }

  # check of alpha
  if (length(x = alpha) != 1) {
    stop("alpha must be a single value.",
         call. = FALSE)
  }
  if (is.numeric(x = alpha) == FALSE) {
    stop("alpha must be numerical.",
         call. = FALSE)
  }
  if (alpha < 0 ||
      alpha > 1) {
    stop("alpha must be choosen between 0 and 1, respectively equal to 0 or 1.",
         call. = FALSE)
  }

  # check of disp.data
  if (length(x = disp.data) != 1) {
    stop("disp.data must be a single value.",
         call. = FALSE)
  }
  if (is.logical(x = disp.data) == FALSE) {
    stop("disp.data must be logical.",
         call. = FALSE)
  }

  # check of plot.km
  if (length(x = plot.km) != 1) {
    stop("plot.km must be a single value.",
         call. = FALSE)
  }
  if (is.logical(x = plot.km) == FALSE) {
    stop("plot.km must be logical.",
         call. = FALSE)
  }
}