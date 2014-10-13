PersMonExp <- function(d, o, n.alive, new.pat, cont.time) {

  ## Calculates the further person months
  ## in the exponential survival, i. e.
  ##   S(t) = exp( - lambda * t), lambda > 0, t >= 0,
  ## from the interim analysis until the end of time of continuation.
  ##
  ## Args:
  ##   d: Parameter.
  ##   o: Parameter.
  ##   n.alive: Number of patients still alive.
  ##   new.pat: Number of patients who will be recruited each time unit.
  ##   cont.time: Period of time of continuing the trial.
  ##
  ## Results:
  ##   Further person months.

  # calculation of further person months
  O.star <- (n.alive * (1 - exp(- d / o * cont.time)) / (d / o)
              + new.pat * cont.time / (d / o)
              - new.pat * (1 - exp(- d / o * cont.time)) / (d / o)^2)

  return(O.star)
}