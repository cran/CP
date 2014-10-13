PersMonNonMixExp <- function(lambda, c, n.alive, new.pat, cont.time) {

  ## Calculates the further person months
  ## in the non-mixture model with exponential survival, i. e.
  ##   S(t) = c^[1 - exp(- lambda * t)], lambda > 0, 0 < c < 1, t >= 0,
  ## from the interim analysis until the end of time of continuation.
  ##
  ## Args:
  ##   lambda: Parameter.
  ##   c: Parameter.
  ##   n.alive: Number of patients still alive.
  ##   new.pat: Number of patients who will be recruited each time unit.
  ##   cont.time: Period of time of continuing the trial.
  ##
  ## Results:
  ##   Further person months.

  # auxiliary functions / integrands
  Integrand1 <- function(x) {
    c^(1 - exp(- lambda * x))
  }
  Integrand2 <- function(x) {
    c^exp(- lambda * x)
  }
  Integrand3 <- function(x) {
    (integrate(f     = Integrand2,
               lower = 0,
               upper = x)$value * c^(1 - exp(- lambda * x)))
  }

  # calculation of further person months
  O.star <- (n.alive * integrate(f     = Integrand1,
                                 lower = 0,
                                 upper = cont.time)$value
              + new.pat / c * integrate(f     = Vectorize(Integrand3),
                                        lower = 0,
                                        upper = cont.time)$value)

  return(O.star)
}