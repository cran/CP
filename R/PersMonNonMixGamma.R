PersMonNonMixGamma <- function(a, b, c, n.alive, new.pat, cont.time) {

  ## Calculates the further person months
  ## in the non-mixture model with Gamma type survival, i. e.
  ##   S(t) = c^Gamma^(0)(a, b * t), a > 0, b > 0, 0 < c < 1, t >= 0,
  ## with Gamma^(0) being the regularized incomplete Gamma function of the upper bound,
  ## from the interim analysis until the end of time of continuation.
  ##
  ## Args:
  ##   a: Parameter.
  ##   b: Parameter.
  ##   c: Parameter.
  ##   n.alive: Number of patients still alive.
  ##   new.pat: Number of patients who will be recruited each time unit.
  ##   cont.time: Period of time of continuing the trial.
  ##
  ## Results:
  ##   Further person months.

  # auxiliary functions / integrands
  Integrand1 <- function(x) {
    c^stats::pgamma(q     = x,
                    shape = a,
                    rate  = b)
  }
  Integrand2 <- function(x) {
    c^(1 - stats::pgamma(q     = x,
                         shape = a,
                         rate  = b))
  }
  Integrand3 <- function(x) {
    (stats::integrate(f     = Integrand2,
                      lower = 0,
                      upper = x)$value * c^stats::pgamma(q     = x,
                                                         shape = a,
                                                         rate  = b))
  }

  # calculation of further person months
  O.star <- (n.alive * stats::integrate(f     = Integrand1,
                                        lower = 0,
                                        upper = cont.time)$value
              + new.pat / c * stats::integrate(f     = Vectorize(Integrand3),
                                               lower = 0,
                                               upper = cont.time)$value)

  return(O.star)
}