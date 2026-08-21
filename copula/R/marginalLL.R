## Independent marginal log-likelihood under an ARBITRARY eta prior.
##
##   l_i(theta) = log INT p(y_i | psi = exp(mu + eta)) p(eta; theta) d eta
##
## Deliberately standalone (no saemix internals) so it can serve as an oracle
## for both arms.  Importance sampling with a multivariate-t proposal centred at
## the conditional mode with the Laplace covariance; the t tails are there
## because a copula prior with tail dependence is exactly the case where a
## Gaussian proposal gives heavy-tailed weights.  ESS is returned so that
## failure is visible rather than silent.

suppressMessages({library(mvtnorm); library(numDeriv)})
source("R/etaCopula.R")

## log p(y_i | eta) for the 2-cmt proportional-error model.
llObsPK <- function(eta, y, dose, tim, mu, propErr) {
  psi <- exp(mu + eta)
  f <- model2cmt(matrix(psi, nrow = 1), rep(1, length(y)), cbind(dose, tim))
  if (any(!is.finite(f)) || any(f <= 0)) return(-1e10)
  g <- propErr * f
  sum(dnorm(y, f, g, log = TRUE))
}

## dPrior(eta) must accept a matrix of etas (rows) and return log densities.
marginalLL <- function(dat, mu, dPrior, propErr, M = 3000, tdf = 5,
                       inflate = 1.3, seed = 1) {
  set.seed(seed)
  ids <- unique(dat$id); d <- length(mu)
  ll <- ess <- numeric(length(ids))
  for (k in seq_along(ids)) {
    di <- dat[dat$id == ids[k], ]
    negp <- function(e) -(llObsPK(e, di$y, di$dose, di$time, mu, propErr) +
                          dPrior(matrix(e, nrow = 1)))
    op <- optim(rep(0, d), negp, method = "BFGS", control = list(maxit = 500))
    H  <- tryCatch(numDeriv::hessian(negp, op$par), error = function(e) NULL)
    S  <- NULL
    if (!is.null(H) && all(is.finite(H))) {
      S <- tryCatch(solve((H + t(H)) / 2), error = function(e) NULL)
      if (!is.null(S) && (any(!is.finite(S)) || min(eigen(S, TRUE, TRUE)$values) <= 0)) S <- NULL
    }
    if (is.null(S)) S <- diag(0.1, d)
    S <- S * inflate
    E <- mvtnorm::rmvt(M, sigma = S, df = tdf, delta = op$par, type = "shifted")
    lq <- mvtnorm::dmvt(E, delta = op$par, sigma = S, df = tdf, log = TRUE)
    lp <- apply(E, 1, function(e) llObsPK(e, di$y, di$dose, di$time, mu, propErr)) +
          dPrior(E)
    lw <- lp - lq
    mx <- max(lw)
    w  <- exp(lw - mx)
    ll[k]  <- mx + log(mean(w))
    ess[k] <- sum(w)^2 / sum(w^2)
  }
  list(ll = sum(ll), perSubject = ll, ess = ess, essMin = min(ess), essMean = mean(ess))
}

## Convenience priors on the eta scale.
priorMVN  <- function(Omega) function(E) mvtnorm::dmvnorm(E, sigma = Omega, log = TRUE)
priorVine <- function(vine, sd) function(E) dEtaVine(E, vine, sd, log = TRUE)
