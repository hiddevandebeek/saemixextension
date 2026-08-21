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

## log p(y_i | eta) for the 2-cmt proportional-error model, VECTORISED over
## draws: E is M x d, returns a length-M vector.  One model call for all M
## draws instead of M calls -- the IS loop is otherwise dominated by R overhead.
llObsPKmat <- function(E, y, dose, tim, mu, propErr) {
  E <- as.matrix(E); M <- nrow(E); nt <- length(y)
  psi <- exp(sweep(E, 2, mu, "+"))
  f <- model2cmt(psi, rep(seq_len(M), each = nt),
                 cbind(rep(dose, length.out = M * nt), rep(tim, times = M)))
  f <- matrix(f, nrow = nt, ncol = M)
  bad <- !is.finite(f) | f <= 0
  f[bad] <- 1e-12
  g <- propErr * f
  ll <- colSums(dnorm(rep(y, times = M), as.vector(f), as.vector(g), log = TRUE) |>
                matrix(nrow = nt))
  ll[apply(bad, 2, any)] <- -1e10
  ll
}

llObsPK <- function(eta, y, dose, tim, mu, propErr)
  llObsPKmat(matrix(eta, nrow = 1), y, dose, tim, mu, propErr)

## dPrior(eta) must accept a matrix of etas (rows) and return log densities.
## DEFENSIVE MIXTURE.  A Laplace-t proposal alone gives heavy-tailed weights
## under a copula prior with tail dependence -- measured ESS 77/2000 on the first
## replicate, against 1273 for the Gaussian model.  Mixing a fraction `defensive`
## of draws FROM THE PRIOR bounds the weight ratio by likelihood/defensive, which
## is exactly the failure mode a tail-dependent prior creates.  rPrior must draw
## from the same law dPrior evaluates.
marginalLL <- function(dat, mu, dPrior, propErr, M = 3000, tdf = 5,
                       inflate = 1.3, seed = 1, rPrior = NULL, defensive = 0.25) {
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
    useMix <- !is.null(rPrior) && defensive > 0
    nD <- if (useMix) max(1L, round(M * defensive)) else 0L
    E <- rbind(mvtnorm::rmvt(M - nD, sigma = S, df = tdf, delta = op$par, type = "shifted"),
               if (nD > 0) rPrior(nD) else NULL)
    lqT <- mvtnorm::dmvt(E, delta = op$par, sigma = S, df = tdf, log = TRUE)
    lq <- if (useMix)
      log((1 - defensive) * exp(lqT) + defensive * exp(dPrior(E))) else lqT
    lpri <- dPrior(E)
    lp <- llObsPKmat(E, di$y, di$dose, di$time, mu, propErr) + lpri
    lw <- lp - lq
    mx <- max(lw)
    w  <- exp(lw - mx)
    ll[k]  <- mx + log(mean(w))
    ess[k] <- sum(w)^2 / sum(w^2)
  }
  list(ll = sum(ll), perSubject = ll, ess = ess, essMin = min(ess), essMean = mean(ess))
}

## Convenience priors on the eta scale, each paired with its sampler so the
## defensive mixture can draw from exactly the law it evaluates.
priorMVN  <- function(Omega) function(E) mvtnorm::dmvnorm(E, sigma = Omega, log = TRUE)
rpriorMVN <- function(Omega) function(n) mvtnorm::rmvnorm(n, sigma = Omega)
priorVine <- function(vine, sd) function(E) dEtaVine(E, vine, sd, log = TRUE)
rpriorVine<- function(vine, sd) function(n) rEtaVine(n, vine, sd)
