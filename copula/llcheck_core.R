## Shared machinery for the "does it still target the true likelihood" check.
## Everything here is a NEW file; nothing in copula/R or the package is touched.
suppressMessages({library(mvtnorm); library(numDeriv); library(rvinecopulib)})

## Generic vectorised observation log-likelihood: any structural model f(psi,id,xidep)
## with proportional error.  E is M x d.
llObsGen <- function(E, y, dose, tim, mu, propErr, mfun) {
  E <- as.matrix(E); M <- nrow(E); nt <- length(y)
  psi <- exp(sweep(E, 2, mu, "+"))
  f <- mfun(psi, rep(seq_len(M), each = nt),
            cbind(rep(dose, length.out = M * nt), rep(tim, times = M)))
  f <- matrix(f, nrow = nt, ncol = M)
  bad <- !is.finite(f) | f <= 0
  f[bad] <- 1e-12
  g <- propErr * f
  ll <- colSums(matrix(dnorm(rep(y, times = M), as.vector(f), as.vector(g), log = TRUE), nrow = nt))
  ll[apply(bad, 2, any)] <- -1e10
  ll
}

## One importance-sampling pass per subject.  Returns BOTH the marginal
## log-likelihood AND the self-normalised posterior moments / weighted PIT
## sample, so a single pass serves the likelihood oracle and the IFM iteration.
## Common random numbers (fixed seed) make everything smooth in theta.
isPass <- function(dat, mu, sdv, vine, propErr, mfun, M = 2000, tdf = 5,
                   inflate = 1.3, seed = 1, defensive = 0.25, keepDraws = TRUE) {
  set.seed(seed)
  ids <- unique(dat$id); d <- length(mu)
  dPrior <- function(E) dEtaVine(E, vine, sdv, log = TRUE)
  rPrior <- function(n) rEtaVine(n, vine, sdv)
  ll <- ess <- numeric(length(ids))
  m1 <- matrix(0, length(ids), d); m2 <- matrix(0, length(ids), d)
  Eall <- if (keepDraws) vector("list", length(ids)) else NULL
  Wall <- if (keepDraws) vector("list", length(ids)) else NULL
  for (k in seq_along(ids)) {
    di <- dat[dat$id == ids[k], ]
    negp <- function(e) -(llObsGen(matrix(e, nrow = 1), di$y, di$dose, di$time, mu, propErr, mfun) +
                          dPrior(matrix(e, nrow = 1)))
    op <- optim(rep(0, d), negp, method = "BFGS", control = list(maxit = 500))
    H <- tryCatch(numDeriv::hessian(negp, op$par), error = function(e) NULL)
    S <- NULL
    if (!is.null(H) && all(is.finite(H))) {
      S <- tryCatch(solve((H + t(H)) / 2), error = function(e) NULL)
      if (!is.null(S) && (any(!is.finite(S)) || min(eigen(S, TRUE, TRUE)$values) <= 0)) S <- NULL
    }
    if (is.null(S)) S <- diag(0.1, d)
    S <- S * inflate
    nD <- max(1L, round(M * defensive))
    E <- rbind(mvtnorm::rmvt(M - nD, sigma = S, df = tdf, delta = op$par, type = "shifted"),
               rPrior(nD))
    lqT <- mvtnorm::dmvt(E, delta = op$par, sigma = S, df = tdf, log = TRUE)
    lq <- log((1 - defensive) * exp(lqT) + defensive * exp(dPrior(E)))
    lp <- llObsGen(E, di$y, di$dose, di$time, mu, propErr, mfun) + dPrior(E)
    lw <- lp - lq; mx <- max(lw); w <- exp(lw - mx)
    ll[k] <- mx + log(mean(w)); ess[k] <- sum(w)^2 / sum(w^2)
    wn <- w / sum(w)
    m1[k, ] <- colSums(wn * E); m2[k, ] <- colSums(wn * E^2)
    if (keepDraws) { Eall[[k]] <- E; Wall[[k]] <- wn }
  }
  list(ll = sum(ll), perSubject = ll, ess = ess, essMin = min(ess),
       m1 = m1, m2 = m2, E = Eall, W = Wall)
}

## Rebuild a d=2 vine from a single Kendall tau, family/rotation held fixed.
vine2FromTau <- function(tau, family, rotation) {
  tau <- sign(tau) * min(abs(tau), 0.95)
  if (family == "indep") return(etaVine(list(rvinecopulib::bicop_dist("indep")), 2))
  tmpl <- rvinecopulib::bicop_dist(family, rotation = rotation,
            parameters = rvinecopulib::ktau_to_par(family, if (rotation %in% c(90,270)) -0.5 else 0.5))
  etaVine(list(rvinecopulib::bicop_dist(family, rotation = rotation,
            parameters = rvinecopulib::ktau_to_par(tmpl, tau))), 2)
}

## Unconstrained parameter vector p = (mu1, mu2, log sd1, log sd2, atanh tau)
unpack <- function(p) list(mu = p[1:2], sdv = exp(p[3:4]), tau = tanh(p[5]))
pack   <- function(mu, sdv, tau) c(mu, log(sdv), atanh(tau))
