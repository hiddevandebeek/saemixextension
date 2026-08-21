## Two-compartment IV bolus popPK model, 4 random effects.
##
## d = 4 gives an R-vine with 3 trees / 6 edges -- enough structure for the vine
## to have content, small enough that a fit runs in seconds.  Rich sampling and
## low residual error are deliberate: this is the BEST case for detecting
## non-Gaussian eta dependence.  If it cannot be seen here it cannot be seen.

source("R/etaCopula.R")
source("R/simEta.R")

PK_NAMES <- c("V1", "CL", "Q", "V2")
PK_TRUE  <- c(V1 = 10, CL = 3, Q = 2, V2 = 30)
PK_SD    <- c(0.30, 0.30, 0.40, 0.40)          # eta SD on the log scale
PK_TIMES <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24)
PK_DOSE  <- 100

## Correlation of the etas.  Moderate, PD, and not near-singular.
PK_R <- matrix(c(1,   .55, .30, .10,
                 .55, 1,   .45, .20,
                 .30, .45, 1,   .40,
                 .10, .20, .40, 1), 4, 4)

## Biexponential IV bolus.  psi columns must be in PK_NAMES order.
model2cmt <- function(psi, id, xidep) {
  dose <- xidep[, 1]; tim <- xidep[, 2]
  V1 <- psi[id, 1]; CL <- psi[id, 2]; Q <- psi[id, 3]; V2 <- psi[id, 4]
  k10 <- CL / V1; k12 <- Q / V1; k21 <- Q / V2
  s <- k10 + k12 + k21
  disc <- sqrt(pmax(s * s - 4 * k10 * k21, 0))
  alpha <- (s + disc) / 2; beta <- (s - disc) / 2
  A <- dose / V1 * (alpha - k21) / (alpha - beta)
  B <- dose / V1 * (k21 - beta) / (alpha - beta)
  A * exp(-alpha * tim) + B * exp(-beta * tim)
}

## Simulate one dataset.  `vine` is an rvinecopulib vinecop_dist on 4 variables
## in PK_NAMES order (see the ordering contract in etaCopula.R).
simPK <- function(N, vine, sdEta = PK_SD, propErr = 0.10,
                  times = PK_TIMES, dose = PK_DOSE, truth = PK_TRUE) {
  eta <- rEtaVine(N, vine, sdEta)
  psi <- sweep(exp(eta), 2, truth, "*")
  colnames(psi) <- PK_NAMES
  nt <- length(times)
  d <- data.frame(id = rep(seq_len(N), each = nt),
                  dose = dose, time = rep(times, N))
  f <- model2cmt(psi, d$id, cbind(d$dose, d$time))
  d$y <- f * (1 + propErr * rnorm(nrow(d)))
  d$y[d$y <= 0] <- 1e-6
  list(data = d, eta = eta, psi = psi)
}

## saemix model object matching simPK.  covariance.model is FULL: the Gaussian
## arm must be allowed to fit all 6 correlations, otherwise any "copula wins"
## result is just an unmodelled-correlation artefact.
pkSaemixModel <- function(fullOmega = TRUE) {
  cm <- if (fullOmega) matrix(1, 4, 4) else diag(4)
  saemix::saemixModel(
    model = model2cmt, modeltype = "structural",
    description = "2-cmt IV bolus",
    psi0 = matrix(PK_TRUE, ncol = 4, byrow = TRUE,
                  dimnames = list(NULL, PK_NAMES)),
    transform.par = c(1, 1, 1, 1),
    covariance.model = cm,
    omega.init = diag(PK_SD^2),
    error.model = "proportional", verbose = FALSE)
}

pkSaemixData <- function(d)
  saemix::saemixData(name.data = d, header = TRUE,
                     name.group = "id", name.predictors = c("dose", "time"),
                     name.response = "y", verbose = FALSE)
