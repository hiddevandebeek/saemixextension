source("copula/tests/helper-load.R")

u <- seq(.001, .999, length.out = 999L)
typical <- matrix(3.5, nrow = length(u), ncol = 1L)
for (margin in list(copulaNaturalMarginLognormal(.45),
                    copulaNaturalMarginGamma(2.5),
                    copulaNaturalMarginWeibull(2.2))) {
  psi <- copulaNaturalMarginsQuantile(matrix(u, ncol = 1L), typical,
    list(margin))
  direct <- margin$log_density(psi[, 1L], typical[, 1L], margin$parameters)

  phiIdentity <- copulaNaturalToWorking(psi, 0L)
  recoveredIdentity <- copulaWorkingToNatural(phiIdentity, 0L)
  logIdentity <- direct + copulaWorkingLogJacobian(phiIdentity, 0L)[, 1L]

  phiLog <- copulaNaturalToWorking(psi, 1L)
  recoveredLog <- copulaWorkingToNatural(phiLog, 1L)
  logMapped <- direct + copulaWorkingLogJacobian(phiLog, 1L)[, 1L]

  stopifnot(max(abs(recoveredIdentity - psi)) < 2e-12,
    max(abs(recoveredLog - psi)) < 2e-12,
    max(abs(logIdentity - direct)) < 2e-12,
    max(abs(logMapped - (direct + log(psi[, 1L])))) < 2e-12,
    max(abs(margin$cdf(psi[, 1L], typical[, 1L], margin$parameters) - u)) <
      2e-10)
}

## The same natural Gamma law integrates to one in identity and log working
## coordinates; only the Jacobian differs.
gamma <- copulaNaturalMarginGamma(2.5)
identityIntegral <- integrate(function(phi)
  exp(gamma$log_density(phi, rep(3.5, length(phi)), gamma$parameters)),
  0, Inf, rel.tol = 1e-8)$value
logIntegral <- integrate(function(phi) {
  psi <- exp(phi)
  exp(gamma$log_density(psi, rep(3.5, length(phi)), gamma$parameters) + phi)
}, -Inf, Inf, rel.tol = 1e-8)$value
stopifnot(abs(identityIntegral - 1) < 2e-8,
  abs(logIntegral - 1) < 2e-8)

## E-step densities under identity and log coordinates are the same natural
## Gamma law after their respective Jacobians.
vine1 <- copulaGaussianRvineFromCor(diag(2),
  rvinecopulib::cvine_structure(1:2))
typicalNatural <- cbind(rep(3.5, length(u)), rep(2.2, length(u)))
psiGamma <- cbind(
  gamma$quantile(u, typicalNatural[, 1L], gamma$parameters),
  gamma$quantile(rev(u), typicalNatural[, 2L], gamma$parameters))
etaIdentity <- psiGamma - typicalNatural
etaLog <- log(psiGamma) - log(typicalNatural)
identityKernel <- copulaNaturalWorkingPriorKernel(vine1, list(gamma, gamma),
  typicalNatural, c(0L, 0L))
logKernel <- copulaNaturalWorkingPriorKernel(vine1, list(gamma, gamma),
  log(typicalNatural), c(1L, 1L))
directGamma <- -rowSums(cbind(
  gamma$log_density(psiGamma[, 1L], typicalNatural[, 1L], gamma$parameters),
  gamma$log_density(psiGamma[, 2L], typicalNatural[, 2L], gamma$parameters)))
stopifnot(max(abs(identityKernel$negative(etaIdentity) - directGamma)) < 2e-11,
  max(abs(logKernel$negative(etaLog) -
    (directGamma - rowSums(log(psiGamma))))) < 2e-11)

cat("natural parameter mapping checks passed for identity and log coordinates\n")
