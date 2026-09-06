## Arbitrary-margin generalized-EM blocks must retain the literal common Q;
## the final simultaneous polish remains available for stationarity.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})

n_fail <- 0L
ok <- function(label, pass, detail = "") {
  if (!isTRUE(pass)) n_fail <<- n_fail + 1L
  cat(sprintf("%-70s %s %s\n", label,
              if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

g <- function(rho) bicop_dist("gaussian", parameters = rho)
truth <- vinecop_dist(list(
  list(bicop_dist("frank", parameters = 3),
       bicop_dist("clayton", rotation = 180, parameters = 1.8), g(.55)),
  list(bicop_dist("joe", rotation = 270, parameters = 1.35),
       bicop_dist("frank", parameters = -2.2)),
  list(g(.18))), dvine_structure(1:4))
margins <- list(copulaMarginNormal(.30), copulaMarginNormal(.42),
  copulaMarginCovariateLognormal(log(72), .20),
  copulaMarginCovariateWeibull(3.2, 98))
set.seed(82451)
E <- copulaMarginsQuantile(rvinecop(1800L, truth), margins)
w <- rep(1 / nrow(E), nrow(E))

start_margins <- margins
start_margins[[1]] <- copulaMarginNormal(.38)
start_margins[[2]] <- copulaMarginNormal(.34)
start_margins[[3]] <- copulaMarginCovariateLognormal(log(68), .28)
start_margins[[4]] <- copulaMarginCovariateWeibull(2.5, 90)
X <- matrix(1, nrow(E), 2)
loc_map <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
beta0 <- c(.08, -.06)
residual <- E - copulaLocation(X, beta0, loc_map)
q0 <- sum(w * copulaLogPrior(residual, truth, margins = start_margins,
                              numericalPolicy = "exact"))

block <- copulaMaximiseJoint(E, w, copulaMarginScales(start_margins), truth,
  4L, maxit = 12L, X = X, locMap = loc_map, beta0 = beta0,
  betaFree = 1:2, cores = 1L, engine = "hybrid", cycles = 1L,
  polishMaxit = 0L, margins = start_margins)
block_residual <- E - copulaLocation(X, block$beta, loc_map)
q_block <- sum(w * copulaLogPrior(block_residual, block$vine,
  margins = block$margins, numericalPolicy = "exact"))
ok("HA1 arbitrary-margin hybrid uses the block engine",
   identical(block$backend, "hybrid-margins"), block$backend)
ok("HA2 every retained block is common-Q nondecreasing", q_block >= q0 - 1e-10,
   sprintf("%.8f -> %.8f", q0, q_block))

polished <- copulaMaximiseJoint(E, w, copulaMarginScales(start_margins), truth,
  4L, maxit = 12L, X = X, locMap = loc_map, beta0 = beta0,
  betaFree = 1:2, cores = 1L, engine = "hybrid", cycles = 1L,
  polishMaxit = 12L, margins = start_margins, diagnostics = TRUE)
polished_residual <- E - copulaLocation(X, polished$beta, loc_map)
q_polished <- sum(w * copulaLogPrior(polished_residual, polished$vine,
  margins = polished$margins, numericalPolicy = "exact"))
ok("HA3 simultaneous polish retains or improves the block result",
   q_polished >= q_block - 2e-7,
   sprintf("block=%.8f polished=%.8f", q_block, q_polished))

cat(sprintf("\n%d failure(s)\n", n_fail))
if (n_fail > 0L) quit(status = 1L)

