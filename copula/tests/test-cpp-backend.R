## Persistent C++ and reference R objectives must define the same fixed-vine Q.

suppressMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-60s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

set.seed(821)
d <- 4L; n <- 1200L
R <- matrix(.35, d, d); diag(R) <- 1
vine <- makeTwinVines(R, "gumbel")$alt
sd0 <- rep(.34, d)
X <- cbind(1, rep(c(-1, 1), length.out = n), 1, 1, 1)
locMap <- matrix(0, nrow = ncol(X), ncol = d)
locMap[1, 1] <- 1; locMap[2, 1] <- 1
locMap[3, 2] <- 1; locMap[4, 3] <- 1; locMap[5, 4] <- 1
beta <- c(.2, .15, -.1, .05, .3)
phi <- copulaLocation(X, beta, locMap) + rEtaVine(n, vine, rep(.3, d))
w <- rep(1 / n, n)
flat <- copulaPadFlat(vine, d); layout <- copulaParameterLayout(flat)
locFree <- seq_along(beta)
p <- c(beta, log(sd0), layout$par)
ctx <- copula_q_context_create(vine, phi, X, locMap, beta, locFree, w)
qCpp <- -copula_q_context_eval(ctx, p, 1L)
qR <- sum(w * copulaLogPrior(phi - copulaLocation(X, beta, locMap), vine, sd0))
ok("C++ fused Q equals the R reference with covariates",
   abs(qCpp - qR) < 1e-10, sprintf("difference=%.3g", qCpp - qR))

r <- copulaMaximiseJoint(phi, w, sd0, vine, d, maxit = 8L,
  X = X, locMap = locMap, beta0 = beta, backend = "r")
c <- copulaMaximiseJoint(phi, w, sd0, vine, d, maxit = 8L,
  X = X, locMap = locMap, beta0 = beta, backend = "cpp")
cLiteral <- sum(w * copulaLogPrior(
  phi - copulaLocation(X, c$beta, locMap), c$vine, c$sd))
rLiteral <- sum(w * copulaLogPrior(
  phi - copulaLocation(X, r$beta, locMap), r$vine, r$sd))
ok("C++ and R finite fits report the same literal joint-Q target",
   abs(c$value - cLiteral) < 1e-10 &&
     abs(r$value - rLiteral) < 1e-10 &&
     c$value >= qR - 1e-10 && r$value >= qR - 1e-10,
   sprintf("R=%.10f C++=%.10f", r$value, c$value))
ok("C++ result preserves the vine class and discrete model",
   inherits(c$vine, "vinecop_dist") &&
     identical(copulaFingerprint(c$vine, d), copulaFingerprint(vine, d)))

## The common PK case has two random effects.  Auto-dispatch must use the same
## compiled objective here too (the old d >= 3 guard caused the 1000-iteration
## benchmark to evaluate the literal R objective at every optimizer step).
set.seed(822)
d2 <- 2L; n2 <- 1800L
v2 <- etaVine(list(rvinecopulib::bicop_dist(
  "clayton", parameters = 1.35)), d2)
sd2 <- c(.28, .42)
X2 <- cbind(1, rep(c(-1, 1), length.out = n2), 1)
L2 <- matrix(0, nrow = ncol(X2), ncol = d2)
L2[1, 1] <- 1; L2[2, 1] <- 1; L2[3, 2] <- 1
b2 <- c(.15, -.08, .32)
phi2 <- copulaLocation(X2, b2, L2) + rEtaVine(n2, v2, sd2)
w2 <- rep(1 / n2, n2)
flat2 <- copulaPadFlat(v2, d2); layout2 <- copulaParameterLayout(flat2)
p2 <- c(b2, log(sd2), layout2$par)
ctx2 <- copula_q_context_create(v2, phi2, X2, L2, b2, seq_along(b2), w2)
qCpp2 <- -copula_q_context_eval(ctx2, p2, 1L)
qR2 <- sum(w2 * copulaLogPrior(
  phi2 - copulaLocation(X2, b2, L2), v2, sd2))
ok("2D C++ fused Q equals the R reference",
   abs(qCpp2 - qR2) < 1e-10, sprintf("difference=%.3g", qCpp2 - qR2))

r2 <- copulaMaximiseJoint(phi2, w2, sd2, v2, d2, maxit = 12L,
  X = X2, locMap = L2, beta0 = b2, backend = "r")
c2 <- copulaMaximiseJoint(phi2, w2, sd2, v2, d2, maxit = 12L,
  X = X2, locMap = L2, beta0 = b2, backend = "cpp")
a2 <- copulaMaximiseJoint(phi2, w2, sd2, v2, d2, maxit = 12L,
  X = X2, locMap = L2, beta0 = b2, backend = "auto")
cLiteral2 <- sum(w2 * copulaLogPrior(
  phi2 - copulaLocation(X2, c2$beta, L2), c2$vine, c2$sd))
rLiteral2 <- sum(w2 * copulaLogPrior(
  phi2 - copulaLocation(X2, r2$beta, L2), r2$vine, r2$sd))
ok("2D C++ and R finite fits report the same literal joint-Q target",
   abs(c2$value - cLiteral2) < 1e-10 &&
     abs(r2$value - rLiteral2) < 1e-10 &&
     c2$value >= qR2 - 1e-10 && r2$value >= qR2 - 1e-10,
   sprintf("R=%.10f C++=%.10f", r2$value, c2$value))
ok("2D non-Gaussian auto-dispatch selects C++",
   identical(a2$backend, "cpp"), paste("backend=", a2$backend))

q1 <- -copula_q_context_eval(ctx2, p2, 1L)
q4 <- -copula_q_context_eval(ctx2, p2, copulaResolveCores("auto"))
ok("Parallel and serial C++ Q evaluations agree",
   abs(q1 - q4) < 1e-10, sprintf("difference=%.3g", q4 - q1))

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
