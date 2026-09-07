## The prior kernels hoist everything that does not depend on eta and remember
## the margin evaluation one column at a time, because the second MCMC kernel
## moves a single coordinate per proposal.
##
## Caching is only correct if it is invisible, so the memo itself is checked
## exactly: a batch the memo has already seen, after another batch moved
## through it, must come back bit-identical. Against the general log-prior the
## joint kernel is exact too, but the conditional one is not -- it splits the
## quadratic form into the block that moves with eta, the block that is linear
## in it and the block that does not involve it, which is a different
## association order for the same sum. That is checked to 1e-12 rather than to
## zero, and the reason is recorded here so the tolerance is not later mistaken
## for a defect being tolerated.

suppressPackageStartupMessages({
  library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})

nFail <- 0L
## Rejected rows come back as Inf from both sides, and Inf - Inf is NaN, so
## equal values are compared as equal before any subtraction happens.
gap <- function(a, b) {
  difference <- abs(a - b)
  difference[a == b] <- 0
  max(difference)
}
ok <- function(label, difference, tolerance = 0) {
  pass <- isTRUE(difference <= tolerance)
  if (!pass) nFail <<- nFail + 1L
  cat(sprintf("%-46s %s %.3e\n", label, if (pass) "PASS" else "**FAIL**",
    difference))
}

set.seed(9)
n <- 60L; dEta <- 2L
predictor <- cbind(rep(log(2.1), n), rep(log(3.5), n))
conditioning <- matrix(rgamma(n, 5, 1), n, 1L)
transform <- c(1L, 1L)
pair <- function(rho) rvinecopulib::bicop_dist("gaussian", 0, rho)
vine <- rvinecopulib::vinecop_dist(
  structure = rvinecopulib::dvine_structure(1:3),
  pair_copulas = list(list(pair(.4), pair(.3)), list(pair(.1))))
margins <- list(copulaNaturalMarginLognormal(.3),
  copulaNaturalMarginGamma(4),
  copulaFitCovariateMargin(as.numeric(conditioning), "gamma"))
eta <- matrix(rnorm(n * dEta, 0, .2), n, dEta)
moved <- eta; moved[, 2L] <- moved[, 2L] + .05

kernel <- copulaNaturalFremConditionalKernel(conditioning, vine, margins,
  dEta, predictor, transform)
reference <- function(value) -copulaNaturalFremLogPrior(value, conditioning,
  vine, margins, dEta, predictor, transform, "conditional")
ok("conditional kernel matches the reference",
  gap(kernel$negative(eta), reference(eta)), 1e-12)
ok("after moving one column",
  gap(kernel$negative(moved), reference(moved)), 1e-12)

## the joint kernel, with no conditioning
jointVine <- rvinecopulib::vinecop_dist(
  structure = rvinecopulib::dvine_structure(1:2),
  pair_copulas = list(list(pair(.4))))
joint <- copulaNaturalWorkingPriorKernel(jointVine, margins[1:2], predictor,
  transform)
jointReference <- function(value) -copulaNaturalFremLogPrior(value,
  matrix(numeric(), n, 0L), jointVine, margins[1:2], dEta, predictor,
  transform, "joint")
ok("joint working kernel", gap(joint$negative(eta), jointReference(eta)))
ok("joint kernel after one column moves",
  gap(joint$negative(moved), jointReference(moved)))

## The memo has to be exactly invisible whatever the arithmetic around it does.
first <- kernel$negative(eta)
invisible(kernel$negative(moved))
ok("the memo is bit-identical on return",
  if (identical(kernel$negative(eta), first)) 0 else 1)
jointFirst <- joint$negative(eta)
invisible(joint$negative(moved))
ok("the joint memo is bit-identical on return",
  if (identical(joint$negative(eta), jointFirst)) 0 else 1)

## Both kernels take a fast path when every row is valid, which is the ordinary
## case and therefore the only one the checks above reach. An eta large enough
## to send the natural parameter to infinity exercises the other branch.
extreme <- eta; extreme[1L, 1L] <- 800
ok("conditional kernel with an invalid row",
  gap(kernel$negative(extreme), reference(extreme)), 1e-12)
ok("joint kernel with an invalid row",
  gap(joint$negative(extreme), jointReference(extreme)))
if (!is.infinite(kernel$negative(extreme)[1L])) {
  nFail <- nFail + 1L
  cat("**FAIL** the invalid row was not rejected\n")
}

cat(if (nFail) sprintf("\n%d FAILURES\n", nFail) else
  "\nall prior-kernel cache checks passed\n")
if (nFail) quit(status = 1L)
