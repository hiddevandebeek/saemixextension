suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})

g <- function(rho) bicop_dist("gaussian", parameters = rho)
vine <- vinecop_dist(
  list(
    list(g(.25), g(-.15), g(.40)),
    list(g(.20), g(-.30)),
    list(g(.10))
  ),
  cvine_structure(1:4)
)
margins <- rep(list(copulaMarginNormal(1)), 4)

stopifnot(copulaConditioningIsTail(vine, 3:4))
subvine <- copulaConditioningSubvine(vine, 3:4)
stopifnot(identical(as.integer(subvine$structure$order), 1:2))

## Midpoint integration on the copula scale verifies both the extracted
## covariate marginal and normalization of p(eta | c).
u <- seq(.0025, .9975, length.out = 200)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
uc <- c(.37, .71)
joint_u <- cbind(grid, matrix(uc, nrow(grid), 2L, byrow = TRUE))
full <- dvinecop(joint_u, vine)
covariate_density <- dvinecop(matrix(uc, nrow = 1L), subvine)
stopifnot(abs(mean(full) / covariate_density - 1) < 1e-3)

E <- qnorm(joint_u)
conditional_lp <- copulaConditionalLogPrior(E, vine, margins, dEta = 2L)
eta_jacobian <- rowSums(dnorm(E[, 1:2, drop = FALSE], log = TRUE))
stopifnot(abs(mean(exp(conditional_lp - eta_jacobian)) - 1) < 1e-3)

set.seed(91)
sample_u <- rvinecop(300, vine)
sample_E <- qnorm(sample_u)
start_vine <- vinecop_dist(
  lapply(vine$pair_copulas, function(tree)
    lapply(tree, function(edge) {
      if (edge$family == "indep") edge else g(.05)
    })), vine$structure)
base_q <- sum(copulaConditionalLogPrior(
  sample_E, start_vine, margins, dEta = 2L))
updated <- copulaMaximiseJoint(
  sample_E, rep(1, nrow(sample_E)), rep(1, 4), start_vine, 4,
  maxit = 5L, withMu = FALSE, margins = margins,
  optimizeMargins = FALSE, optimizeVine = TRUE,
  likelihoodTarget = "conditional", dEta = 2L)
stopifnot(is.finite(updated$value), updated$value >= base_q - 1e-8)

api_margins <- list(copulaMarginNormal(.3), copulaMarginNormal(.4),
                    copulaMarginCovariateNormal(70, 12),
                    copulaMarginCovariateNormal(95, 20))
conditioning <- cbind(c(65, 75), c(80, 110))
copulaSet(vine, margins = api_margins, conditioning = conditioning,
          likelihoodTarget = "conditional", warmStartOnActivate = FALSE)
stopifnot(identical(copulaGet()$likelihoodTarget, "conditional"))

bad <- vinecop_dist(vine$pair_copulas, dvine_structure(c(1, 3, 2, 4)))
stopifnot(!copulaConditioningIsTail(bad, 3:4))
failed <- try(copulaConditionalLogPrior(E, bad, margins, dEta = 2L),
              silent = TRUE)
stopifnot(inherits(failed, "try-error"))
gaussian_api <- try(copulaSet(
  bad, margins = api_margins, conditioning = conditioning,
  likelihoodTarget = "conditional", warmStartOnActivate = FALSE), silent = TRUE)
stopifnot(!inherits(gaussian_api, "try-error"))
stopifnot(all(is.finite(copulaGaussianFremLogPrior(
  cbind(matrix(c(.1, -.2, .2, .1), nrow = 2, byrow = TRUE), conditioning),
  bad, api_margins, 2L, likelihoodTarget = "conditional"))))

cat("conditional likelihood target checks passed\n")
