suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})

g <- function(rho) bicop_dist("gaussian", parameters = rho)
vine <- vinecop_dist(
  list(
    list(g(.25), g(.45), g(.30)),
    list(g(-.20), g(.15)),
    list(g(.10))
  ),
  dvine_structure(1:4)
)
conditioning <- matrix(c(.4, -.3, -.6, .2), nrow = 2, byrow = TRUE)
margins <- lapply(c(.35, .45, 1, 1), copulaMarginNormal)
margins[[3]]$free[] <- FALSE
margins[[4]]$free[] <- FALSE

implicit_warm_start <- try(copulaSet(
  vine, margins = margins,
  conditioning = list(values = conditioning,
                      variableName = c("WT", "eGFR"))
), silent = TRUE)
stopifnot(inherits(implicit_warm_start, "try-error"))

copulaSet(
  vine, margins = margins,
  conditioning = list(values = conditioning,
                      variableName = c("WT", "eGFR")),
  warmStartOnActivate = FALSE
)

eta <- matrix(c(.1, -.2, -.15, .35), nrow = 2, byrow = TRUE)
conditioning_test <- conditioning[rep(1, 2), , drop = FALSE]
joint <- cbind(eta, conditioning_test)
u_joint <- copulaUeta(joint)

R <- copulaGaussianDvineCor(vine, 4L)
S <- diag(c(.35, .45, 1, 1)) %*% R %*% diag(c(.35, .45, 1, 1))
Spp <- S[1:2, 1:2, drop = FALSE]
Spc <- S[1:2, 3:4, drop = FALSE]
Scc <- S[3:4, 3:4, drop = FALSE]
Scond <- Spp - Spc %*% solve(Scc, t(Spc))
analytic <- vapply(seq_len(nrow(eta)), function(i) {
  mu <- as.vector(Spc %*% solve(Scc, conditioning_test[i, ]))
  z <- eta[i, ] - mu
  .5 * drop(t(z) %*% solve(Scond, z))
}, numeric(1))

stopifnot(max(abs((u_joint - u_joint[1]) -
                  (analytic - analytic[1]))) < 1e-9)
stopifnot(max(abs(copulaOmega() - Spp)) < 1e-9)

copulaPoolUpdate(
  eta, gamma = 1, nchains = 1,
  phiM = eta, xM = matrix(1, 2, 1),
  locMap = matrix(c(1, 1), 1, 2), beta = 0,
  subject = 1:2, iteration = 1L
)
state <- copulaGet()
stopifnot(identical(unname(state$poolConditioning), unname(conditioning)))

cat("joint conditional copula checks passed\n")
