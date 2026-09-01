source("copula/tests/helper-load.R")

R <- matrix(c(1, .20, .45, .10,
              .20, 1, .25, .55,
              .45, .25, 1, .15,
              .10, .55, .15, 1), 4L, 4L, byrow = TRUE)
structure <- rvinecopulib::cvine_structure(c(3, 1, 4, 2))
vine <- copulaGaussianRvineFromCor(R, structure)
conditioning <- cbind(WT = seq(55, 85, length.out = 30L),
  eGFR = seq(60, 120, length.out = 30L))
margins <- list(copulaMarginNormal(.2), copulaMarginNormal(.3),
  copulaMarginCovariateNormal(70, 10),
  copulaMarginCovariateNormal(90, 20))
copulaSet(vine, margins = margins,
  conditioning = list(values = conditioning,
    variableName = c("WT", "eGFR")),
  warmStartOnActivate = FALSE, guard = FALSE)
state <- copulaSnapshot(copulaGet(), variableName = c("eta_V", "eta_CL"))

ffem <- copulaFremToFfem(state)
D <- diag(c(.2, .3, 10, 20)); S <- D %*% R %*% D
B <- S[1:2, 3:4, drop = FALSE] %*% solve(S[3:4, 3:4, drop = FALSE])
Omega <- S[1:2, 1:2, drop = FALSE] - B %*% S[3:4, 1:2, drop = FALSE]
stopifnot(identical(ffem$representation, "linear"),
  isTRUE(ffem$linearAvailable),
  max(abs(ffem$native$coefficient - B)) < 1e-12,
  max(abs(ffem$omega - Omega)) < 1e-12)
z <- copulaFfemTransform(ffem, data.frame(WT = 80, eGFR = 110))
stopifnot(max(abs(z - 1)) < 1e-12,
  max(abs(copulaFfemTransform(ffem, c(WT = 80, eGFR = 110)) - 1)) < 1e-12,
  max(abs(copulaFfemLocation(ffem, data.frame(WT = 80, eGFR = 110)) -
    z %*% t(ffem$coefficient))) < 1e-12)

subsetFfem <- copulaFremToFfem(state, "WT")
Ssub <- S[c(1, 2, 3), c(1, 2, 3), drop = FALSE]
Bsub <- Ssub[1:2, 3, drop = FALSE] / Ssub[3, 3]
Osub <- Ssub[1:2, 1:2, drop = FALSE] -
  Bsub %*% Ssub[3, 1:2, drop = FALSE]
stopifnot(max(abs(subsetFfem$native$coefficient - Bsub)) < 1e-12,
  max(abs(subsetFfem$omega - Osub)) < 1e-12)

gammaMargins <- margins
gammaMargins[[3L]] <- copulaMarginCovariateGamma(4, 14)
copulaSet(vine, margins = gammaMargins,
  conditioning = list(values = conditioning,
    variableName = c("WT", "eGFR")),
  warmStartOnActivate = FALSE, guard = FALSE)
gammaState <- copulaSnapshot(copulaGet(), variableName = c("eta_V", "eta_CL"))
gammaFfem <- copulaFremToFfem(gammaState, "WT")
raw <- matrix(qgamma(c(.2, .8), 4, scale = 14), ncol = 1L,
  dimnames = list(NULL, "WT"))
stopifnot(isTRUE(gammaFfem$linearAvailable), is.null(gammaFfem$native),
  max(abs(copulaFfemTransform(gammaFfem, raw) -
    qnorm(c(.2, .8)))) < 1e-10)

generalMargins <- gammaMargins
generalMargins[[1L]] <- copulaMarginCenteredGamma(4, .2)
copulaSet(vine, margins = generalMargins,
  conditioning = list(values = conditioning,
    variableName = c("WT", "eGFR")),
  warmStartOnActivate = FALSE, guard = FALSE)
generalState <- copulaSnapshot(copulaGet(), variableName = c("eta_V", "eta_CL"))
generalFfem <- copulaFremToFfem(generalState, "WT")
stopifnot(identical(generalFfem$representation, "distributional"),
  !isTRUE(generalFfem$linearAvailable))
set.seed(481L)
translated <- copulaFfemSimulate(generalFfem, raw, n = 3L)
full <- cbind(WT = raw[rep(seq_len(nrow(raw)), each = 3L), 1L],
  eGFR = NA_real_)
set.seed(481L)
direct <- copulaGaussianFremRandEta(full, generalState$vine,
  generalState$margins, generalState$dEta)
stopifnot(max(abs(translated - direct)) == 0)

vine2 <- copulaGaussianRvineFromCor(matrix(c(1, .5, .5, 1), 2L),
  rvinecopulib::cvine_structure(c(2, 1)))
categorical <- matrix(rep(c(0, 1), 15L), ncol = 1L,
  dimnames = list(NULL, "SEX"))
copulaSet(vine2, margins = list(copulaMarginNormal(.25),
    copulaMarginBernoulli(.5)), conditioning = categorical,
  warmStartOnActivate = FALSE, guard = FALSE)
categoricalState <- copulaSnapshot(copulaGet(), variableName = "eta_CL")
categoricalFfem <- copulaFremToFfem(categoricalState, "SEX")
stopifnot(identical(categoricalFfem$representation, "distributional"),
  inherits(try(copulaFfemTransform(categoricalFfem, c(0, 1)), silent = TRUE),
    "try-error"),
  all(is.finite(copulaFfemSimulate(categoricalFfem, c(0, 1), n = 10L,
    seed = 482L))))

copulaClear()
cat("FREM-to-FFEM translation checks passed\n")
