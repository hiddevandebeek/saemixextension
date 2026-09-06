suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

path <- "C:/package/saemix-copula/copula/examples/gamma-tail-vpc/out/example_fit.rds"
if (!file.exists(path)) stop("fitted Gamma-tail example is unavailable")
x <- readRDS(path)

standard <- copulaFremToFfem(x$standardFit, "CRP")
flexible <- copulaFremToFfem(x$flexibleFit, "CRP")
stopifnot(isTRUE(standard$linearAvailable), !is.null(standard$native),
  isTRUE(flexible$linearAvailable), is.null(flexible$native),
  identical(flexible$covariates, "CRP"))

crp <- matrix(qgamma(c(.95, .99, .999), x$truth$gammaShape,
  scale = x$truth$gammaScale), ncol = 1L, dimnames = list(NULL, "CRP"))
z <- copulaFfemTransform(flexible, crp)
expected <- copulaFfemLocation(flexible, crp)
stopifnot(max(abs(expected - z %*% t(flexible$coefficient))) < 1e-12)

fullStandard <- cbind(CRP = crp[, 1L], ALB = NA_real_)
standardConditional <- copulaGaussianFremConditional(fullStandard,
  x$standardState$vine, x$standardState$margins, x$standardState$dEta)
standardDirect <- standardConditional$mean %*%
  diag(copulaMarginScales(x$standardState$margins)[1:2], 2L)
standardScoreFfem <- copulaFfemLocation(standard, crp)
standardRawFfem <- sweep(crp, 2L, standard$native$centre, "-") %*%
  t(standard$native$coefficient)
stopifnot(max(abs(standardDirect - standardScoreFfem)) < 1e-10,
  max(abs(standardDirect - standardRawFfem)) < 1e-12)

set.seed(5001L)
translated <- copulaFfemSimulate(flexible, crp, n = 2000L)
full <- cbind(CRP = crp[rep(seq_len(nrow(crp)), each = 2000L), 1L],
  ALB = NA_real_)
set.seed(5001L)
direct <- copulaGaussianFremRandEta(full, x$flexibleState$vine,
  x$flexibleState$margins, x$flexibleState$dEta)
stopifnot(max(abs(translated - direct)) == 0)

sourceRow <- attr(translated, "sourceRow")
for (i in seq_len(nrow(crp))) {
  local <- translated[sourceRow == i, , drop = FALSE]
  stopifnot(max(abs(colMeans(local) - expected[i, ])) < .02,
    max(abs(stats::cov(local) - flexible$omega)) < .012)
}

times <- c(.5, 2, 8)
predict_response <- function(eta) {
  V <- x$flexibleFit@results@fixed.effects[1L] * exp(eta[, 1L])
  CL <- x$flexibleFit@results@fixed.effects[2L] * exp(eta[, 2L])
  vapply(times, function(time) 100 / V * exp(-(CL / V) * time),
    numeric(nrow(eta)))
}
stopifnot(max(abs(predict_response(translated) - predict_response(direct))) == 0)

copulaClear()
cat("fitted FREM-to-FFEM translation checks passed\n")
