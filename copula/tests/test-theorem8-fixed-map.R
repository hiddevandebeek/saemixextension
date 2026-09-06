## Executable obligations for the sufficient MAX1--MAX2 certificate used in
## the formal proof.  This tests the algebra and the diagnostic, not MAX1/MAX2
## for every fitted marginal family.

source(file.path("R", "theorem8.R"))

failures <- character()
check <- function(label, condition, detail = "") {
  status <- if (isTRUE(condition)) "PASS" else "FAIL"
  cat(sprintf("%-78s %s %s\n", label, status, detail))
  if (!isTRUE(condition)) failures <<- c(failures, label)
}

H <- matrix(c(4, 1, 1, 2), 2, 2)
eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
m <- min(eigenvalues); M <- max(eigenvalues)
optimum <- c(0.4, -0.8)
Q <- function(theta) -0.5 * sum((theta - optimum) *
  as.numeric(H %*% (theta - optimum)))
gradient <- function(theta) -as.numeric(H %*% (theta - optimum))
step <- 2 / (m + M)
rho <- copulaTheorem8ContractionRate(m, M, step)
map <- function(theta) copulaTheorem8GradientMap(theta, gradient(theta), step)
starts <- rbind(c(-3, 2), c(2, -4), c(0.5, -0.7), c(10, 10))
audit <- copulaTheorem8FixedMapAudit(map, optimum, starts,
  objective = Q, rate = rho, iterations = 40L)

check("T8-1 analytical quadratic rate lies strictly below one",
  rho > 0 && rho < 1, sprintf("rho=%.6g", rho))
check("T8-2 the fixed map is deterministic", audit$deterministic)
check("T8-3 every fixed-map step increases the common objective",
  audit$objectiveNondecreasing)
check("T8-4 every trajectory satisfies the analytical geometric bound",
  audit$certifiedBoundSatisfied,
  sprintf("max observed ratio=%.6g", audit$maximumObservedRatio))
check("T8-5 repeated maps converge to the declared unique maximizer",
  audit$maximumFinalDistance < 1e-9,
  sprintf("max distance=%.3g", audit$maximumFinalDistance))

bad <- try(copulaTheorem8ContractionRate(m, M, 3 / m), silent = TRUE)
check("T8-6 a noncontractive proposed step is rejected",
  inherits(bad, "try-error"))

exact <- copulaTheorem8BackendStatus(TRUE, TRUE, exactMstep = TRUE)
conditional <- copulaTheorem8BackendStatus(TRUE, TRUE, fixedMap = TRUE,
  note = "requires family-specific certificate")
score <- copulaTheorem8BackendStatus(FALSE, FALSE, fixedMap = TRUE)
check("T8-7 an exact finite-statistic M-step reports the M5/Theorem5 route",
  exact$theoremAligned && identical(exact$route, "M5/Theorem5"))
check("T8-8 an uncertified nonlinear fixed map is not labelled theorem-aligned",
  !conditional$theoremAligned &&
    identical(conditional$route, "conditional-until-MAX1-MAX2-certified"))
check("T8-9 a non-finite-statistic backend is directed to augmentation/score-SA",
  !score$theoremAligned &&
    identical(score$route, "augmented-state-or-score-SA"))

cat("\n", length(failures), " failure(s)\n", sep = "")
if (length(failures)) quit(status = 1L)
