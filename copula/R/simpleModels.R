## A ladder of increasingly demanding models, so that "where does it break" can
## be attributed to a CAUSE rather than just observed.  The two candidate causes
## are confounded in the 2-cmt model used so far: dimension (d) and
## identifiability (shrinkage).  This varies them as independently as possible.
##
##   iv1   d=2  V, CL          1-cmt IV bolus, rich profile   -- both strong
##   oral1 d=3  ka, V, CL      1-cmt oral                     -- all three strong
##   iv2   d=4  V1, CL, Q, V2  2-cmt IV bolus                 -- Q, V2 weak
##
## Everything else is held fixed: eta SDs, pairwise correlation structure,
## residual error, number of subjects, observations per subject.

MODELS <- list(

  iv1 = list(
    d = 2, names = c("V", "CL"), true = c(V = 10, CL = 3),
    times = c(0.25, .5, 1, 2, 3, 4, 6, 8, 10, 12), dose = 100,
    f = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      V <- psi[id, 1]; CL <- psi[id, 2]
      dose / V * exp(-CL / V * tim)
    }),

  oral1 = list(
    d = 3, names = c("ka", "V", "CL"), true = c(ka = 1.5, V = 10, CL = 3),
    times = c(0.25, .5, 1, 1.5, 2, 3, 4, 6, 8, 12), dose = 100,
    f = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      ka <- psi[id, 1]; V <- psi[id, 2]; CL <- psi[id, 3]
      k <- CL / V
      dose * ka / (V * (ka - k)) * (exp(-k * tim) - exp(-ka * tim))
    }),

  iv2 = list(
    d = 4, names = c("V1", "CL", "Q", "V2"), true = c(V1 = 10, CL = 3, Q = 2, V2 = 30),
    times = c(0.25, .5, 1, 2, 4, 6, 8, 12, 16, 24), dose = 100,
    f = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      V1 <- psi[id, 1]; CL <- psi[id, 2]; Q <- psi[id, 3]; V2 <- psi[id, 4]
      k10 <- CL / V1; k12 <- Q / V1; k21 <- Q / V2
      s <- k10 + k12 + k21
      disc <- sqrt(pmax(s * s - 4 * k10 * k21, 0))
      al <- (s + disc) / 2; be <- (s - disc) / 2
      A <- dose / V1 * (al - k21) / (al - be)
      B <- dose / V1 * (k21 - be) / (al - be)
      A * exp(-al * tim) + B * exp(-be * tim)
    })
)

## Same eta SD for every parameter, and an AR(1)-style correlation, so the
## dependence structure is comparable across d.
ETA_SD  <- 0.30
ETA_RHO <- 0.55

etaR <- function(d) ETA_RHO^abs(outer(seq_len(d), seq_len(d), "-"))

simModel <- function(m, N, vine, sdEta = rep(ETA_SD, m$d), propErr = 0.10) {
  eta <- rEtaVine(N, vine, sdEta)
  psi <- sweep(exp(eta), 2, m$true, "*"); colnames(psi) <- m$names
  nt <- length(m$times)
  dd <- data.frame(id = rep(seq_len(N), each = nt), dose = m$dose,
                   time = rep(m$times, N))
  f <- m$f(psi, dd$id, cbind(dd$dose, dd$time))
  dd$y <- f * (1 + propErr * rnorm(nrow(dd)))
  dd$y[dd$y <= 0] <- 1e-6
  list(data = dd, eta = eta, psi = psi)
}

saemixModelFor <- function(m) saemix::saemixModel(
  model = m$f, modeltype = "structural", description = "",
  psi0 = matrix(m$true, ncol = m$d, byrow = TRUE, dimnames = list(NULL, m$names)),
  transform.par = rep(1, m$d), covariance.model = matrix(1, m$d, m$d),
  omega.init = diag(rep(ETA_SD^2, m$d)), error.model = "proportional", verbose = FALSE)

saemixDataFor <- function(dd) saemix::saemixData(name.data = dd, header = TRUE,
  name.group = "id", name.predictors = c("dose", "time"),
  name.response = "y", verbose = FALSE)
