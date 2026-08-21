## Calibrated power of the route-C diagnostic.
## Null threshold is the 95th percentile of the statistic under the GAUSSIAN
## truth -- not a chi-square quantile, because family selection, the parameter
## boundary and within-subject draw correlation all break the chi-square.
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args)) args[1] else "out/power_N100_gumbel.rds"
r <- readRDS(f)
cat(sprintf("file: %s   replicates: gauss=%d alt=%d   N=%d\n\n", f,
            sum(r$arm == "gauss" & r$how == "param"),
            sum(r$arm == "alt" & r$how == "param"), r$N[1]))
cat(sprintf("eta shrinkage: mean=%.3f  max=%.3f\n\n",
            mean(r$shrinkMean), max(r$shrinkMax)))

cat(sprintf("%-7s %-8s %9s %9s %9s %7s %9s\n",
            "PIT", "level", "medNull", "medAlt", "crit95", "power", "sepAUC"))
for (hh in unique(r$how)) for (lv in c("Oracle", "Draws", "Ebe")) {
  col <- paste0("stat", lv)
  a <- r[[col]][r$how == hh & r$arm == "gauss"]
  b <- r[[col]][r$how == hh & r$arm == "alt"]
  crit <- stats::quantile(a, 0.95, names = FALSE)
  pw <- mean(b > crit)
  auc <- mean(outer(b, a, ">")) + 0.5 * mean(outer(b, a, "=="))
  cat(sprintf("%-7s %-8s %9.2f %9.2f %9.2f %7.2f %9.3f\n",
              hh, lv, median(a), median(b), crit, pw, auc))
}

cat("\n-- type-I error if the naive chi-square(6) 95% cutoff (12.59) were used --\n")
for (hh in unique(r$how)) for (lv in c("Oracle", "Draws", "Ebe")) {
  a <- r[[paste0("stat", lv)]][r$how == hh & r$arm == "gauss"]
  cat(sprintf("%-7s %-8s  %.3f\n", hh, lv, mean(a > 12.59)))
}

## NOTE ON THE REFERENCE DISTRIBUTION.  There is no standard one: swapping a
## Gaussian pair copula for Gumbel/Clayton/Frank/Joe changes the FAMILY at equal
## parameter count (0 df), while t adds 1, so the compared models are largely
## non-nested and the "obvious" chi-square(6) has ill-defined df.  That
## ill-definedness IS the finding; the number above shows what taking the obvious
## reference at face value would cost.  The cleanest statement is the next block:
## under a GAUSSIAN truth, how often does AIC select a non-Gaussian pair copula?
cat("\n-- P(AIC selects >=1 non-Gaussian pair copula) --\n")
for (hh in unique(r$how)) for (lv in c("Oracle", "Draws", "Ebe"))
  cat(sprintf("%-7s %-8s  gauss-truth=%.3f   vine-truth=%.3f\n", hh, lv,
      mean(r[[paste0("ng", lv)]][r$how == hh & r$arm == "gauss"] > 0),
      mean(r[[paste0("ng", lv)]][r$how == hh & r$arm == "alt"] > 0)))

cat("\n-- non-Gaussian pair copulas selected (out of 6) --\n")
for (hh in unique(r$how)) for (lv in c("Oracle", "Draws", "Ebe"))
  cat(sprintf("%-7s %-8s gauss-truth=%.2f  vine-truth=%.2f\n", hh, lv,
      mean(r[[paste0("ng", lv)]][r$how == hh & r$arm == "gauss"]),
      mean(r[[paste0("ng", lv)]][r$how == hh & r$arm == "alt"])))
