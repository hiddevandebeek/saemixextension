## Profiling is opt-in, append-only, and numerically inert.

suppressMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-60s %s %s\n", lbl,
    if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

v <- rvinecopulib::vinecop_dist(
  list(list(rvinecopulib::bicop_dist("clayton", parameters = 1.2))),
  rvinecopulib::dvine_structure(1:2))
profilePath <- tempfile(fileext = ".csv")
copulaSet(v, sd = c(.3, .4), freezeSd = TRUE, freezeVine = TRUE,
  profileFile = profilePath, profileEvery = 2L)
copulaProfileEvent(1L, "iteration_start", 1)
copulaProfileEvent(2L, "iteration_start", .5)
copulaProfileEvent(2L, "iteration_end", .5, .01)
p <- utils::read.csv(profilePath, stringsAsFactors = FALSE)
ok("Profiler respects profileEvery", nrow(p) == 2L)
ok("Profiler writes stable iteration/event fields",
   identical(p$iteration, c(2L, 2L)) &&
     identical(p$event, c("iteration_start", "iteration_end")))
ok("Automatic core resolver is conservative",
   copulaResolveCores("auto") >= 1L && copulaResolveCores("auto") <= 4L)

unlink(profilePath)
copulaClear()
cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
