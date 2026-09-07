## Summarise the test-first study with exactly the calculations the published
## summariser uses, then add the counts that only the test-first workflow has.
##
##   Rscript summarize_testfirst.R INPUT_DIR OUTPUT_DIR
##
## summarize.R fixes its own input and output from a pilot/final switch; its
## body is reused verbatim with those header lines removed, so the endpoint
## tables of the two studies are produced by the same code.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: summarize_testfirst.R INPUT OUTPUT")
input <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
output <- args[2L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)

base <- Sys.getenv("SUMMARIZE_BASE", "summarize.R")
if (!file.exists(base)) stop("cannot find summarize.R at ", base)
src <- readLines(base)
header <- grep("^(args <-|mode <-|if \\(!mode|root <-|input <-|output <-|dir.create\\(output)",
  src)
stopifnot(length(header) >= 6L)
eval(parse(text = src[-header]))

## Test-first specific record: what the score test decided, how often the
## ranking was consulted at all, and whether it agreed with the test.
tests <- do.call(rbind, lapply(results, function(r) {
  table <- r$shapeTest
  table <- table[table$test %in% c("He3", "He4"), ]
  data.frame(replicate = r$replicate, coordinate = table$coordinate,
    test = table$test, statistic = table$statistic.opg, score = table$score,
    p_value = table$p.value.boot)
}))
write.csv(tests, file.path(output, "shape_test_table.csv"), row.names = FALSE)

level <- results[[1L]]$shapeTestLevel %||% .05
rejectFrame <- do.call(rbind, lapply(results, function(r) data.frame(
  replicate = r$replicate,
  volume = isTRUE(r$shapeTestRejects[1L]),
  clearance = isTRUE(r$shapeTestRejects[2L]),
  ranking_run = isTRUE(r$rankingRun))))
write.csv(rejectFrame, file.path(output, "shape_test_decisions.csv"),
  row.names = FALSE)

flexible <- summary[summary$arm == "Flexible FREM", ]
flexible <- flexible[order(flexible$replicate), ]
rejectFrame <- rejectFrame[order(rejectFrame$replicate), ]
stopifnot(identical(flexible$replicate, rejectFrame$replicate))

quality$test_level <- level
quality$clearance_rejected <- sum(rejectFrame$clearance)
quality$volume_rejected <- sum(rejectFrame$volume)
quality$either_rejected <- sum(rejectFrame$ranking_run)
quality$neither_rejected <- sum(!rejectFrame$ranking_run)
## Of the datasets where the test opened the clearance coordinate, how often
## the ranking then chose a non-lognormal family for it.
quality$clearance_changed_after_reject <-
  sum(rejectFrame$clearance & flexible$familyCL != "lognormal")
quality$clearance_changed_without_reject <-
  sum(!rejectFrame$clearance & flexible$familyCL != "lognormal")
write.csv(quality, file.path(output, "quality_summary.csv"), row.names = FALSE)

byDirection <- aggregate(p_value ~ coordinate + test, data = tests,
  FUN = function(p) sum(p < level))
names(byDirection)[3L] <- "rejections"
byDirection$datasets <- length(results)
write.csv(byDirection, file.path(output, "shape_test_rejections.csv"),
  row.names = FALSE)
print(quality)
print(byDirection)
