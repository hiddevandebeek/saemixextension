suppressPackageStartupMessages({
  library(devtools)
  copulaTestRepo <- normalizePath(".", mustWork = TRUE)
  load_all(copulaTestRepo, quiet = TRUE)
})
