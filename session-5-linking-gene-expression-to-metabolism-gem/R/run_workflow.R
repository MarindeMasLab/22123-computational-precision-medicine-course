# Run from hands-on_session/: Rscript R/run_workflow.R
source("R/model_io.R")
source("R/ko_workflow.R")
source("R/exercise5_workflow.R")
result <- exercise5_write_outputs()
cat("PASS: Exercise 5 student inputs prepared; live native pair count:",
  result$panel_status$pair_count, "\n")