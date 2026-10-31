# src/09_seed_sensitivity.R
#
# Seed sensitivity of the bias bootstrap (review item, 2026-10-30; Deviation 7).
# The primary results use seed SEED = 20261029 (helpers.R). This script re-runs
# every NON-matching cell (naive, regression, IPW, AIPW, DiD; primary, placebo,
# pre-trend and robustness grid) under two more seeds, with the same B = 2000.
# Matching cells are not re-run: their verdicts use the Abadie-Imbens SE
# (src/10_verdicts.R), which does not depend on the bootstrap seed.
#
# Output: results/seed_sensitivity.csv, one row per (cell, seed). Rows are
# appended as each cell finishes, so an interrupted run resumes where it stopped.
# Primary files (boot_*.rds, the primary CSVs) are not touched.
# Run from the repository root: Rscript src/09_seed_sensitivity.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
S <- load_samples()
CORES <- default_cores()
OUT <- file.path("results", "seed_sensitivity.csv")
EXTRA_SEEDS <- c(20361029, 20461029)

all_cells <- list()
for (label in c("PSID", "CPS")) {
  cl <- c(primary_cells(label), placebo_cells(label), list(pretrend_cell(label)), grid_cells(label))
  cl <- lapply(cl, function(x) { x$group <- label; x })
  all_cells <- c(all_cells, cl)
}
all_cells <- Filter(function(x) x$method != "match", all_cells)
message("Non-matching cells per seed: ", length(all_cells))

done <- character(0)
if (file.exists(OUT)) {
  prev <- read.csv(OUT, stringsAsFactors = FALSE)
  done <- paste(prev$cell_id, prev$seed)
}

for (sd in EXTRA_SEEDS) {
  SEED <- sd   # run_cell() reads SEED from the global environment
  for (cl in all_cells) {
    key <- paste(cl$id, sd)
    if (key %in% done) next
    t0 <- proc.time()[3]
    r <- run_cell(cl, S, S[[if (cl$group == "PSID") "psid" else "cps"]], B = B_BOOT, cores = CORES)
    row <- r$row
    row$seed <- sd
    first <- !file.exists(OUT)
    write.table(row, OUT, sep = ",", row.names = FALSE, col.names = first,
                append = !first, qmethod = "double")
    message(sprintf("seed %d %-22s bias=%9.1f CI[%9.1f, %9.1f] pass=%s (%.0fs)",
                    sd, cl$id, row$bias, row$bias_lo, row$bias_hi, row$pass,
                    proc.time()[3] - t0))
  }
}
message("done")
