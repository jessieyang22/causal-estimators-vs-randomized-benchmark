# src/04b_placebo_re74.R
#
# Supplementary placebo check added 2026-10-30 after the primary results (see
# DESIGN.md Deviation 7). It is NOT one of the protocol placebo tests in 3e.
#
# The 3e placebo removes re75 from the covariates, so it tests a model whose only
# earnings control is re74. This check keeps the primary covariate structure
# except the outcome: the placebo outcome is re74 (1974 earnings), and re75 stays
# in the covariates, as it does in the primary model. Estimators: the five
# cross-sectional ones. DiD is not applicable (no pre-period for re74).
# The benchmark is the NSW experimental difference in re74.
# Output: results/placebo_re74.csv
# Run from the repository root: Rscript src/04b_placebo_re74.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
res_dir <- "results"
proc_dir <- file.path("data", "processed")
S <- load_samples()
CORES <- default_cores()

OUT_RE74 <- function(d) d$re74
COVS_RE74PLACEBO <- setdiff(COVS8, "re74")   # re75 stays in

out <- list()
for (label in c("PSID", "CPS")) {
  comp <- if (label == "PSID") S$psid else S$cps
  base <- placebo_cells(label)
  for (cl in base) {
    if (cl$method == "did") next
    cl$id <- sub("_placebo_", "_placebo74_", cl$id)
    cl$obs <- OUT_RE74
    cl$bench <- OUT_RE74
    cl$covs <- COVS_RE74PLACEBO
    cl$grid_type <- "placebo_re74"
    t0 <- proc.time()[3]
    r <- run_cell(cl, S, comp, B = B_BOOT, cores = CORES)
    out[[cl$id]] <- r$row
    saveRDS(r$reps, file.path(proc_dir, paste0("boot_", cl$id, ".rds")))
    message(sprintf("%-24s att=%9.1f bias=%9.1f CI[%9.1f, %9.1f] informative=%s pass=%s (%.0fs)",
                    cl$id, r$row$att, r$row$bias, r$row$bias_lo, r$row$bias_hi,
                    r$row$informative, r$row$pass, proc.time()[3] - t0))
  }
}
tbl <- do.call(rbind, out)
rownames(tbl) <- NULL
write.csv(tbl, file.path(res_dir, "placebo_re74.csv"), row.names = FALSE)
print(tbl[, c("cell_id", "att", "bias", "bias_lo", "bias_hi", "informative", "pass")], digits = 5)
