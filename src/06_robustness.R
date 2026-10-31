# src/06_robustness.R
#
# Component F: robustness grid (DESIGN.md 3f, with Deviation 4(m)).
# Grid variants (fixed before results), defined in src/helpers.R (grid_cells):
#   - propensity cap 0.05 instead of 0.01 (IPW, AIPW)
#   - matching caliper 0.1 SD (secondary) and no caliper (matching)
#   - covariate set without re74 (regression, matching, IPW, AIPW)
# Every grid cell gets the same bias bootstrap and pass rule as the primary cells.
# The standard-error variants (Abadie-Imbens vs bootstrap, IF vs bootstrap) are
# reported in src/08_figures.R and in the ATT columns, not as separate cells.
#
# Also writes, per estimator, the fragility rule (DESIGN.md 5, item 5): an
# estimator that passes the primary benchmark test but fails in more than half
# of its grid cells is marked "fragile". Also writes the RMSE of the bias across
# the primary and grid cells. The revised matching rule is in src/10_verdicts.R.
# Run from the repository root: Rscript src/06_robustness.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
res_dir <- "results"
proc_dir <- file.path("data", "processed")
S <- load_samples()
CORES <- default_cores()

rows <- list()
for (g in c("psid", "cps")) {
  label <- if (g == "psid") "PSID" else "CPS"
  comp <- S[[g]]
  for (cl in grid_cells(label)) {
    cl$group <- label
    t0 <- proc.time()[3]
    r <- run_cell(cl, S, comp, B = B_BOOT, cores = CORES)
    rows[[cl$id]] <- r$row
    saveRDS(r$reps, file.path(proc_dir, paste0("boot_", cl$id, ".rds")))
    message(sprintf("%-22s att=%9.1f bias=%9.1f CI[%9.1f, %9.1f] informative=%s pass=%s (%.0fs)",
                    cl$id, r$row$att, r$row$bias, r$row$bias_lo, r$row$bias_hi,
                    r$row$informative, r$row$pass, proc.time()[3] - t0))
  }
}
grid <- do.call(rbind, rows)
rownames(grid) <- NULL

# Primary cells (from 02) are needed for the fragility rule and the RMSE.
prim <- read.csv(file.path(res_dir, "observational_vs_benchmark.csv"), stringsAsFactors = FALSE)
prim$grid_type <- "primary"
both <- rbind(prim[, names(grid)], grid)

summ <- do.call(rbind, lapply(unique(both$method), function(m) {
  do.call(rbind, lapply(c("PSID", "CPS"), function(gk) {
    p <- both[both$group == gk & both$method == m & both$grid_type == "primary", ]
    gg <- both[both$group == gk & both$method == m & both$grid_type != "primary", ]
    if (nrow(p) == 0) return(NULL)
    n_grid <- nrow(gg)
    n_fail <- sum(!gg$pass)
    data.frame(
      group = gk, method = m, estimator = p$estimator[1],
      primary_pass = p$pass[1],
      n_grid_cells = n_grid,
      n_grid_fail = n_fail,
      grid_fail_share = if (n_grid > 0) n_fail / n_grid else NA_real_,
      fragile = if (n_grid > 0) (p$pass[1] && n_fail > n_grid / 2) else NA,
      rmse_bias_all_cells = sqrt(mean(both$bias[both$group == gk & both$method == m]^2)),
      rmse_bias_grid_only = if (n_grid > 0) sqrt(mean(gg$bias^2)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(summ) <- NULL

write.csv(grid, file.path(res_dir, "robustness_grid.csv"), row.names = FALSE)
write.csv(summ, file.path(res_dir, "robustness_summary.csv"), row.names = FALSE)
print(grid[, c("cell_id", "grid_type", "att", "bias", "bias_lo", "bias_hi", "informative", "pass")], digits = 5)
print(summ, digits = 4)
