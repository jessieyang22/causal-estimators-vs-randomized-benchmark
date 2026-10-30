# src/02_observational.R
#
# Component B: six observational estimators for each comparison group (PSID,
# CPS), each compared with the NSW experimental benchmark (DESIGN.md 3b to 3d).
# Run from the repository root: Rscript src/02_observational.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
res_dir <- "results"
proc_dir <- file.path("data", "processed")
S <- load_samples()
CORES <- 8L

out <- list()
for (g in c("psid", "cps")) {
  label <- if (g == "psid") "PSID" else "CPS"
  comp <- S[[g]]
  cells <- tag(primary_cells(label), label)
  for (cl in cells) {
    t0 <- proc.time()[3]
    r <- run_cell(cl, S, comp, B = B_BOOT, cores = CORES)
    out[[cl$id]] <- r$row
    saveRDS(r$reps, file.path(proc_dir, paste0("boot_", cl$id, ".rds")))
    message(sprintf("%-18s att=%9.1f bias=%9.1f CI[%9.1f, %9.1f] informative=%s pass=%s (%.0fs)",
                    cl$id, r$row$att, r$row$bias, r$row$bias_lo, r$row$bias_hi,
                    r$row$informative, r$row$pass, proc.time()[3] - t0))
  }
}

tbl <- do.call(rbind, out)
rownames(tbl) <- NULL
write.csv(tbl, file.path(res_dir, "observational_vs_benchmark.csv"), row.names = FALSE)
print(tbl[, c("cell_id", "att", "att_lo", "att_hi", "bench", "bias", "bias_pct",
              "bias_lo", "bias_hi", "H_o", "H_b", "informative", "pass",
              "mean_dropped_treated")], digits = 5)
