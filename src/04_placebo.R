# src/04_placebo.R
#
# Component D: placebo and pre-trend tests (DESIGN.md 3e, with Deviation 4(c)).
#   - Placebo outcome re75 (1975 earnings, treated as pre-program) for the five
#     cross-sectional estimators. re75 is removed from the covariate set.
#   - Pre-trend test for difference-in-differences: change from 1974 to 1975
#     (re75 minus re74), against the same change in the NSW experiment.
# It also writes the summary-claim check (DESIGN.md 3d), which needs both the
# benchmark result (02) and these placebo results.
# Run from the repository root: Rscript src/04_placebo.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
res_dir <- "results"
proc_dir <- file.path("data", "processed")
S <- load_samples()
CORES <- 8L

placebo_out <- list()
pretrend_out <- list()

for (g in c("psid", "cps")) {
  label <- if (g == "psid") "PSID" else "CPS"
  comp <- S[[g]]

  # Placebo: re75 outcome, covariates without re75.
  pl_cells <- list(
    list(id = paste0(label, "_placebo_naive"), estimator = "1 Naive difference",
         method = "naive", obs = OUT_RE75, bench = OUT_RE75, covs = COVS_PLACEBO,
         trim = 0.01, caliper = 0.2, grid_type = "placebo", group = label),
    list(id = paste0(label, "_placebo_reg"), estimator = "2 Regression adjustment",
         method = "reg", obs = OUT_RE75, bench = OUT_RE75, covs = COVS_PLACEBO,
         trim = 0.01, caliper = 0.2, grid_type = "placebo", group = label),
    list(id = paste0(label, "_placebo_match"), estimator = "3 Nearest-neighbour matching",
         method = "match", obs = OUT_RE75, bench = OUT_RE75, covs = COVS_PLACEBO,
         trim = 0.01, caliper = 0.2, grid_type = "placebo", group = label),
    list(id = paste0(label, "_placebo_ipw"), estimator = "4 Inverse probability weighting",
         method = "ipw", obs = OUT_RE75, bench = OUT_RE75, covs = COVS_PLACEBO,
         trim = 0.01, caliper = 0.2, grid_type = "placebo", group = label),
    list(id = paste0(label, "_placebo_aipw"), estimator = "5 Doubly robust (AIPW)",
         method = "aipw", obs = OUT_RE75, bench = OUT_RE75, covs = COVS_PLACEBO,
         trim = 0.01, caliper = 0.2, grid_type = "placebo", group = label)
  )
  for (cl in pl_cells) {
    t0 <- proc.time()[3]
    r <- run_cell(cl, S, comp, B = B_BOOT, cores = CORES)
    placebo_out[[cl$id]] <- r$row
    saveRDS(r$reps, file.path(proc_dir, paste0("boot_", cl$id, ".rds")))
    message(sprintf("%-22s att=%9.1f bias=%9.1f CI[%9.1f, %9.1f] informative=%s pass=%s (%.0fs)",
                    cl$id, r$row$att, r$row$bias, r$row$bias_lo, r$row$bias_hi,
                    r$row$informative, r$row$pass, proc.time()[3] - t0))
  }

  # Pre-trend test for DiD: change 1974 to 1975, unadjusted (DiD has no covariates).
  pt_cell <- list(id = paste0(label, "_pretrend_did"), estimator = "6 Difference-in-differences",
                  method = "did", obs = OUT_CHG75_74, bench = OUT_CHG75_74,
                  covs = COVS8, trim = 0.01, caliper = 0.2, grid_type = "pretrend", group = label)
  t0 <- proc.time()[3]
  r <- run_cell(pt_cell, S, comp, B = B_BOOT, cores = CORES)
  pretrend_out[[pt_cell$id]] <- r$row
  saveRDS(r$reps, file.path(proc_dir, paste0("boot_", pt_cell$id, ".rds")))
  message(sprintf("%-22s att=%9.1f bias=%9.1f CI[%9.1f, %9.1f] informative=%s pass=%s (%.0fs)",
                  pt_cell$id, r$row$att, r$row$bias, r$row$bias_lo, r$row$bias_hi,
                  r$row$informative, r$row$pass, proc.time()[3] - t0))
}

placebo_tbl <- do.call(rbind, placebo_out)
rownames(placebo_tbl) <- NULL
write.csv(placebo_tbl, file.path(res_dir, "placebo_re75.csv"), row.names = FALSE)

pretrend_tbl <- do.call(rbind, pretrend_out)
rownames(pretrend_tbl) <- NULL
write.csv(pretrend_tbl, file.path(res_dir, "pretrend_did.csv"), row.names = FALSE)

# ---- Summary-claim check (DESIGN.md 3d) ---------------------------------------
# Claim allowed only if, for at least one comparison group, at least 4 of the
# 6 estimators pass the benchmark test AND pass their placebo test.
bench_tbl <- read.csv(file.path(res_dir, "observational_vs_benchmark.csv"), stringsAsFactors = FALSE)
plac_all <- rbind(placebo_tbl, pretrend_tbl)
plac_all$method_key <- sub(".*_", "", plac_all$cell_id)   # naive, reg, match, ipw, aipw, did
plac_all$method_key <- ifelse(plac_all$grid_type == "pretrend", "did", plac_all$method_key)
plac_all$group_key <- ifelse(grepl("^PSID", plac_all$cell_id), "PSID", "CPS")
bench_tbl$method_key <- bench_tbl$method
bench_tbl$group_key <- bench_tbl$group

summ <- do.call(rbind, lapply(c("PSID", "CPS"), function(gk) {
  b <- bench_tbl[bench_tbl$group_key == gk, c("method_key", "pass", "informative")]
  p <- plac_all[plac_all$group_key == gk, c("method_key", "pass", "informative")]
  names(p) <- c("method_key", "placebo_pass", "placebo_informative")
  m <- merge(b, p, by = "method_key", sort = FALSE)
  m$both_pass <- m$pass & m$placebo_pass
  data.frame(group = gk,
             n_estimators = nrow(m),
             n_pass_benchmark = sum(m$pass),
             n_pass_benchmark_and_placebo = sum(m$both_pass),
             claim_met = sum(m$both_pass) >= 4,
             passing_estimators = paste(m$method_key[m$both_pass], collapse = ";"),
             stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(res_dir, "summary_claim.csv"), row.names = FALSE)

print(placebo_tbl[, c("cell_id", "att", "bias", "bias_lo", "bias_hi", "informative", "pass")], digits = 5)
print(pretrend_tbl[, c("cell_id", "att", "bias", "bias_lo", "bias_hi", "informative", "pass")], digits = 5)
print(summ)
