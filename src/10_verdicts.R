# src/10_verdicts.R
#
# Revised verdicts (review items 1 and 3; DESIGN.md Deviation 7, 2026-10-30).
# This script only reads committed or generated CSV files. It does not re-run
# any bootstrap. Its outputs are the ones the README reports.
#
#   1. Matching rule. For matching cells, the bias interval is built from the
#      Abadie-Imbens SE of the ATT (att_se in the primary CSVs), not from the
#      bootstrap percentiles, because the bootstrap is not valid for matching
#      with replacement (Abadie and Imbens 2008). The half-width is
#      1.96 * sqrt(att_se^2 + bench_se^2). This is conservative: the treated mean
#      cancels in the bias, but the two SEs are added as if independent, so
#      the shared treated-group noise is counted twice. Pass = informative
#      (half-width <= 2 * H_b) and interval contains zero. Non-matching cells keep
#      the bootstrap rule of DESIGN.md 3d.
#   2. Seed sensitivity. results/seed_sensitivity.csv holds the non-matching cells
#      re-run under two extra seeds (src/09). Each seed's verdicts are reported
#      separately; the primary seed is 20261029.
#
# Outputs:
#   results/verdicts_revised.csv          one row per cell, pass under each seed
#   results/summary_claim_revised.csv     DESIGN.md 3d claim check, per seed
#   results/robustness_summary_revised.csv fragility rule, per seed
# Run from the repository root: Rscript src/10_verdicts.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
res <- "results"
rd <- function(f) read.csv(file.path(res, f), stringsAsFactors = FALSE)
PRIMARY_SEED <- 20261029
EXTRA <- c(20361029, 20461029)

prim <- rd("observational_vs_benchmark.csv")
plac <- rd("placebo_re75.csv")
pret <- rd("pretrend_did.csv")
grid <- rd("robustness_grid.csv")
p74  <- if (file.exists(file.path(res, "placebo_re74.csv"))) rd("placebo_re74.csv") else NULL
seedf <- if (file.exists(file.path(res, "seed_sensitivity.csv"))) rd("seed_sensitivity.csv") else NULL

base_cols <- c("cell_id", "group", "grid_type", "method", "estimator", "att", "bias",
               "bias_lo", "bias_hi", "H_b", "att_se_method", "att_se", "bench_se")
all <- rbind(prim[, base_cols], plac[, base_cols], pret[, base_cols],
             grid[, base_cols])
all$boot_pass_seed1 <- c(prim$pass, plac$pass, pret$pass, grid$pass)
all$boot_lo <- all$bias_lo
all$boot_hi <- all$bias_hi
all$is_match <- all$method == "match"

# Revised rule for matching: Abadie-Imbens SE, conservative quadrature sum.
half_ai <- 1.96 * sqrt(all$att_se^2 + all$bench_se^2)
all$rule_src <- ifelse(all$is_match, "abadie_imbens", "bootstrap")
all$rule_lo <- ifelse(all$is_match, all$bias - half_ai, all$boot_lo)
all$rule_hi <- ifelse(all$is_match, all$bias + half_ai, all$boot_hi)
all$H_o_rule <- ifelse(all$is_match, half_ai, (all$boot_hi - all$boot_lo) / 2)
all$informative_rule <- all$H_o_rule <= 2 * all$H_b
all$pass_rule <- all$informative_rule & all$rule_lo <= 0 & all$rule_hi >= 0
all$half_over_limit <- all$H_o_rule / (2 * all$H_b)

# Seed columns. Primary seed uses the primary run; matching cells do not depend on
# the bootstrap seed and carry their (single) verdict under every seed.
seed_pass <- function(sd) {
  if (sd == PRIMARY_SEED) return(all$pass_rule)
  if (is.null(seedf)) return(rep(NA, nrow(all)))
  s <- seedf[seedf$seed == sd, ]
  idx <- match(all$cell_id, s$cell_id)
  out <- ifelse(all$is_match, all$pass_rule, s$pass[idx])
  out
}
seeds_all <- c(PRIMARY_SEED, EXTRA)
for (sd in seeds_all) all[[paste0("pass_seed_", sd)]] <- as.logical(seed_pass(sd))
all$seeds_complete <- rowSums(!is.na(all[, paste0("pass_seed_", seeds_all)])) == length(seeds_all)
all$seed_flip <- apply(all[, paste0("pass_seed_", seeds_all)], 1, function(x) {
  x <- x[!is.na(x)]; length(unique(x)) > 1 })

verdicts <- all[, c("cell_id", "group", "grid_type", "method", "estimator", "bias",
                    "boot_lo", "boot_hi", "rule_src", "rule_lo", "rule_hi",
                    "H_o_rule", "H_b", "informative_rule", "pass_rule",
                    "boot_pass_seed1", paste0("pass_seed_", seeds_all),
                    "seeds_complete", "seed_flip")]
verdicts$boot_pass_seed1 <- as.logical(verdicts$boot_pass_seed1)
write.csv(verdicts, file.path(res, "verdicts_revised.csv"), row.names = FALSE)

# ---- Summary claim (DESIGN.md 3d) under the revised rule, per seed ----------
claim_rows <- list()
for (sd in seeds_all) {
  pc <- paste0("pass_seed_", sd)
  for (gk in c("PSID", "CPS")) {
    bm <- all[all$group == gk & all$grid_type == "primary", ]
    pl <- all[all$group == gk & all$grid_type %in% c("placebo", "pretrend"), ]
    bm_pass <- as.logical(bm[[pc]])
    # placebo rows: the method key is the last token of the cell id; the DiD
    # pre-trend row (grid_type "pretrend") is the DiD key.
    key <- ifelse(pl$grid_type == "pretrend", "did", sub(".*_", "", pl$cell_id))
    pl_pass <- as.logical(pl[[pc]])[match(bm$method, key)]
    both <- bm_pass & pl_pass
    claim_rows[[paste(sd, gk)]] <- data.frame(
      seed = sd, group = gk,
      n_pass_benchmark = sum(bm_pass, na.rm = TRUE),
      n_pass_placebo = sum(pl_pass, na.rm = TRUE),
      n_pass_both = sum(both, na.rm = TRUE),
      claim_met = sum(both, na.rm = TRUE) >= 4,
      passing_both = paste(bm$method[which(both)], collapse = ";"),
      benchmark_passing = paste(bm$method[which(bm_pass)], collapse = ";"),
      stringsAsFactors = FALSE)
  }
}
claim <- do.call(rbind, claim_rows); rownames(claim) <- NULL
write.csv(claim, file.path(res, "summary_claim_revised.csv"), row.names = FALSE)

# ---- Fragility (DESIGN.md 5, item 5) under the revised rule, per seed --------
frag_rows <- list()
for (sd in seeds_all) {
  pc <- paste0("pass_seed_", sd)
  for (gk in c("PSID", "CPS")) {
    for (m in c("naive", "reg", "match", "ipw", "aipw", "did")) {
      p <- all[all$group == gk & all$method == m & all$grid_type == "primary", ]
      gg <- all[all$group == gk & all$method == m & !(all$grid_type %in% c("primary", "placebo", "pretrend")), ]
      if (nrow(p) == 0) next
      n_grid <- nrow(gg)
      gpass <- as.logical(gg[[pc]])
      n_fail <- sum(!gpass, na.rm = TRUE)
      ppass <- as.logical(p[[pc]])
      frag_rows[[paste(sd, gk, m)]] <- data.frame(
        seed = sd, group = gk, method = m,
        primary_pass = ppass, n_grid_cells = n_grid, n_grid_fail = n_fail,
        fragile = if (n_grid > 0) isTRUE(ppass) && n_fail > n_grid / 2 else NA,
        stringsAsFactors = FALSE)
    }
  }
}
frag <- do.call(rbind, frag_rows); rownames(frag) <- NULL
write.csv(frag, file.path(res, "robustness_summary_revised.csv"), row.names = FALSE)

# ---- Printed summary --------------------------------------------------------
cat("\nVerdict count (primary-seed and extra-seed pass columns, pooled):\n")
cat("  primary cells (benchmark):", nrow(all[all$grid_type == "primary", ]),
    "| placebo re75:", nrow(all[all$grid_type == "placebo", ]),
    "| pre-trend:", nrow(all[all$grid_type == "pretrend", ]),
    "| grid:", sum(!(all$grid_type %in% c("primary", "placebo", "pretrend"))), "\n")
print(verdicts[, c("cell_id", "bias", "boot_lo", "boot_hi", "rule_src", "rule_lo", "rule_hi",
                   "H_o_rule", "H_b", "pass_rule", paste0("pass_seed_", seeds_all), "seed_flip")],
      digits = 5)
print(claim, digits = 4)
print(frag, digits = 4)
