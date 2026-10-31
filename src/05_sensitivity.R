# src/05_sensitivity.R
#
# Component E: sensitivity analysis for hidden confounding (DESIGN.md 4, E).
#   1. Rosenbaum bounds (Gamma 1 to 3, step 0.25) on the matched pairs from
#      the primary matching estimator (Wilcoxon signed-rank test, rbounds::psens).
#   2. Robustness value (sensemakr) for the pooled regression of re78 on the
#      treatment and the eight covariates (Deviation 4(k): sensemakr needs a
#      treatment coefficient, so this is the pooled OLS, not the comparison-only
#      regression estimate). Includes the value needed to move the estimate to
#      the NSW benchmark.
#   3. Benchmarking: partial R-squared of each observed covariate, to compare
#      with the strength needed to explain the gap. Run in both directions:
#      the default (reduce = TRUE) and upward (reduce = FALSE, Deviation 6,
#      written to results/sensemakr_benchmark_upward.csv).
# Run from the repository root: Rscript src/05_sensitivity.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
suppressMessages({
  library(rbounds)
  library(sensemakr)
})
res_dir <- "results"
S <- load_samples()
nswT <- S$nswT
bench_est <- read.csv(file.path(res_dir, "benchmark.csv"), stringsAsFactors = FALSE)
bench_re78 <- bench_est$est[bench_est$outcome == "re78 (primary benchmark, unadjusted)"]

# Matched pairs taken directly from Matching::Match, so they are exactly the
# pairs behind the matching estimate in src/helpers.R. Match() lists one row per
# (treated, matched control) pair; a treated unit with tied controls has several
# rows. The matched outcome for each treated unit is the mean over its rows, the
# same equal weighting Match() uses (Deviation 4(k)).
match_pairs <- function(yT, yC, lT, lC, caliper = 0.2) {
  nT <- length(yT); nC <- length(yC)
  m <- Match(Y = c(yT, yC), Tr = c(rep(1, nT), rep(0, nC)),
             X = c(lT, lC), M = 1, caliper = caliper, replace = TRUE,
             estimand = "ATT")
  tr <- m$index.treated
  co <- m$index.control - nT
  units <- unique(tr)
  x <- yT[units]
  y <- sapply(units, function(i) mean(yC[co[tr == i]]))
  list(x = x, y = y, dropped = nT - length(units), match_est = m$est)
}

gamma_grid <- seq(1, 3, by = 0.25)
gamma_rows <- list()
rv_rows <- list()
up_rows <- list()

for (g in c("psid", "cps")) {
  label <- if (g == "psid") "PSID" else "CPS"
  comp <- S[[g]]

  # ---- 1. Rosenbaum bounds on the matched pairs ----
  ps <- fit_ps(nswT, comp, COVS8)
  lT <- qlogis(ps$pT); lC <- qlogis(ps$pC)
  mp <- match_pairs(nswT$re78, comp$re78, lT, lC, caliper = 0.2)
  # Sanity check: the pair-based ATT should equal the matching estimator.
  pair_att <- mean(mp$x - mp$y)
  match_att <- estimate("match", nswT$re78, comp$re78, nswT, comp, COVS8, caliper = 0.2)$att
  message(label, ": pair-based ATT = ", round(pair_att, 3),
          " ; Match() ATT = ", round(match_att, 3),
          " ; dropped treated = ", mp$dropped)
  if (abs(pair_att - match_att) > 1e-6) stop(label, ": pair ATT differs from Match() ATT")

  ps_res <- rbounds::psens(mp$x, mp$y, Gamma = max(gamma_grid), GammaInc = 0.25)
  pt <- as.data.frame(ps_res$bounds)
  # psens returns one row per Gamma. Column names: "Gamma", "Lower bound", "Upper bound".
  names(pt) <- make.names(names(pt))
  upper_col <- grep("Upper", names(pt), value = TRUE)[1]
  lower_col <- grep("Lower", names(pt), value = TRUE)[1]
  pt$group <- label
  gamma_rows[[label]] <- data.frame(group = label, Gamma = pt$Gamma,
                                    p_lower = pt[[lower_col]],
                                    p_upper = pt[[upper_col]],
                                    stringsAsFactors = FALSE)
  crossed <- which(pt[[upper_col]] > 0.05)
  first_gamma <- if (length(crossed) > 0) pt$Gamma[crossed[1]] else NA
  message(label, ": matched pairs n = ", length(mp$x),
          " ; first Gamma with upper p > 0.05 = ",
          ifelse(is.na(first_gamma), "none in [1,3]", first_gamma))

  # ---- 2. sensemakr: robustness value, pooled regression ----
  d <- rbind(
    data.frame(nswT[, c("re78", COVS8)], treat = 1),
    data.frame(comp[, c("re78", COVS8)], treat = 0)
  )
  fit <- lm(re78 ~ treat + age + education + black + hispanic + married +
              nodegree + re74 + re75, data = d)
  est_pooled <- unname(coef(fit)["treat"])
  # Moving the estimate to the benchmark. If the benchmark is closer to zero
  # (or past zero), sensemakr reduces the estimate by q; otherwise reduce = FALSE
  # moves it away from zero by q (sensemakr's "reduce" argument).
  q_raw <- (est_pooled - bench_re78) / est_pooled
  reduce_b <- q_raw > 0
  q_bench <- abs(q_raw)
  sens_zero <- sensemakr(model = fit, treatment = "treat", q = 1, alpha = 0.05,
                         benchmark_covariates = "re75")
  sens_bench <- sensemakr(model = fit, treatment = "treat", q = q_bench, alpha = 0.05,
                          reduce = reduce_b, benchmark_covariates = "re75")
  st_zero <- sens_zero$sensitivity_stats
  st_bench <- sens_bench$sensitivity_stats
  rv_rows[[paste(label, "summary")]] <- data.frame(
    group = label, quantity = c("pooled OLS treat coefficient", "NSW benchmark",
                                "fraction of estimate to remove to reach benchmark (negative = increase)",
                                "robustness value to reduce estimate to zero (RV)",
                                "robustness value to reduce estimate to zero, alpha 0.05 (RV_qa)",
                                "robustness value to move estimate to benchmark (RV_q)"),
    value = c(est_pooled, bench_re78, q_raw,
              st_zero$rv_q, st_zero$rv_qa, st_bench$rv_q),
    r2_with_outcome = NA_real_,
    estimate_if_confounder_as_strong = NA_real_,
    stringsAsFactors = FALSE
  )

  # ---- 3. Benchmarking against each observed covariate ----
  for (cv in COVS8) {
    sb <- sensemakr(model = fit, treatment = "treat", q = 1, alpha = 0.05,
                    benchmark_covariates = cv)
    bnd <- sb$bounds
    row <- bnd[bnd$bound_label == paste0("1x ", cv), , drop = FALSE]
    rv_rows[[paste(label, cv)]] <- data.frame(
      group = label,
      quantity = paste0("benchmark covariate: ", cv, " (value = partial R2 with treatment; r2_with_outcome = partial R2 with outcome)"),
      value = row$r2dz.x,
      r2_with_outcome = row$r2yz.dx,
      estimate_if_confounder_as_strong = row$adjusted_estimate,
      stringsAsFactors = FALSE
    )

    # Deviation 6: the loop above only moves the estimate down (sensemakr's
    # default reduce = TRUE). The gap to the benchmark is upward, so the same
    # covariate-strength confounder is also computed in the upward direction.
    sb_up <- sensemakr(model = fit, treatment = "treat", q = 1, alpha = 0.05,
                       reduce = FALSE, benchmark_covariates = cv)
    row_up <- sb_up$bounds[sb_up$bounds$bound_label == paste0("1x ", cv), , drop = FALSE]
    up_rows[[paste(label, cv)]] <- data.frame(
      group = label,
      covariate = cv,
      r2_with_treatment = row_up$r2dz.x,
      r2_with_outcome = row_up$r2yz.dx,
      estimate_if_confounder_as_strong_upward = row_up$adjusted_estimate,
      nsw_benchmark = bench_re78,
      reaches_benchmark = row_up$adjusted_estimate >= bench_re78,
      stringsAsFactors = FALSE
    )
  }
  message(label, ": pooled est = ", round(est_pooled, 2), " ; RV = ", round(st_zero$rv_q, 4),
          " ; RV to benchmark = ", round(st_bench$rv_q, 4), " ; q = ", round(q_raw, 4),
          " ; reduce = ", reduce_b)
}

gam <- do.call(rbind, gamma_rows)
rownames(gam) <- NULL
write.csv(gam, file.path(res_dir, "sensitivity_gamma.csv"), row.names = FALSE)

rv <- do.call(rbind, rv_rows)
rownames(rv) <- NULL
write.csv(rv, file.path(res_dir, "sensemakr_rv.csv"), row.names = FALSE)

up <- do.call(rbind, up_rows)
rownames(up) <- NULL
write.csv(up, file.path(res_dir, "sensemakr_benchmark_upward.csv"), row.names = FALSE)

print(gam)
print(rv)
print(up)
