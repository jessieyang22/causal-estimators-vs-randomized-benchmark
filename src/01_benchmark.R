# src/01_benchmark.R
#
# Component A: benchmark construction and randomization check (DESIGN.md 3a).
# Run from the repository root: Rscript src/01_benchmark.R
#
# Uses pre-treatment covariates for the balance check, and re78 / re75 / re74
# for the benchmark estimates. Nothing here is a comparison with an observational
# estimator.

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
res_dir <- "results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

S <- load_samples()
nswT <- S$nswT
nswN <- S$nswN
nsw <- rbind(nswT, nswN)
cat("NSW treated:", nrow(nswT), " NSW control:", nrow(nswN), "\n")

# ---- A1. Benchmark (primary): unadjusted difference in mean re78 -----------
bench_re78 <- neyman(nswT$re78, nswN$re78)

# ---- A2. Secondary: regression adjusted, experimental sample only ----------
fit_adj <- lm(re78 ~ treat + age + education + black + hispanic + married +
                nodegree + re74 + re75, data = nsw)
cf <- summary(fit_adj)$coefficients["treat", ]
adj_est <- unname(cf["Estimate"]); adj_se <- unname(cf["Std. Error"])

# ---- A3. Permutation p-value for the primary benchmark (10,000 re-assignments)
set.seed(SEED)
B_PERM <- 10000
obs_diff <- bench_re78["est"]
n1 <- nrow(nswT)
perm_diff <- replicate(B_PERM, {
  idx <- sample.int(nrow(nsw), n1)
  mean(nsw$re78[idx]) - mean(nsw$re78[-idx])
})
perm_p <- (1 + sum(abs(perm_diff) >= abs(obs_diff))) / (1 + B_PERM)

# ---- A4. Placebo benchmark (re75) and pre-trend benchmark (re75 - re74) -----
bench_re75 <- neyman(nswT$re75, nswN$re75)
bench_pre <- neyman(OUT_CHG75_74(nswT), OUT_CHG75_74(nswN))

bench_tbl <- rbind(
  data.frame(outcome = "re78 (primary benchmark, unadjusted)",
             estimator = "NSW experimental difference in means",
             est = bench_re78["est"], se = bench_re78["se"],
             lo = bench_re78["lo"], hi = bench_re78["hi"],
             n_treated = nrow(nswT), n_control = nrow(nswN),
             perm_p = perm_p, stringsAsFactors = FALSE),
  data.frame(outcome = "re78 (secondary, regression adjusted)",
             estimator = "OLS coefficient on treat, 8 covariates",
             est = adj_est, se = adj_se,
             lo = adj_est - 1.96 * adj_se, hi = adj_est + 1.96 * adj_se,
             n_treated = nrow(nswT), n_control = nrow(nswN),
             perm_p = NA, stringsAsFactors = FALSE),
  data.frame(outcome = "re75 (placebo benchmark)",
             estimator = "NSW experimental difference in means",
             est = bench_re75["est"], se = bench_re75["se"],
             lo = bench_re75["lo"], hi = bench_re75["hi"],
             n_treated = nrow(nswT), n_control = nrow(nswN),
             perm_p = NA, stringsAsFactors = FALSE),
  data.frame(outcome = "re75 minus re74 (pre-trend benchmark)",
             estimator = "NSW experimental difference in change scores",
             est = bench_pre["est"], se = bench_pre["se"],
             lo = bench_pre["lo"], hi = bench_pre["hi"],
             n_treated = nrow(nswT), n_control = nrow(nswN),
             perm_p = NA, stringsAsFactors = FALSE)
)
rownames(bench_tbl) <- NULL
write.csv(bench_tbl, file.path(res_dir, "benchmark.csv"), row.names = FALSE)

# ---- A5. Balance table, NSW treated vs NSW control (pre-treatment only) ----
sd_nsw <- sapply(COVS8, function(v) pooled_sd(nswT[[v]], nswN[[v]]))
bal_nsw <- data.frame(
  covariate = COVS8,
  mean_treated = sapply(COVS8, function(v) mean(nswT[[v]])),
  mean_control = sapply(COVS8, function(v) mean(nswN[[v]])),
  smd_treated_minus_control = sapply(COVS8, function(v) smd(nswT[[v]], nswN[[v]], sd_nsw[[v]])),
  stringsAsFactors = FALSE
)
rownames(bal_nsw) <- NULL
write.csv(bal_nsw, file.path(res_dir, "balance_nsw.csv"), row.names = FALSE)
print(bal_nsw, digits = 4)

# ---- A6. Omnibus randomization test (stop rule, DESIGN 3a) -----------------
# Likelihood-ratio statistic for a logit of treatment on the eight covariates,
# compared with 10,000 random re-assignments of treatment.
X <- cbind(1, as.matrix(nsw[, COVS8]))
lr_stat <- function(y) {
  f1 <- glm.fit(x = X, y = y, family = binomial())
  f0 <- glm.fit(x = X[, 1, drop = FALSE], y = y, family = binomial())
  f0$deviance - f1$deviance
}
lr_obs <- lr_stat(nsw$treat)
set.seed(SEED)
lr_perm <- replicate(B_PERM, lr_stat(sample(nsw$treat)))
lr_p <- (1 + sum(lr_perm >= lr_obs)) / (1 + B_PERM)
omni <- data.frame(statistic = "LR for treat ~ 8 covariates (NSW)",
                   observed_LR = lr_obs, df = 8,
                   perm_p = lr_p, n_perm = B_PERM, seed = SEED,
                   stop_threshold = 0.01,
                   stop_triggered = lr_p < 0.01)
write.csv(omni, file.path(res_dir, "omnibus_randomization_test.csv"), row.names = FALSE)

cat("\nBenchmark (re78, Neyman):\n")
print(bench_tbl[, c("outcome", "est", "se", "lo", "hi", "perm_p")], digits = 6)
cat("\nOmnibus LR =", round(lr_obs, 3), " perm p =", round(lr_p, 4), "\n")
cat("Permutation p for benchmark =", round(perm_p, 4), "\n")
