# src/helpers.R
#
# Shared code for Components A to F. Logged as Deviation 4(a) in DESIGN.md.
# Nothing in this file reads or uses a result. It defines:
#   - how the data are loaded,
#   - the observational estimators (naive, regression adjustment, matching,
#     IPW, AIPW, and difference-in-differences),
#   - the bootstrap that gives the bias confidence interval and the precision
#     guard, and
#   - the balance (standardized mean difference) helpers.
#
# Plain-English summary of the main idea:
#   The benchmark is the randomized NSW difference in 1978 earnings. Each
#   observational method gives its own estimate of the same effect using the
#   NSW treated people plus a comparison group (PSID or CPS). The "bias" is the
#   observational estimate minus the benchmark. A bootstrap (re-sampling the
#   people with replacement and repeating the whole calculation 2,000 times)
#   gives the uncertainty of that bias.

suppressMessages({
  library(Matching)
  library(parallel)
})

COVS8 <- c("age", "education", "black", "hispanic", "married",
           "nodegree", "re74", "re75")
COVS_PLACEBO <- setdiff(COVS8, "re75")   # Section 3e: re75 is the placebo outcome
COVS_NO74    <- setdiff(COVS8, "re74")   # Section 3f: covariate-set variation
SEED <- 20261029
B_BOOT <- 2000

# Outcome definitions. Each takes a data frame with the NSW or comparison rows.
OUT_RE78 <- function(d) d$re78
OUT_RE75 <- function(d) d$re75
OUT_CHG78_75 <- function(d) d$re78 - d$re75   # DiD primary: change 1975 to 1978
OUT_CHG75_74 <- function(d) d$re75 - d$re74   # DiD pre-trend: change 1974 to 1975

# Load the cleaned samples written by src/00_prep_data.R.
load_samples <- function(root = ".") {
  proc <- file.path(root, "data", "processed")
  nsw <- readRDS(file.path(proc, "nsw_experimental.rds"))
  list(
    nswT = nsw[nsw$treat == 1, ],
    nswN = nsw[nsw$treat == 0, ],
    psid = readRDS(file.path(proc, "psid_comparison.rds")),
    cps  = readRDS(file.path(proc, "cps_comparison.rds"))
  )
}

# Difference in means with the Neyman (independent groups) standard error and a
# normal 95% interval. Used for the benchmark and for the naive estimator.
neyman <- function(yT, yC) {
  est <- mean(yT) - mean(yC)
  se <- sqrt(var(yT) / length(yT) + var(yC) / length(yC))
  c(est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se)
}

# Pooled standard deviation used as the denominator of every SMD.
pooled_sd <- function(xT, xC) sqrt((var(xT) + var(xC)) / 2)

# Standardized mean difference (treated minus comparison). If weights are given,
# the comparison mean is weighted. The denominator is always the unadjusted
# pooled SD passed in (Deviation 4(j)).
smd <- function(xT, xC, sd_pool, wC = NULL) {
  mC <- if (is.null(wC)) mean(xC) else sum(wC * xC) / sum(wC)
  (mean(xT) - mC) / sd_pool
}

# Logit propensity score for "treated (NSW) vs comparison", main effects only.
fit_ps <- function(dT, dC, covs) {
  d <- rbind(dT[, covs, drop = FALSE], dC[, covs, drop = FALSE])
  d$treat <- c(rep(1, nrow(dT)), rep(0, nrow(dC)))
  fml <- as.formula(paste("treat ~", paste(covs, collapse = " + ")))
  fit <- glm(fml, data = d, family = binomial())
  p <- as.numeric(fitted(fit))
  list(pT = p[seq_len(nrow(dT))],
       pC = p[nrow(dT) + seq_len(nrow(dC))])
}

# One observational estimate of the ATT.
#   method : "naive", "did", "reg", "match", "ipw", "aipw"
#   yT, yC : outcome vectors for the NSW treated (T) and the comparison group (C)
#   dT, dC : covariate data frames for the same rows
#   covs   : covariate names used for adjustment
#   trim   : propensity cap for IPW and AIPW (Deviation 4(h))
#   caliper: matching caliper in SD of the logit propensity score; NULL = none
# Returns a list with att, se (the primary SE for the method, NA if none),
# dropped (treated units dropped by a caliper), and, for IPW, AIPW and
# matching, the comparison weights, propensity scores and matched-control
# counts used for diagnostics.
estimate <- function(method, yT, yC, dT, dC, covs,
                     trim = 0.01, caliper = 0.2) {
  nT <- length(yT)
  nC <- length(yC)
  out <- list(att = NA_real_, se = NA_real_, dropped = 0L,
              w = NULL, pT = NULL, pC = NULL, ctrl_count = NULL)

  if (method %in% c("naive", "did")) {
    # "did" is the same difference in means, applied to change scores by the caller.
    ne <- neyman(yT, yC)
    out$att <- unname(ne["est"])
    out$se <- unname(ne["se"])
    return(out)
  }

  if (method == "reg") {
    dfC <- data.frame(y = yC, dC[, covs, drop = FALSE])
    fit <- lm(y ~ ., data = dfC)
    pr <- predict(fit, newdata = dT[, covs, drop = FALSE])
    out$att <- mean(yT) - mean(pr)
    return(out)
  }

  ps <- fit_ps(dT, dC, covs)
  out$pT <- ps$pT
  out$pC <- ps$pC

  if (method == "match") {
    lp <- qlogis(c(ps$pT, ps$pC))
    Y <- c(yT, yC)
    Tr <- c(rep(1, nT), rep(0, nC))
    m <- Match(Y = Y, Tr = Tr, X = lp, M = 1, caliper = caliper,
               replace = TRUE, estimand = "ATT")
    out$att <- unname(m$est)
    out$se <- unname(m$se)
    out$dropped <- nT - length(unique(m$index.treated))
    # Match() gives each matched row a weight (1 / number of tied controls for
    # that treated unit). Summing the row weights per comparison unit gives the
    # weight each control carries in the matched sample; the total is nT.
    stopifnot(length(m$weights) == length(m$index.control))
    cc <- numeric(nC)
    agg <- rowsum(m$weights, m$index.control - nT)
    cc[as.integer(rownames(agg))] <- agg[, 1]
    out$ctrl_count <- cc
    return(out)
  }

  pCc <- pmin(pmax(ps$pC, trim), 1 - trim)
  W <- pCc / (1 - pCc)           # ATT odds weights for the comparison group
  w <- W / sum(W)                # Hajek normalization (sums to one)

  if (method == "ipw") {
    out$att <- mean(yT) - sum(w * yC)
    out$w <- w
    return(out)
  }

  if (method == "aipw") {
    dfC <- data.frame(y = yC, dC[, covs, drop = FALSE])
    fit <- lm(y ~ ., data = dfC)
    m0C <- as.numeric(fitted(fit))
    m0T <- as.numeric(predict(fit, newdata = dT[, covs, drop = FALSE]))
    R <- yC - m0C
    r <- sum(w * R)
    out$att <- mean(yT - m0T) - r
    # Influence-function SE (ratio linearization). Treats the fitted propensity
    # and outcome models as known; the bootstrap is the cross-check.
    u <- W * (R - r) / mean(W)
    out$se <- sqrt(var(yT - m0T) / nT + var(u) / nC)
    out$w <- w
    return(out)
  }
  stop("unknown method: ", method)
}

# Effective sample size (Kish) of normalized weights.
kish_ess <- function(w) 1 / sum(w^2)

# One observational cell (method + outcome + covariates) and its bias bootstrap.
#   cell : list(id, group, estimator, method, obs, bench, covs, trim, caliper, covset, grid_type)
#          obs and bench are outcome functions (see OUT_* above).
#   S    : list(nswT, nswN) (NSW treated and control)
#   comp : comparison group data frame
# Bias = observational ATT minus benchmark ATT (Deviation 4(b)).
# Bootstrap: each replicate re-samples NSW treated, NSW control and the
# comparison group with replacement, then re-runs the whole estimator and the
# benchmark. Replicate b uses seed SEED + b, so results do not depend on the
# number of cores.
run_cell <- function(cell, S, comp, B = B_BOOT, cores = 1L) {
  dT <- S$nswT
  dN <- S$nswN
  yT <- cell$obs(dT)
  yC <- cell$obs(comp)
  bT <- cell$bench(dT)
  bN <- cell$bench(dN)

  pt <- estimate(cell$method, yT, yC, dT, dC = comp, covs = cell$covs,
                 trim = cell$trim, caliper = cell$caliper)
  bench <- neyman(bT, bN)
  bias_pt <- pt$att - unname(bench["est"])

  one_rep <- function(b) {
    set.seed(SEED + b)
    iT <- sample.int(nrow(dT), replace = TRUE)
    iN <- sample.int(nrow(dN), replace = TRUE)
    iC <- sample.int(nrow(comp), replace = TRUE)
    e <- suppressWarnings(estimate(cell$method, yT[iT], yC[iC],
                                   dT[iT, , drop = FALSE],
                                   comp[iC, , drop = FALSE],
                                   covs = cell$covs, trim = cell$trim,
                                   caliper = cell$caliper))
    bb <- mean(bT[iT]) - mean(bN[iN])
    c(att = e$att, bench = bb, bias = e$att - bb, dropped = e$dropped)
  }
  reps <- do.call(rbind, mclapply(seq_len(B), one_rep, mc.cores = cores))
  ok <- is.finite(reps[, "bias"])

  bias_ci <- unname(quantile(reps[ok, "bias"], c(0.025, 0.975)))
  H_o <- (bias_ci[2] - bias_ci[1]) / 2
  H_b <- 1.96 * unname(bench["se"])
  informative <- H_o <= 2 * H_b
  covers0 <- bias_ci[1] <= 0 && bias_ci[2] >= 0
  pass <- informative && covers0

  # ATT interval: method's own SE, or the bootstrap SE for DiD, regression, IPW.
  se_method <- switch(cell$method,
    naive = "neyman", did = "bootstrap", reg = "bootstrap", ipw = "bootstrap",
    match = "abadie_imbens", aipw = "influence_function")
  att_se <- if (se_method == "bootstrap") sd(reps[ok, "att"]) else pt$se
  att_boot_se <- sd(reps[ok, "att"])

  list(
    row = data.frame(
      cell_id = cell$id, group = cell$group, estimator = cell$estimator,
      method = cell$method, grid_type = cell$grid_type,
      covs = paste(cell$covs, collapse = ";"),
      trim = cell$trim, caliper = ifelse(is.null(cell$caliper), NA, cell$caliper),
      n_treated = nrow(dT), n_comp = nrow(comp),
      att = pt$att, att_se_method = se_method, att_se = att_se,
      att_boot_se = att_boot_se,
      att_lo = pt$att - 1.96 * att_se, att_hi = pt$att + 1.96 * att_se,
      bench = unname(bench["est"]), bench_se = unname(bench["se"]),
      bench_lo = unname(bench["lo"]), bench_hi = unname(bench["hi"]),
      H_b = H_b,
      bias = bias_pt, bias_lo = bias_ci[1], bias_hi = bias_ci[2],
      bias_pct = 100 * bias_pt / unname(bench["est"]),
      H_o = H_o, informative = informative, bias_ci_covers_zero = covers0,
      pass = pass,
      mean_dropped_treated = mean(reps[ok, "dropped"]),
      point_dropped_treated = pt$dropped,
      n_boot_ok = sum(ok),
      stringsAsFactors = FALSE
    ),
    reps = reps,
    point = pt,
    bench = bench
  )
}

# Standard cell definitions for one comparison group. "group" is the label.
primary_cells <- function(group) {
  list(
    list(id = paste0(group, "_naive"), estimator = "1 Naive difference",
         method = "naive", obs = OUT_RE78, bench = OUT_RE78, covs = COVS8,
         trim = 0.01, caliper = 0.2, grid_type = "primary"),
    list(id = paste0(group, "_reg"), estimator = "2 Regression adjustment",
         method = "reg", obs = OUT_RE78, bench = OUT_RE78, covs = COVS8,
         trim = 0.01, caliper = 0.2, grid_type = "primary"),
    list(id = paste0(group, "_match"), estimator = "3 Nearest-neighbour matching",
         method = "match", obs = OUT_RE78, bench = OUT_RE78, covs = COVS8,
         trim = 0.01, caliper = 0.2, grid_type = "primary"),
    list(id = paste0(group, "_ipw"), estimator = "4 Inverse probability weighting",
         method = "ipw", obs = OUT_RE78, bench = OUT_RE78, covs = COVS8,
         trim = 0.01, caliper = 0.2, grid_type = "primary"),
    list(id = paste0(group, "_aipw"), estimator = "5 Doubly robust (AIPW)",
         method = "aipw", obs = OUT_RE78, bench = OUT_RE78, covs = COVS8,
         trim = 0.01, caliper = 0.2, grid_type = "primary"),
    list(id = paste0(group, "_did"), estimator = "6 Difference-in-differences",
         method = "did", obs = OUT_CHG78_75, bench = OUT_RE78, covs = COVS8,
         trim = 0.01, caliper = 0.2, grid_type = "primary")
  )
}

# Adds the group label and the estimator label to a cell list.
tag <- function(cells, group) lapply(cells, function(x) { x$group <- group; x })

# Balance helper for one comparison group, using the unadjusted pooled SD of
# NSW treated and the full comparison group (Deviation 4(j)).
sd_full <- function(S, comp, covs) {
  sapply(covs, function(v) pooled_sd(S$nswT[[v]], comp[[v]]))
}

# Core count for mclapply. Results do not depend on it (each replicate has its
# own seed). Windows has no fork, so mclapply must run on one core there.
default_cores <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  min(8L, parallel::detectCores())
}

# Placebo cells (DESIGN.md 3e, Deviation 4(c)): re75 outcome, re75 removed from
# the covariates. Used by 04_placebo.R and the seed run (09).
placebo_cells <- function(label) {
  list(
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
}

# Pre-trend cell for DiD (DESIGN.md 3e): change 1974 to 1975, unadjusted.
pretrend_cell <- function(label) {
  list(id = paste0(label, "_pretrend_did"), estimator = "6 Difference-in-differences",
       method = "did", obs = OUT_CHG75_74, bench = OUT_CHG75_74,
       covs = COVS8, trim = 0.01, caliper = 0.2, grid_type = "pretrend", group = label)
}

# Robustness grid cells (DESIGN.md 3f, Deviation 4(m)). Used by 06 and 09.
grid_cells <- function(label) {
  base <- function(id, est, method, covs, trim = 0.01, caliper = 0.2, variant) {
    list(id = paste0(label, "_", id), estimator = est, method = method,
         obs = OUT_RE78, bench = OUT_RE78, covs = covs, trim = trim,
         caliper = caliper, grid_type = variant)
  }
  list(
    base("ipw_trim05", "4 Inverse probability weighting", "ipw", COVS8, trim = 0.05, variant = "trim 0.05"),
    base("aipw_trim05", "5 Doubly robust (AIPW)", "aipw", COVS8, trim = 0.05, variant = "trim 0.05"),
    base("match_cal01", "3 Nearest-neighbour matching", "match", COVS8, caliper = 0.1, variant = "caliper 0.1 SD"),
    base("match_nocal", "3 Nearest-neighbour matching", "match", COVS8, caliper = NULL, variant = "no caliper"),
    base("reg_no74", "2 Regression adjustment", "reg", COVS_NO74, variant = "without re74"),
    base("match_no74", "3 Nearest-neighbour matching", "match", COVS_NO74, variant = "without re74"),
    base("ipw_no74", "4 Inverse probability weighting", "ipw", COVS_NO74, variant = "without re74"),
    base("aipw_no74", "5 Doubly robust (AIPW)", "aipw", COVS_NO74, variant = "without re74")
  )
}
