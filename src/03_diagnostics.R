# src/03_diagnostics.R
#
# Component C: diagnostics for why an estimator fails (DESIGN.md 4, Component C).
#   - covariate balance before and after IPW, AIPW and matching (SMD, with the
#     unweighted pooled-SD denominator, Deviation 4(j));
#   - SMD of the logit propensity score before and after;
#   - Kish effective sample size and largest-weight ratio (IPW and AIPW);
#   - matching counts (matched, dropped, distinct controls);
#   - overlap plots of the propensity score, treated versus comparison;
#   - how much of the IPW and AIPW weight comes from propensity scores that are
#     capped (raw score below 0.01, or below 0.05 for the grid variant), and the
#     glm warnings about fitted probabilities of 0 or 1 (review item, Deviation 7).
#     Written to results/overlap_cap_shares.csv.
# Run from the repository root: Rscript src/03_diagnostics.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
suppressMessages(library(ggplot2))
res_dir <- "results"
S <- load_samples()
nswT <- S$nswT

bal_rows <- list()
diag_rows <- list()

for (g in c("psid", "cps")) {
  label <- if (g == "psid") "PSID" else "CPS"
  comp <- S[[g]]
  sdv <- sd_full(S, comp, COVS8)
  yT <- nswT$re78
  yC <- comp$re78

  # Propensity score scale for the overlap check (logit), pooled SD unadjusted.
  glm_warn <- character(0)
  pt_ps <- withCallingHandlers(fit_ps(nswT, comp, COVS8),
    warning = function(w) {
      glm_warn <<- c(glm_warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  lT <- qlogis(pt_ps$pT)
  lC <- qlogis(pt_ps$pC)
  sd_lp <- pooled_sd(lT, lC)

  # Overlap plot (primary propensity model, logit scale).
  plot_df <- rbind(
    data.frame(logit_ps = lT, group = "NSW treated"),
    data.frame(logit_ps = lC, group = paste(label, "comparison"))
  )
  p <- ggplot(plot_df, aes(x = logit_ps, colour = group, fill = group)) +
    geom_density(alpha = 0.25, linewidth = 0.6) +
    labs(title = paste("Overlap of propensity scores:", label),
         subtitle = "Logit of the primary propensity score (8 covariates, main effects)",
         x = "Logit of propensity score", y = "Density", colour = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(res_dir, paste0("overlap_", g, ".png")), p,
         width = 7, height = 4.5, dpi = 150)

  # Balance before adjustment (unweighted comparison group).
  before <- sapply(COVS8, function(v) smd(nswT[[v]], comp[[v]], sdv[[v]]))
  before_ps <- smd(lT, lC, sd_lp)

  for (m in c("reg", "match", "ipw", "aipw")) {
    e <- estimate(m, yT, yC, nswT, comp, COVS8, trim = 0.01, caliper = 0.2)
    if (m == "reg") {
      bal_rows[[paste(g, m)]] <- data.frame(
        group = label, method = m, covariate = c(COVS8, "logit_ps"),
        smd_before = c(before, before_ps), smd_after = NA_real_,
        note = "Not applicable: regression adjustment has no weights",
        stringsAsFactors = FALSE)
      next
    }
    if (m %in% c("ipw", "aipw")) {
      wC <- e$w
    } else {
      wC <- e$ctrl_count / sum(e$ctrl_count)
    }
    after <- sapply(COVS8, function(v) smd(nswT[[v]], comp[[v]], sdv[[v]], wC = wC))
    after_ps <- smd(lT, lC, sd_lp, wC = wC)
    bal_rows[[paste(g, m)]] <- data.frame(
      group = label, method = m, covariate = c(COVS8, "logit_ps"),
      smd_before = c(before, before_ps), smd_after = c(after, after_ps),
      note = if (m == "match") "Weights = each control's total Match() weight (tied controls share one unit; the weights sum to the number of treated units)" else "",
      stringsAsFactors = FALSE)

    max_i <- which.max(abs(c(after, after_ps)))
    cov_names <- c(COVS8, "logit_ps")
    row <- data.frame(group = label, method = m,
                      n_treated = length(lT), n_comp = length(lC),
                      matched_treated = NA_integer_, dropped_treated = NA_integer_,
                      distinct_controls = NA_integer_,
                      ess = NA_real_, ess_fraction = NA_real_,
                      max_weight_ratio = NA_real_, overlap_failure = NA,
                      stringsAsFactors = FALSE)
    if (m %in% c("ipw", "aipw")) {
      row$ess <- kish_ess(wC)
      row$ess_fraction <- row$ess / length(wC)
      row$max_weight_ratio <- max(wC) / mean(wC)
      row$overlap_failure <- row$ess < 50 || row$max_weight_ratio > 10
    } else {
      row$matched_treated <- length(lT) - e$dropped
      row$dropped_treated <- e$dropped
      row$distinct_controls <- sum(e$ctrl_count > 0)
    }
    row$largest_abs_smd_after <- max(abs(c(after, after_ps)))
    row$largest_smd_covariate <- cov_names[max_i]
    row$balance_failure <- max(abs(c(after, after_ps))) > 0.1
    row$primary_att <- e$att
    diag_rows[[paste(g, m)]] <- row
  }
}

# ---- Cap shares for IPW and AIPW (review item, Deviation 7) ------------------
# A comparison unit with raw propensity score below the cap gets the cap value,
# so all such units share one weight. Shares below are of the 15,992 or 2,490
# comparison units, and of the total Hajek weight.
cap_rows <- list()
for (g in c("psid", "cps")) {
  label <- if (g == "psid") "PSID" else "CPS"
  comp <- S[[g]]
  cap_warn <- character(0)
  ps <- withCallingHandlers(fit_ps(nswT, comp, COVS8),
    warning = function(w) {
      cap_warn <<- c(cap_warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  pC <- ps$pC
  nC <- length(pC)
  for (cap in c(0.01, 0.05)) {
    pCc <- pmin(pmax(pC, cap), 1 - cap)
    W <- pCc / (1 - pCc)
    w <- W / sum(W)
    capped <- pC < cap
    cap_rows[[paste(g, cap)]] <- data.frame(
      group = label, cap = cap,
      n_comp = nC,
      n_capped_low = sum(pC < cap),
      share_capped_low = mean(pC < cap),
      n_capped_high = sum(pC > 1 - cap),
      share_weight_from_capped_low = sum(w[capped]) / sum(w),
      max_weight_ratio = max(w) / mean(w),
      ess = kish_ess(w),
      n_raw_ps_exact_0_or_1 = sum(pC == 0 | pC == 1),
      glm_warning = if (length(cap_warn) > 0) paste(unique(cap_warn), collapse = " | ") else "",
      stringsAsFactors = FALSE)
  }
}
cap_tbl <- do.call(rbind, cap_rows)
rownames(cap_tbl) <- NULL
write.csv(cap_tbl, file.path(res_dir, "overlap_cap_shares.csv"), row.names = FALSE)
print(cap_tbl, digits = 4)

bal <- do.call(rbind, bal_rows)
rownames(bal) <- NULL
write.csv(bal, file.path(res_dir, "balance_after_adjustment.csv"), row.names = FALSE)

diag <- do.call(rbind, diag_rows)
rownames(diag) <- NULL
write.csv(diag, file.path(res_dir, "diagnostics_weights_matching.csv"), row.names = FALSE)

print(diag, digits = 4)
print(bal[which(bal$covariate == "logit_ps" | abs(bal$smd_after) > 0.1), ], digits = 3)
