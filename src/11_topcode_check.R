# src/11_topcode_check.R
#
# Point-estimate sensitivity to the CPS top-coding pattern (DATA.md Section 5),
# added 2026-10-30 (Deviation 7). POINT ESTIMATES ONLY, no bootstrap: this shows
# whether the six estimates move when the rows at the file maximum are removed.
# It is not a protocol test and no pass/fail verdict is taken from it.
#   full      : all 15,992 CPS rows (primary)
#   drop_any  : drop comparison rows at the file maximum of re74, re75 or re78
#   drop_re78 : drop comparison rows at the file maximum of re78 only
# Output: results/cps_topcode_point_check.csv
# Run from the repository root: Rscript src/11_topcode_check.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
source(file.path("src", "helpers.R"))
S <- load_samples()
nswT <- S$nswT
nswN <- S$nswN
cps <- S$cps
mx <- sapply(c("re74", "re75", "re78"), function(v) max(cps[[v]]))
at_max_any <- cps$re74 == mx["re74"] | cps$re75 == mx["re75"] | cps$re78 == mx["re78"]
at_max_78 <- cps$re78 == mx["re78"]
samples <- list(full = cps, drop_any = cps[!at_max_any, ], drop_re78 = cps[!at_max_78, ])

bench <- neyman(nswT$re78, nswN$re78)
rows <- list()
for (nm in names(samples)) {
  comp <- samples[[nm]]
  for (m in c("naive", "reg", "match", "ipw", "aipw", "did")) {
    yT <- if (m == "did") nswT$re78 - nswT$re75 else nswT$re78
    yC <- if (m == "did") comp$re78 - comp$re75 else comp$re78
    e <- estimate(m, yT, yC, nswT, comp, COVS8, trim = 0.01, caliper = 0.2)
    rows[[paste(nm, m)]] <- data.frame(
      sample = nm, n_comp = nrow(comp), method = m, att = e$att,
      bias_point = e$att - unname(bench["est"]), stringsAsFactors = FALSE)
  }
}
tbl <- do.call(rbind, rows); rownames(tbl) <- NULL
tbl$share_comp_dropped <- 1 - tbl$n_comp / nrow(cps)
write.csv(tbl, file.path("results", "cps_topcode_point_check.csv"), row.names = FALSE)
cat(sprintf("CPS rows at any earnings maximum: %d (%.1f%%); at re78 maximum: %d (%.1f%%)\n",
            sum(at_max_any), 100 * mean(at_max_any), sum(at_max_78), 100 * mean(at_max_78)))
print(tbl[, c("sample", "n_comp", "method", "att", "bias_point")], digits = 6)
