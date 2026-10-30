# src/08_figures.R
#
# Figures from results/ (DESIGN.md Component F and Section 6). Also writes
# results/session_info.txt for reproducibility (DESIGN.md Section 6).
#   results/robustness_plot.png : ATT estimate with 95% interval for each
#                                  estimator and grid variant, against the NSW
#                                  benchmark (shaded band = benchmark interval).
# Run from the repository root: Rscript src/08_figures.R

if (!file.exists("DESIGN.md")) stop("Run this script from the repository root.")
suppressMessages(library(ggplot2))
res_dir <- "results"

grid <- read.csv(file.path(res_dir, "robustness_grid.csv"), stringsAsFactors = FALSE)
prim <- read.csv(file.path(res_dir, "observational_vs_benchmark.csv"), stringsAsFactors = FALSE)
prim$grid_type <- "primary"
both <- rbind(prim[, names(grid)], grid)

both$spec <- factor(paste0(both$method, ifelse(both$grid_type == "primary", "", paste0(" [", both$grid_type, "]"))),
                    levels = rev(unique(paste0(both$method, ifelse(both$grid_type == "primary", "", paste0(" [", both$grid_type, "]"))))))
both$spec_lab <- both$spec
both$method_lab <- factor(both$method,
                          levels = c("naive", "reg", "match", "ipw", "aipw", "did"),
                          labels = c("1 Naive", "2 Regression", "3 Matching", "4 IPW", "5 AIPW", "6 DiD"))
bench <- prim$bench[1]
bench_lo <- prim$bench_lo[1]
bench_hi <- prim$bench_hi[1]
bench_df <- data.frame(group = c("PSID", "CPS"), lo = bench_lo, hi = bench_hi, est = bench)

both$variant <- ifelse(both$grid_type == "primary", "primary", both$grid_type)

p <- ggplot(both, aes(y = spec, x = att, colour = variant)) +
  geom_rect(data = bench_df, inherit.aes = FALSE,
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.5) +
  geom_vline(data = bench_df, aes(xintercept = est), linetype = 2) +
  geom_errorbar(aes(xmin = att_lo, xmax = att_hi), orientation = "y", width = 0.2) +
  geom_point(size = 1.8) +
  facet_wrap(~ group, ncol = 2, scales = "free_x") +
  labs(title = "Observational ATT estimates against the NSW experimental benchmark",
       subtitle = "Points: estimate. Bars: 95% interval (method's own SE). Dashed line and grey band: NSW benchmark and its 95% interval.",
       x = "Estimated effect on 1978 earnings (dollars)", y = NULL, colour = "Specification") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))

ggsave(file.path(res_dir, "robustness_plot.png"), p, width = 10, height = 7, dpi = 150)

# Session information for reproducibility.
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session_info.txt"))
message("Wrote results/robustness_plot.png and results/session_info.txt")
