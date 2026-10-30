# src/00_prep_data.R
#
# Reads the raw LaLonde / Dehejia-Wahba files from data/raw/, checks them
# against each other, and writes one cleaned file per sample to data/processed/
# (gitignored). It also writes results/data_summary.csv, which records the row
# counts, value ranges, and missing-value counts used in DATA.md.
#
# This script does not estimate any treatment effect and does not use the
# outcome to build any sample.
#
# Run from the repository root, after data/download_lalonde.R:
#   Rscript src/00_prep_data.R

if (!file.exists("DESIGN.md")) {
  stop("Run this script from the repository root.")
}
suppressMessages(library(foreign))

raw_dir  <- file.path("data", "raw")
proc_dir <- file.path("data", "processed")
res_dir  <- "results"
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Column order of the text files, as stated on the source page.
txt_cols <- c("treat", "age", "education", "black", "hispanic", "married",
              "nodegree", "re74", "re75", "re78")

read_txt <- function(name) {
  df <- read.table(file.path(raw_dir, name), header = FALSE, col.names = txt_cols)
  df
}

# Experimental sample (Stata file). It has 445 rows: 185 treated and 260 control.
nsw <- read.dta(file.path(raw_dir, "nsw_dw.dta"))
nsw <- nsw[, txt_cols]
nsw$sample <- "NSW experimental"

# Text copies of the same experimental sample. They should match the Stata file.
nsw_tr_txt <- read_txt("nswre74_treated.txt")
nsw_co_txt <- read_txt("nswre74_control.txt")
nsw_txt    <- rbind(nsw_tr_txt, nsw_co_txt)

# Consistency check: the Stata file and the two text files should agree.
# Text files store values with about 8 significant digits, so the comparison
# uses a small tolerance.
stopifnot(nrow(nsw_txt) == nrow(nsw))
max_abs_diff <- max(abs(as.matrix(nsw[, txt_cols]) - as.matrix(nsw_txt[, txt_cols])))
rows_match_within_tol <- max_abs_diff < 1e-2
message(sprintf("Max absolute difference, Stata vs text copies of NSW: %.6f", max_abs_diff))
if (!rows_match_within_tol) {
  stop("The Stata file and the text files for the NSW sample disagree.")
}

# Comparison groups. They have no treatment (all treat = 0).
psid <- read_txt("psid_controls.txt")
cps  <- read_txt("cps_controls.txt")
psid$sample <- "PSID comparison"
cps$sample  <- "CPS comparison"
stopifnot(all(psid$treat == 0), all(cps$treat == 0))

# Treatment indicator checks.
tab_nsw <- table(nsw$treat)
stopifnot(tab_nsw[["1"]] == 185, tab_nsw[["0"]] == 260)

# Save the cleaned files. Each is a plain data frame with the same ten
# variables plus a 'sample' label. Missing values are counted and reported,
# but none are imputed or dropped here.
saveRDS(nsw,  file.path(proc_dir, "nsw_experimental.rds"))
saveRDS(psid, file.path(proc_dir, "psid_comparison.rds"))
saveRDS(cps,  file.path(proc_dir, "cps_comparison.rds"))

# Summary table: counts, ranges, missingness, and a top-coding check.
# No outcome means are computed, and the outcome (re78) is not counted by arm.
summ_one <- function(df, label) {
  out <- data.frame(
    sample = label,
    variable = txt_cols,
    n_rows = nrow(df),
    n_missing = sapply(df[, txt_cols], function(x) sum(is.na(x))),
    min = sapply(df[, txt_cols], function(x) min(x, na.rm = TRUE)),
    max = sapply(df[, txt_cols], function(x) max(x, na.rm = TRUE)),
    # Rows sitting exactly at the column maximum, for the earnings variables
    # only. A large count would suggest top-coding (a cap on recorded earnings).
    n_at_max = sapply(txt_cols, function(v) {
      if (v %in% c("re74", "re75", "re78")) sum(df[[v]] == max(df[[v]], na.rm = TRUE)) else NA_integer_
    }),
    # Zeros in the pre-treatment earnings variables only. Zeros in re78 are
    # deliberately not counted here, so no outcome counts are produced.
    n_zero_pretreatment = sapply(txt_cols, function(v) {
      if (v %in% c("re74", "re75")) sum(df[[v]] == 0) else NA_integer_
    }),
    n_treated = sum(df$treat == 1),
    n_control = sum(df$treat == 0),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}
summary_tbl <- rbind(
  summ_one(nsw,  "NSW experimental"),
  summ_one(psid, "PSID comparison"),
  summ_one(cps,  "CPS comparison")
)
write.csv(summary_tbl, file.path(res_dir, "data_summary.csv"), row.names = FALSE)

# Sample sizes and missingness, printed for the log.
print(summary_tbl[summary_tbl$variable == "treat", c("sample", "n_rows", "n_treated", "n_control")], row.names = FALSE)
message("Total missing values across all three samples: ",
        sum(summary_tbl$n_missing))
message("Wrote processed files to ", proc_dir, " and summary to ", file.path(res_dir, "data_summary.csv"))

# Save session info (DESIGN.md Section 6 requires it for reproducibility).
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session_info_prep.txt"))
