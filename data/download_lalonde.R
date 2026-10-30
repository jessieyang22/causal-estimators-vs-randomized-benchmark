# data/download_lalonde.R
#
# Downloads the LaLonde / Dehejia-Wahba National Supported Work (NSW) files and
# the PSID and CPS comparison groups from Rajeev Dehejia's data page, which is
# hosted by NBER: https://users.nber.org/~rdehejia/nswdata2.html
#
# Raw files are saved to data/raw/ (listed in .gitignore, never committed).
# A file that is already present with the expected row count is not downloaded
# again. After download, the script counts rows, computes MD5 checksums, and
# writes data/raw/manifest.csv so the exact files used can be identified later.
#
# Run from the repository root:
#   Rscript data/download_lalonde.R

if (!file.exists("DESIGN.md")) {
  stop("Run this script from the repository root (the folder that contains DESIGN.md).")
}

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

base_url <- "https://users.nber.org/~rdehejia/data/"

# Row counts stated on the source page (and checked on 2026-10-29).
files <- data.frame(
  name = c("nsw_dw.dta",
           "nswre74_treated.txt",
           "nswre74_control.txt",
           "psid_controls.txt",
           "cps_controls.txt"),
  expected_rows = c(445, 185, 260, 2490, 15992),
  stringsAsFactors = FALSE
)

count_text_rows <- function(path) {
  # Text files have no header line, one person per line.
  length(readLines(path, warn = FALSE))
}

count_dta_rows <- function(path) {
  suppressMessages(library(foreign))
  nrow(read.dta(path))
}

counts <- integer(nrow(files))
for (i in seq_len(nrow(files))) {
  nm <- files$name[i]
  dest <- file.path(raw_dir, nm)
  is_dta <- grepl("\\.dta$", nm)
  have_good_copy <- FALSE
  if (file.exists(dest)) {
    n_now <- if (is_dta) count_dta_rows(dest) else count_text_rows(dest)
    have_good_copy <- identical(as.integer(n_now), as.integer(files$expected_rows[i]))
  }
  if (!have_good_copy) {
    message("Downloading ", nm, " ...")
    download.file(paste0(base_url, nm), destfile = dest, mode = "wb", quiet = TRUE)
  } else {
    message("Already present, skipping download: ", nm)
  }
  counts[i] <- if (is_dta) count_dta_rows(dest) else count_text_rows(dest)
  if (counts[i] != files$expected_rows[i]) {
    stop(sprintf("Row count mismatch for %s: found %d, expected %d",
                 nm, counts[i], files$expected_rows[i]))
  }
}

manifest <- data.frame(
  file = files$name,
  rows = counts,
  bytes = file.size(file.path(raw_dir, files$name)),
  md5 = as.character(tools::md5sum(file.path(raw_dir, files$name))),
  source_url = paste0(base_url, files$name),
  downloaded_or_checked = format(Sys.Date(), "%Y-%m-%d"),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(raw_dir, "manifest.csv"), row.names = FALSE)
print(manifest)
message("All NSW and comparison-group files are present and row counts match the source page.")
