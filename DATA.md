# Data record: causal-estimators-vs-randomized-benchmark

Status: data acquired and checked on 2026-10-29. No treatment effect, p-value on an outcome, or comparison of estimates has been computed.

**Disclosure of an outcome count.** During the data check on 2026-10-29, after the design was locked, one cross-tabulation of the outcome was printed once: `re78 == 0` by treatment arm in the NSW experimental sample gave 45 of 185 treated and 92 of 260 controls with zero 1978 earnings. This is a descriptive tally, not an effect estimate. No protocol choice was made or changed in response to it. It is recorded as Deviation 3 in `DESIGN.md`. No further `re78` statistic has been computed by arm. The ranges and counts in this file, other than the item above, are not split by arm.

This file records what was actually observed in the downloaded files. Every number below comes from `data/download_lalonde.R`, `src/00_prep_data.R`, and the files they write (`data/raw/manifest.csv`, `results/data_summary.csv`).

---

## 1. What was used

| Role | Dataset | Used? |
|---|---|---|
| Benchmark (randomized) and comparison groups | LaLonde / Dehejia-Wahba National Supported Work (NSW) files, with PSID and CPS comparison groups | Yes |
| Second randomized dataset (Hillstrom e-mail test) | MineThatData 2008 e-mail challenge file | No. Dropped because its terms could not be confirmed. See DESIGN.md, Section 9 (Deviations). |

---

## 2. Source and terms

Source page: https://users.nber.org/~rdehejia/nswdata2.html (Rajeev Dehejia's data page, hosted by NBER).

Terms, as printed on the page: "These data are distributed for attributable non-commercial use (CC by NC)." The page asks that these three papers be cited: Dehejia and Wahba (1999, *JASA*), Dehejia and Wahba (2002, *Review of Economics and Statistics*), and LaLonde (1986, *AER*). Full titles should be taken from the journals before the README is published; they were not copied from the source page.

The page does not state the price year for the dollar amounts. Amounts are therefore reported as "dollars" only.

---

## 3. Files downloaded

Downloaded by `Rscript data/download_lalonde.R` into `data/raw/` (gitignored). Checksums are MD5 values from `data/raw/manifest.csv`.

| File | Rows (observed) | Bytes | MD5 | Contents |
|---|---|---|---|---|
| `nsw_dw.dta` | 445 | 28,598 | e3367ddd449e25cb4c54c04263978800 | NSW Dehejia-Wahba experimental sample, Stata format (185 treated, 260 control) |
| `nswre74_treated.txt` | 185 | 29,785 | 48b20ef531d899761fc8bf3235d3d70b | Same treated units, text format |
| `nswre74_control.txt` | 260 | 41,860 | bb289d7436c31eb8ea270f0c41b3aa38 | Same control units, text format |
| `psid_controls.txt` | 2,490 | 400,890 | 228b82bd417c3d8a9252026963929304 | PSID comparison group (all treat = 0) |
| `cps_controls.txt` | 15,992 | 2,574,712 | f90c33aa193e45056f487929bb432a3a | CPS comparison group (all treat = 0) |

Total raw size: about 3.0 MB. None of the raw files is committed to git.

Row counts match the numbers on the source page (185, 260, 2,490, 15,992).

Consistency check: the Stata file and the two text files were compared row by row (same order). The largest absolute difference in any variable was 0.00125. This comes from the text files storing about 8 significant digits, which is a rounding difference and not a data error. The tolerance used was 0.01.

---

## 4. Samples used in the analysis

| Sample | Rows | Treated (treat = 1) | Control or comparison (treat = 0) |
|---|---|---|---|
| NSW experimental (randomized) | 445 | 185 | 260 |
| PSID comparison | 2,490 | 0 | 2,490 |
| CPS comparison | 15,992 | 0 | 15,992 |

Each sample is saved as a cleaned file in `data/processed/` (gitignored): `nsw_experimental.rds`, `psid_comparison.rds`, `cps_comparison.rds`.

---

## 5. Data dictionary

Variable names and order follow the source page. The `.dta` file also has a `data_id` label (the text "Dehejia-Wahba Sample" on all rows), which is not used.

| Variable | Meaning | Type and coding | Observed range, NSW experimental | Observed range, PSID | Observed range, CPS |
|---|---|---|---|---|---|
| `treat` | Treatment indicator: 1 = received NSW job training, 0 = not treated | 0/1 | 0 to 1 | 0 to 0 | 0 to 0 |
| `age` | Age in years | integer | 17 to 55 | 18 to 55 | 16 to 55 |
| `education` | Years of schooling | integer | 3 to 16 | 0 to 17 | 0 to 18 |
| `black` | 1 = Black | 0/1 | 0 to 1 | 0 to 1 | 0 to 1 |
| `hispanic` | 1 = Hispanic | 0/1 | 0 to 1 | 0 to 1 | 0 to 1 |
| `married` | 1 = married | 0/1 | 0 to 1 | 0 to 1 | 0 to 1 |
| `nodegree` | 1 = no high school degree | 0/1 | 0 to 1 | 0 to 1 | 0 to 1 |
| `re74` | Earnings in 1974 (pre-program), dollars | real number | 0 to 39,570.68 | 0 to 137,148.68 | 0 to 25,862.32 |
| `re75` | Earnings in 1975 (used as the placebo outcome), dollars | real number | 0 to 25,142.24 | 0 to 156,653.23 | 0 to 25,243.55 |
| `re78` | Earnings in 1978 (the outcome), dollars | real number | 0 to 60,307.93 | 0 to 121,173.58 | 0 to 25,564.67 |

Note on `education` in the PSID and CPS files: the minimum is 0, which the source page does not explain further. It is kept as stated in the file.

Note on earnings maxima in the CPS file (top-coding check, from `results/data_summary.csv`, column `n_at_max`): the number of rows exactly at each column's maximum is:

| Sample | `re74` max (rows at max) | `re75` max (rows at max) | `re78` max (rows at max) |
|---|---|---|---|
| NSW experimental | 39,570.68 (1) | 25,142.24 (1) | 60,307.93 (1) |
| PSID comparison | 137,148.68 (1) | 156,653.23 (1) | 121,173.58 (1) |
| CPS comparison | 25,862.32 (**3,052**) | 25,243.55 (**3,007**) | 25,564.67 (**4,145**) |

In the CPS file about 19% of rows sit at the `re74` and `re75` maximum, and about 26% at the `re78` maximum, while the NSW and PSID files have one row at each maximum. This pattern looks like a cap (top-coding) on recorded earnings in the CPS file. The source page does not explain it, and this has not been checked against any other documentation. Treat it as an observed pattern, not a confirmed cause.

---

## 6. Dates and time coverage

- The files contain no date variable and no program start or end date.
- The only time information is in the variable names: `re74` (1974), `re75` (1975), `re78` (1978). These labels come from the source page's variable list.
- The source page does not give the program dates. The placebo design (DESIGN.md, Section 3e) treats `re75` as pre-program. That is an assumption, not something verified from the files.
- The Dehejia-Wahba subsample keeps only people with 1974 earnings data. Its time coverage is therefore 1974 to 1978 earnings.

---

## 7. Missingness

Checked with `is.na()` on every variable in every sample (`results/data_summary.csv`).

| Sample | Variables checked | Missing values |
|---|---|---|
| NSW experimental | 10 | 0 |
| PSID comparison | 10 | 0 |
| CPS comparison | 10 | 0 |

No missing values were found, so no rows were dropped or imputed. Zero values are coded as zero in the source files and are kept as zero. They are not missing values. Zeros in the pre-treatment earnings variables (counted by `src/00_prep_data.R`, not split by arm) are:

| Sample | `re74` = 0 | `re75` = 0 |
|---|---|---|
| NSW experimental (of 445) | 326 | 289 |
| PSID comparison (of 2,490) | 215 | 249 |
| CPS comparison (of 15,992) | 1,913 | 1,748 |

The source page does not say whether zeros for earnings mean no earnings or non-response. This is a limit to state in the write-up.

---

## 8. Known limits of these data

- The NSW benchmark uses 185 treated people. Its estimate has sampling error, and the write-up treats it as a sample estimate with a confidence interval, not as an exact truth.
- The comparison groups are from 1970s survey data (PSID and CPS). Results may not carry over to other programs or periods.
- The source page does not describe how the comparison groups were drawn beyond the variable list.
- Price base of the dollar amounts is not stated on the source page.
- The CPS comparison group has many rows at its earnings maximum (see Section 5). If the file is top-coded, earnings differences in the upper range are understated for CPS. The effect on the estimators has not been checked.
- The placebo variable `re75` is assumed to be pre-program. Enrollment timing is not in the files.

---

## 9. Reproducing the data

From the repository root:

```
Rscript data/download_lalonde.R   # downloads and checks row counts (about 3 MB)
Rscript src/00_prep_data.R        # cleans the files and writes results/data_summary.csv
```

R version used: 4.6.1 (see `results/session_info_prep.txt`). The `foreign` package (version 0.8.91) is used to read the Stata file.
