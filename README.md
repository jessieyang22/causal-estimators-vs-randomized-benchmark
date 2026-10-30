# Do observational causal estimators recover a randomized answer?

This project asks one question with real data. The LaLonde job-training experiment (National Supported Work, NSW) randomly assigned people to job training or to no training. Because of the random assignment, the difference in their 1978 earnings is a trustworthy answer for the people who got the training. Then we pretend we never ran the experiment. We take the same trainees and compare them with people who were never in the program, drawn from survey data (PSID or CPS). We ask: do standard observational methods get close to the experimental answer?

**Short answer: no, not reliably.** Some methods land inside the experimental interval for one comparison group but not the other, and every method that passes the benchmark test fails a placebo test. The full numbers are below. The summary claim that the observational methods "recovered the experimental answer" is **not** supported.

---

## 1. Question

- **Main question:** When observational causal estimators are applied to a non-experimental comparison group, how close do they get to the answer from a randomized experiment that uses the same treated people?
- **Estimand:** the average effect of the training program on 1978 real earnings for the people who were in the program (the average treatment effect on the treated, ATT).
- **Sub-questions** (from `DESIGN.md`):
  1. Is each observational estimate inside the experimental confidence interval?
  2. Do the estimators pass a placebo test, using 1975 earnings (before the program)?
  3. Where a method fails, which covariate imbalance or overlap problem is responsible?
  4. How strong would an unobserved confounder have to be to explain the gap?

The protocol (estimators, pass/fail rules, robustness grid) was written in `DESIGN.md` before any outcome was compared. Implementation rules the design left open were fixed in `DESIGN.md` Deviation 4 before the benchmark was run, and the other deviations are listed there with dates. During the run, a few code errors were fixed (for example, the Rosenbaum pairs were rebuilt from `Matching::Match`'s own output after a sanity check showed a mismatch, and the `sensemakr` direction argument was corrected after an error). Those fixes changed no protocol choice and no pass/fail rule.

---

## 2. Data

All data are real, public, and were downloaded by script. No simulated data are used.

| Sample | Rows | Role |
|---|---|---|
| NSW experimental sample (Dehejia-Wahba subsample) | 445: 185 treated, 260 control | Benchmark (randomized) |
| PSID comparison group | 2,490 | Observational comparison group 1 |
| CPS comparison group | 15,992 | Observational comparison group 2 |

- **Source:** Rajeev Dehejia's data page, hosted by NBER: https://users.nber.org/~rdehejia/nswdata2.html
- **Terms (as printed on the source page):** "distributed for attributable non-commercial use (CC by NC)."
- **Citations requested by the source page:** Dehejia and Wahba (1999, *Journal of the American Statistical Association*); Dehejia and Wahba (2002, *Review of Economics and Statistics*); LaLonde (1986, *American Economic Review*). `[Full article titles to be added after checking the journals; they were not copied from the source page.]`
- **Variables:** treatment indicator, age, education (years), black, hispanic, married, no high-school degree, 1974 earnings (re74), 1975 earnings (re75), 1978 earnings (re78, the outcome). Dollar amounts are as stated in the files; the price year is not given.
- **Not committed:** the raw files (about 3 MB in total) are downloaded by `data/download_lalonde.R` into `data/raw/`, which is gitignored. Data notes and checks are in `DATA.md`.
- **Not used:** a second randomized dataset (Hillstrom e-mail test). Its terms could not be confirmed, so it was dropped before any analysis. See `DESIGN.md`, Deviation 1.

---

## 3. Method, step by step

Plain-English version of each step. Each script is in `src/`.

**Step 1: Benchmark (`src/01_benchmark.R`).** Take the average 1978 earnings of the 185 treated people and subtract the average of the 260 randomized controls. This is the experimental answer. Its 95% interval uses the standard formula for two independent groups. A permutation test (10,000 random reshuffles of who was treated) gives a p-value. We also check that the randomization looks sound: the eight pre-program covariates are compared between the two arms (standardized mean differences), and a joint test asks whether treatment is related to the covariates more than chance allows.

**Step 2: Six observational estimators (`src/02_observational.R`).** Each one uses the same 185 treated people and one comparison group at a time.

1. *Naive:* difference in average 1978 earnings, with no adjustment.
2. *Regression adjustment:* fit a model of 1978 earnings on the eight covariates using only the comparison group. Predict what each treated person would have earned without the program, and average the gap.
3. *Nearest-neighbour matching:* for each treated person, find the comparison person with the closest propensity score (the estimated chance of being in the program, from the eight covariates). Caliper 0.2 standard deviations, matching with replacement, using the `Matching` package.
4. *Inverse probability weighting (IPW):* weight comparison people by their odds of being treated, so that the comparison group looks more like the treated group.
5. *Doubly robust (AIPW):* combine the regression and the weighting, so that the answer is right if either model is right.
6. *Difference-in-differences (DiD):* the treated group's change from 1975 to 1978 minus the comparison group's change over the same years.

**Step 3: Uncertainty.** The "bias" of each estimator is its estimate minus the benchmark. To get a 95% interval for the bias, the whole calculation (including the propensity model) is repeated 2,000 times on resampled data, with the treated, NSW control and comparison people each resampled with replacement. A method passes the benchmark test if this bias interval contains zero and the interval is not too wide to be informative (its half-width is at most twice the benchmark's half-width).

**Step 4: Diagnostics (`src/03_diagnostics.R`).** For each weighted or matched estimator: how balanced the covariates are before and after adjustment, the effective sample size (how many comparison people still count after weighting), the largest weight, and overlap plots of the propensity scores.

**Step 5: Placebo tests (`src/04_placebo.R`).** The same estimators are run on 1975 earnings, which happened before the program, with 1975 earnings removed from the covariates. An estimator that finds a large "effect" on 1975 earnings is picking up differences that were already there. For DiD, the test is the 1974-to-1975 change. A method has to pass both the benchmark test and its placebo test to count toward the summary claim.

**Step 6: Sensitivity to hidden confounding (`src/05_sensitivity.R`).** Rosenbaum bounds ask how large a hidden bias would have to be, in odds terms, before the matched result stops looking non-zero. The robustness value (from `sensemakr`) asks how strong an unobserved confounder would have to be to shrink the estimate to zero, or to move it to the benchmark.

**Step 7: Robustness grid (`src/06_robustness.R`).** Re-run the estimators with a looser propensity cap (0.05), a 0.1 SD or no caliper, and without the 1974 earnings covariate. Every grid cell gets the same bootstrap rule as the main results. A method is labelled "fragile" if it passes the main test but fails in more than half of its grid cells.

**Step 8: Figures (`src/08_figures.R`).** Plots of each estimate against the benchmark, and the session information file.

---

## 4. Results

All dollar amounts are 1978 real earnings, as stated in the source files. Bias is estimate minus benchmark. All numbers below come from the CSV files in `results/`, produced by the scripts above.

### 4a. Benchmark (experimental answer)

| Quantity | Estimate | 95% interval |
|---|---|---|
| NSW experimental difference in 1978 earnings (primary) | **$1,794** (SE $671) | $479 to $3,109 |
| Same, regression-adjusted (secondary, reported only) | $1,676 (SE $639) | $425 to $2,928 |
| Permutation p-value (10,000 reshuffles) | 0.0035 | |
| Placebo benchmark, 1975 earnings | $265 (SE $305) | −$333 to $863 |
| Pre-trend benchmark, 1975 minus 1974 | $277 (SE $384) | −$476 to $1,029 |

Randomization check. The standardized mean differences (treated minus control, pooled SD) are: age 0.11, education 0.14, black 0.04, hispanic −0.17, married 0.09, no degree −0.30, 1974 earnings −0.00, 1975 earnings 0.08. The joint test (likelihood-ratio statistic 16.98 on 8 degrees of freedom, permutation p = 0.034) is above the pre-set stop threshold of 0.01, so the project continued, as `DESIGN.md` specifies. Imbalance on "no degree" is the largest, and it is the one the adjusted estimators have to handle.

### 4b. Six estimators against the benchmark (primary specification)

Pass = the bias interval contains zero and the cell is informative.

**PSID comparison (2,490 people)**

| Estimator | ATT estimate | 95% interval | Bias | Bias 95% interval | Pass? |
|---|---|---|---|---|---|
| 1 Naive | −$15,205 | −$16,493 to −$13,917 | −$16,999 | −$17,909 to −$16,103 | No |
| 2 Regression | $791 | −$769 to $2,350 | −$1,004 | −$2,342 to $340 | **Yes** |
| 3 Matching | $2,022 | −$757 to $4,801 | +$228 | −$2,552 to $1,255 | **Yes** |
| 4 IPW | −$183 | −$2,097 to $1,731 | −$1,977 | −$3,823 to −$265 | No |
| 5 AIPW | $1,925 | $273 to $3,577 | +$130 | −$1,379 to $1,539 | **Yes** |
| 6 DiD | $2,327 | $1,085 to $3,568 | +$532 | −$398 to $1,447 | **Yes** |

**CPS comparison (15,992 people)**

| Estimator | ATT estimate | 95% interval | Bias | Bias 95% interval | Pass? |
|---|---|---|---|---|---|
| 1 Naive | −$8,498 | −$9,641 to −$7,354 | −$10,292 | −$10,953 to −$9,604 | No |
| 2 Regression | $690 | −$523 to $1,903 | −$1,104 | −$1,952 to −$280 | No |
| 3 Matching | $1,693 | $233 to $3,153 | −$101 | −$1,673 to $856 | **Yes** |
| 4 IPW | −$3,549 | −$4,707 to −$2,391 | −$5,343 | −$6,076 to −$4,579 | No |
| 5 AIPW | $1,022 | −$185 to $2,229 | −$773 | −$1,629 to $74 | **Yes** |
| 6 DiD | $3,621 | $2,437 to $4,806 | +$1,827 | +$1,027 to +$2,639 | No |

Bias as a percent of the benchmark (PSID / CPS): naive −947% / −574%; regression −56% / −62%; matching +13% / −6%; IPW −110% / −298%; AIPW +7% / −43%; DiD +30% / +102%.

All six primary cells were informative. Two other cells were not, and they are never counted as passes: the PSID placebo matching cell (Section 4c) and the PSID matching cell without 1974 earnings (Section 4g).

### 4c. Placebo tests (1975 earnings, before the program)

The benchmark for these tests is $265 (1975 earnings, experimental) and $277 (1974 to 1975 change, experimental).

| Estimator | PSID: placebo bias (95% interval) | PSID: pass? | CPS: placebo bias (95% interval) | CPS: pass? |
|---|---|---|---|---|
| 1 Naive | −$17,796 (−$18,424 to −$17,120) | No | −$12,384 (−$12,780 to −$11,963) | No |
| 2 Regression | −$1,938 (−$2,982 to −$950) | No | −$1,635 (−$2,389 to −$902) | No |
| 3 Matching | −$4,167 (−$5,498 to −$1,875), uninformative | No | −$1,405 (−$2,675 to −$611) | No |
| 4 IPW | −$4,945 (−$6,161 to −$3,780) | No | −$6,958 (−$7,437 to −$6,434) | No |
| 5 AIPW | −$2,662 (−$3,814 to −$1,557) | No | −$1,608 (−$2,358 to −$890) | No |
| 6 DiD pre-trend | −$475 (−$1,091 to $83) | **Yes** | −$474 (−$1,018 to $8) | **Yes** |

Every cross-sectional estimator (1 to 5) fails the 1975 placebo in both comparison groups. A negative placebo bias means that, after the method's adjustment, the comparison group's 1975 earnings are higher than the trainees' actual 1975 earnings, so the methods predict the trainees should have earned more before the program than they did. Only the DiD pre-trend test passes, and DiD's benchmark test passes only for PSID.

### 4d. Summary claim (rule from `DESIGN.md` 3d)

The claim "observational methods recovered the experimental answer" needs at least 4 of 6 estimators to pass both the benchmark test and the placebo test, for one comparison group.

- **PSID:** 4 of 6 pass the benchmark test (regression, matching, AIPW, DiD), but only 1 of 6 passes both tests (DiD).
- **CPS:** 2 of 6 pass the benchmark test (matching, AIPW), and none passes both.

**The summary claim is not met for either group.** The benchmark passes are not supported by the placebo tests.

### 4e. Why methods fail (diagnostics, `results/diagnostics_weights_matching.csv`, `results/balance_after_adjustment.csv`)

- **Weighting (IPW and AIPW) has an overlap problem in both groups.** Effective sample size is 51 of 2,490 comparison people for PSID (2.1%), and the largest weight is 133 times the average. For CPS it is 1,209 of 15,992 (7.6%), and the largest weight is 48 times the average. Both fail the largest-weight rule in `DESIGN.md` 5 (a weight above 10 times the mean). The PSID effective sample size (51) is just above the 50 threshold. After IPW, the largest remaining imbalance is on the propensity score itself (SMD 0.26 for PSID, 1.39 for CPS). Among the eight covariates, the largest remaining imbalance is no-degree (0.43 for CPS). The weighted PSID estimate rests on very few comparison people.
- **Matching balances the covariates better but not perfectly.** For PSID, all 185 treated people found a match within the caliper, using 171 distinct comparison people; the largest remaining imbalance is married (SMD 0.32). For CPS, all 185 found a match, using 1,046 distinct people; the largest remaining imbalance is black (SMD 1.12). The CPS matches are still far from balanced on race, which is the main reason to distrust the CPS matching estimate even though it passes the benchmark test. The caliper never dropped a treated person, so the caliper and no-caliper results are identical.
- **Regression and DiD have no weights.** Their failures come from the model: regression depends on the outcome model being right, and DiD depends on comparison people changing like the treated people, which the placebo test does not support.

### 4f. Sensitivity to hidden confounding (`src/05_sensitivity.R`)

- **Rosenbaum bounds (matched pairs, 185 pairs).** For PSID, the matched result is significant at Gamma = 1 (upper-bound p = 0.004) and stops being significant at the 5% level between Gamma = 1 and Gamma = 1.25 (upper-bound p = 0.076 at 1.25). For CPS, the matched-pair signed-rank test is not significant even at Gamma = 1 (upper-bound p = 0.07), so the matched CPS estimate is not robust to any hidden bias. These are the p-values returned by `rbounds::psens`. Matching with replacement reuses comparison people, so these pairs are not independent, which means the bounds are approximate.
- **Robustness value (pooled regression, `sensemakr`).** The pooled regression is used because `sensemakr` needs a treatment coefficient. Its estimate is $752 for PSID and $699 for CPS. An unobserved confounder would need to explain only 1.6% (PSID) and 1.0% (CPS) of the remaining variation in both treatment and earnings to shrink that estimate to zero. Moving it to the benchmark ($1,794) needs 2.2% (PSID) and 1.6% (CPS). The estimate is not significant at the 5% level even before any confounding (the significance robustness value is 0 in both groups).
- **Benchmarking against observed covariates.** The 1975 earnings variable has the largest partial R-squared with the outcome (0.16 in PSID, 0.09 in CPS). Its partial R-squared with treatment is small (0.003 PSID, 0.001 CPS). If a confounder as strong as 1975 earnings were left out, the pooled PSID estimate would move to about −$266, and the pooled CPS estimate to about +$46. This shows that the pooled PSID estimate can change sign under a confounder of 1975 earnings' observed strength. It does not show that such a confounder exists.

### 4g. Robustness grid (`results/robustness_grid.csv`, `results/robustness_summary.csv`)

Each row is one variant of one estimator. "Pass" uses the same rule as the primary cells (bias interval contains zero and the cell is informative). Naive and DiD have no tuning choice or covariate-set variant, so they have no grid rows.

| Comparison | Variant | Estimator | ATT (95% interval) | Bias (95% bootstrap interval) | Pass? |
|---|---|---|---|---|---|
| PSID | Propensity cap 0.05 | IPW | −$5,566 (−$7,488 to −$3,643) | −$7,360 (−$9,200 to −$5,517) | No |
| PSID | Propensity cap 0.05 | AIPW | $1,541 ($92 to $2,991) | −$253 (−$1,559 to $1,044) | Yes |
| PSID | Caliper 0.1 SD | Matching | $2,022 (−$757 to $4,801) | +$228 (−$2,552 to $1,242) | Yes |
| PSID | No caliper | Matching | $2,022 (−$757 to $4,801) | +$228 (−$2,542 to $1,255) | Yes |
| PSID | Without 1974 earnings | Regression | $443 (−$1,156 to $2,042) | −$1,351 (−$2,761 to $25) | Yes |
| PSID | Without 1974 earnings | Matching | −$717 (−$4,135 to $2,701) | −$2,511 (−$5,818 to $504), uninformative | No |
| PSID | Without 1974 earnings | IPW | −$1,259 (−$3,555 to $1,036) | −$3,054 (−$5,400 to −$1,093) | No |
| PSID | Without 1974 earnings | AIPW | $637 (−$1,357 to $2,632) | −$1,157 (−$3,233 to $585) | Yes |
| CPS | Propensity cap 0.05 | IPW | −$7,148 (−$8,280 to −$6,016) | −$8,943 (−$9,611 to −$8,244) | No |
| CPS | Propensity cap 0.05 | AIPW | $788 (−$392 to $1,968) | −$1,007 (−$1,850 to −$195) | **No** |
| CPS | Caliper 0.1 SD | Matching | $1,693 ($233 to $3,153) | −$101 (−$1,673 to $856) | Yes |
| CPS | No caliper | Matching | $1,693 ($233 to $3,153) | −$101 (−$1,673 to $856) | Yes |
| CPS | Without 1974 earnings | Regression | $625 (−$574 to $1,824) | −$1,169 (−$1,989 to −$353) | No |
| CPS | Without 1974 earnings | Matching | $1,189 (−$271 to $2,649) | −$605 (−$1,710 to $865) | Yes |
| CPS | Without 1974 earnings | IPW | −$3,575 (−$4,735 to −$2,416) | −$5,370 (−$6,113 to −$4,620) | No |
| CPS | Without 1974 earnings | AIPW | $936 (−$259 to $2,131) | −$859 (−$1,692 to −$6) | **No** |

The robustness summary (`results/robustness_summary.csv`) applies the fragility rule from `DESIGN.md` 5:

- **PSID matching:** passes in the primary specification and in 2 of 3 grid cells. The third cell (no 1974 earnings) is uninformative, so it counts as not passing, but that is 1 of 3, not more than half. Not fragile.
- **PSID AIPW:** passes in the primary specification and in both grid cells. Not fragile.
- **PSID regression:** passes in the primary specification and in its one grid cell. Not fragile.
- **CPS matching:** passes in the primary specification and in all three grid cells. Not fragile.
- **CPS AIPW:** passes in the primary specification, then fails in both grid cells (propensity cap 0.05, and no 1974 earnings). **Fragile.**
- **Failing estimators** (naive, IPW in both groups, CPS regression and CPS DiD) are not labelled fragile, because they already fail the primary test.

Root mean squared error of the bias across the primary and grid cells: PSID regression $1,190, matching $1,271, IPW $4,740, AIPW $688, DiD $532 (primary only), naive $16,999. CPS regression $1,137, matching $315, IPW $6,766, AIPW $885, DiD $1,827 (primary only), naive $10,292.

The grid does not change the conclusion. Matching is the estimator that passes the benchmark test in both groups and in almost every grid cell, but it fails the 1975 placebo in both groups (Section 4c), so it does not meet the summary-claim rule. Its covariate balance on race is also poor (Section 4e).

### 4h. Figures

- `results/overlap_psid.png`, `results/overlap_cps.png`: propensity score overlap, treated versus comparison.
- `results/robustness_plot.png`: each estimate and its interval against the benchmark, for each grid cell.

---

## 5. What failed

- **The summary claim.** Neither comparison group has 4 of 6 estimators passing both the benchmark test and the placebo test.
- **The naive estimator** is off by about $17,000 (PSID) and $10,000 (CPS), which is roughly 6 to 10 times the size of the benchmark.
- **Inverse probability weighting** fails the benchmark test in both groups, and its effective sample sizes are very small (see 4e).
- **Every cross-sectional estimator fails the 1975 placebo** in both groups, even the estimators that pass the 1978 benchmark. A method that passes the benchmark by chance while failing the placebo cannot be trusted for this question.
- **The CPS doubly robust (AIPW) estimate is fragile.** It passes the benchmark test in the primary specification, but it fails with a propensity cap of 0.05 and without 1974 earnings. Its bias is about −$770 to −$1,000 in those variants, and the intervals exclude zero.
- **The CPS difference-in-differences estimate** is about $1,800 too high, and its interval does not contain zero.
- **Regression adjustment for CPS** is about $1,100 too low, and its interval does not contain zero.
- **The CPS matching estimate** passes the benchmark test, but its covariate balance on race is poor (SMD 1.12), and the matched-pair test is not significant even at Gamma = 1.
- **The randomization check** showed an omnibus imbalance (p = 0.034) in the experimental sample. It did not reach the stop threshold, but the adjusted estimators must handle it.
- **Synthetic control** was not applied. The data have only two pre-program years and no panel of units, so the design could not be built as specified (`DESIGN.md`, Section 3c).

---

## 6. Limits

- The benchmark is one experiment with 185 treated people. Its interval is wide ($479 to $3,109), so a "pass" means the observational estimate is inside that interval, not that it matches exactly.
- The Dehejia-Wahba subsample keeps people with 1974 earnings. Results apply to that subsample.
- The comparison groups come from 1970s survey data. Results may not carry over to other programs or periods.
- The placebo outcome (1975 earnings) is assumed to be before enrollment. The source page does not give program dates, and any overlap with enrollment would weaken the placebo test.
- The bootstrap intervals for matching with replacement are approximate. The influence-function standard error for AIPW treats the fitted models as known; the bootstrap intervals in the tables do not rely on that.
- The Rosenbaum bounds assume independent pairs. Matching with replacement breaks that assumption.
- The robustness value uses a pooled regression, not the comparison-only regression estimate. It shows sensitivity of the pooled coefficient only.
- Sensitivity analyses bound the effect of unobserved confounders relative to observed ones. They do not rule out every hidden confounder.
- Dollar amounts have no stated price year in the source.
- The CPS file has many rows at its earnings maximum (19% at the 1974 and 1975 maxima, 26% at the 1978 maximum), which looks like top-coding. This has not been checked against any other documentation. Effects on the estimators have not been checked (`DATA.md`, Section 5).
- Zeros in earnings are kept as zeros. The source does not say whether zero means no earnings or non-response.

---

## 7. Reproduce

From the repository root:

```
Rscript data/download_lalonde.R      # downloads the raw files (about 3 MB) into data/raw/
Rscript src/00_prep_data.R           # cleans the files, writes data/processed/ and results/data_summary.csv
Rscript src/01_benchmark.R           # Component A
Rscript src/02_observational.R       # Component B (about 10 minutes on 8 cores)
Rscript src/03_diagnostics.R         # Component C
Rscript src/04_placebo.R             # Component D (about 25 minutes on 8 cores)
Rscript src/05_sensitivity.R         # Component E
Rscript src/06_robustness.R          # Component F (about 40 minutes on 8 cores)
Rscript src/08_figures.R             # figures and results/session_info.txt
```

Component G (`src/07_hillstrom.R`) is not included (Deviation 1). The random seed is 20261029. Each bootstrap replicate uses seed 20261029 plus its replicate number, so results do not depend on the number of cores used. R version 4.6.1. Package versions are in `results/session_info.txt` and `results/session_info_prep.txt`.

---

## 8. Files

- `DESIGN.md`: the locked protocol, with all dated deviations.
- `DATA.md`: data record, file checksums, variable ranges, and missingness.
- `data/README.md`: data sources and how to re-download.
- `src/helpers.R`: shared estimator, bootstrap and balance code.
- `src/00_prep_data.R` to `src/08_figures.R`: the analysis, in order.
- `results/`: all CSV tables and PNG figures.
