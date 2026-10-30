# Design: causal-estimators-vs-randomized-benchmark

Status: design locked on 2026-10-29, before any estimator was run.
Any change to this document after results exist must go in the "Deviations" section at the bottom, with the reason and the date.

Work already done for this design (data checks only, no estimators run):
- Downloaded the LaLonde/Dehejia-Wahba files and counted rows and treatment indicators.
- Downloaded the Hillstrom e-mail file and counted rows and arm sizes.
- Read the data-use terms on the source pages.
- Computed the pre-treatment covariate balance of the NSW experimental sample (Section 3a). This uses no outcome data.
- Installed and loaded the R packages named in Section 3c (versions in Section 2d).

No effect estimate, p-value on an outcome, or comparison of estimates has been computed yet.

---

## 1. Research question

Main question: When observational causal estimators are applied to a non-experimental comparison group, how close do they get to the answer from a randomized experiment that uses the same treated people?

Sub-questions:
1. For each estimator, is the observational estimate of the effect of job training on 1978 earnings inside the experimental confidence interval?
2. Do the estimators pass a placebo test, where the "effect" is measured on earnings from 1975, before the program started? An estimator that finds a large "effect" on 1975 earnings is picking up pre-existing differences between the groups.
3. Where a method fails, which covariate imbalance or overlap problem is responsible?
4. How strong would an unobserved confounder have to be to explain the gap between the observational estimate and the experimental benchmark?

The estimand is the average treatment effect on the treated (ATT): the average effect of the training program on 1978 real earnings for the people who were actually in the program.

---

## 2. Data

### 2a. Primary dataset: LaLonde National Supported Work (NSW) data, Dehejia-Wahba subsample

Source page: https://users.nber.org/~rdehejia/nswdata2.html (Rajeev Dehejia's data page, hosted by NBER).

Files used (direct URLs, all checked with HTTP GET on 2026-10-29):

| File | URL | Verified size | Verified rows |
|---|---|---|---|
| NSW experimental sample, Stata | https://users.nber.org/~rdehejia/data/nsw_dw.dta | 28,598 bytes | 445 (185 treated, 260 control) |
| NSW treated, text | https://users.nber.org/~rdehejia/data/nswre74_treated.txt | not committed | 185 |
| NSW control, text | https://users.nber.org/~rdehejia/data/nswre74_control.txt | not committed | 260 |
| PSID comparison group, text | https://users.nber.org/~rdehejia/data/psid_controls.txt | 400,890 bytes | 2,490 (all treat = 0) |
| CPS comparison group, text | https://users.nber.org/~rdehejia/data/cps_controls.txt | 2,574,712 bytes | 15,992 (all treat = 0) |

Row counts match the numbers stated on the source page (185, 260, 2,490, 15,992).

Variables (order as stated on the source page for the text files): treatment indicator, age, education (years), black, hispanic, married, nodegree, re74 (1974 earnings), re75 (1975 earnings), re78 (1978 earnings, the outcome). The .dta file has the same variables with a label column (`data_id`).

Total size of the three files needed for the main analysis is about 3.0 MB. This is under the 5 MB limit, but the raw files will still be downloaded by a script and not committed, so the repository stays small and the data stays with its source.

Terms of use (from the source page): "distributed for attributable non-commercial use (CC by NC)". The page asks that three papers be cited: Dehejia and Wahba (1999, JASA), Dehejia and Wahba (2002, REStat), and LaLonde (1986, AER). This is compatible with a public, non-commercial portfolio repository, provided the citations and the license note appear in the README.

Why this dataset is the right benchmark: The 185 treated people and 260 control people were randomly assigned by the NSW experiment. The difference in average 1978 earnings between these two randomized groups is an unbiased estimate of the ATT. That experimental estimate is the "known answer" for the observational comparisons. The comparison groups (PSID and CPS) were not randomized, so they are the observational data where bias can occur.

Limit to state in the write-up: The experimental estimate is itself a sample estimate with sampling error. It is not an exact truth. The confidence interval around it is part of the benchmark.

### 2b. Secondary dataset (gated): Hillstrom e-mail marketing experiment, 2008

Source page: http://blog.minethatdata.com/2008/03/minethatdata-e-mail-analytics-and-data.html (MineThatData, Kevin Hillstrom).
Data file: http://www.minethatdata.com/Kevin_Hillstrom_MineThatData_E-MailAnalytics_DataMiningChallenge_2008.03.20.csv

Verified on 2026-10-29: HTTP 200, `text/csv`, 3,964,977 bytes, 64,001 lines (one header and 64,000 rows). The https URL on the bare domain (minethatdata.com) did not respond, so the http://www. link above is the one that was verified.

Arm sizes (randomized): Mens E-Mail 21,307; No E-Mail 21,306; Womens E-Mail 21,387. Conversions: 578 out of 64,000.

Variables: recency, history_segment, history, mens, womens, zip_code, newbie, channel, segment (treatment arm), visit, conversion, spend.

Terms: The source page did not show a license statement in the portion checked. Before this component is run or committed, the license or terms must be confirmed in writing from the source. If they cannot be confirmed, this component is dropped and the README says so. The file is under 5 MB but is still downloaded by script, not committed.

### 2d. R packages (installed 2026-10-29 on R 4.6.1, all loaded successfully)

Matching 4.10.15; sensemakr 0.1.6; sensitivitymv 1.4.4; rbounds 2.2; MatchIt 4.8.1; WeightIt 2.1.0; cobalt 5.0.0; EValue 4.1.4. MatchIt, WeightIt and cobalt are used only for cross-checks of the matching and weighting code. No substitute package is needed for any named method.

### 2c. Datasets considered and not used

- Criteo uplift data: much larger than 5 MB and not needed for the question.
- Packaged copies of the LaLonde data inside R packages (for example `causaldata`; one of its files was reachable on 2026-10-29, the others were not checked): these are copies of the same underlying data. The NBER source above is used instead, so the data's origin and terms are clear.
- Synthetic control: not applicable (see Section 4). No public panel of many units over many periods with a randomized arm was identified for this design.

---

## 3. Evaluation protocol (fixed before any results)

### 3a. The benchmark

- Benchmark estimate (primary): unadjusted difference in mean re78 between the 185 NSW treated and the 260 NSW experimental controls. Standard error: Neyman (independent-groups) formula. Confidence interval: 95%. This is the estimate that randomization supports without any model.
- Benchmark estimate (secondary, reported only): regression-adjusted version of the same comparison, with the eight covariates, using the experimental sample only.
- Randomization check (pre-treatment covariates only, no outcome). Standardized mean differences (SMD) are computed for each of the eight covariates and reported. They are diagnostics, not a stop rule. The stop rule is an omnibus test: the likelihood-ratio statistic for a logit of treatment on all eight covariates, compared with its permutation distribution from 10,000 random re-assignments of treatment within the experimental sample. The project stops for investigation only if this permutation p-value is below 0.01.
- Permutation p-value for the benchmark itself, from 10,000 random re-assignments of treatment within the experimental sample.

Result of the pre-treatment check, computed on 2026-10-29 from `nsw_dw.dta` (no outcome used):

| Covariate | SMD (treated minus control, pooled SD) |
|---|---|
| age | 0.107 |
| education | 0.141 |
| black | 0.044 |
| hispanic | -0.175 |
| married | 0.094 |
| nodegree | -0.304 |
| re74 | -0.002 |
| re75 | 0.084 |

Several SMDs are above 0.1 (education, hispanic, nodegree, age, married are above or near it). This is why the per-covariate rule was dropped. Omnibus permutation test (10,000 permutations, seed 20261029): observed LR statistic 16.98 on 8 degrees of freedom, permutation p = 0.034. This is above the 0.01 stop threshold, so the project continues. It is reported as a real imbalance in the experimental sample that the adjusted benchmark (secondary) and all adjusted estimators must handle. The unadjusted benchmark remains the primary one because randomization still supports it.

### 3b. The observational estimates

- Treated group: the 185 NSW treated people (fixed in every analysis).
- Comparison groups: two, run separately: (i) PSID comparison group (2,490 people); (ii) CPS comparison group (15,992 people). Comparison groups are never mixed.
- Covariates for adjustment, fixed: age, education, black, hispanic, married, nodegree, re74, re75.
- Outcome: re78.

### 3c. Estimators (six per comparison group, fixed)

1. Naive difference in means (no adjustment). This is the reference point for how large the bias is before any adjustment.
2. Regression adjustment: OLS of re78 on the covariates, fit only on the comparison group. Predicted re78 for each treated person is the counterfactual; the ATT is the mean of (observed treated re78 minus predicted re78). This is used instead of a pooled regression with a treatment dummy, because the pooled coefficient is not an ATT when effects differ across people.
3. Nearest-neighbour matching: 1-to-1 matching with replacement on the logit of the propensity score. Caliper: 0.2 standard deviations of the logit propensity score. Units outside the caliper are dropped. ATT estimate from the matched sample. Standard error: Abadie and Imbens (2006) variance estimator, computed in R package `Matching`, as well as a 2,000-replicate bootstrap, both reported. The bootstrap is only an approximation for matching with replacement, so the Abadie-Imbens result is the primary one.
4. Inverse probability weighting (IPW), ATT form: propensity score from logistic regression on the covariates (main effects only). Comparison-group weights are odds p/(1-p). Weights normalized to sum to one (Hajek form). Propensity scores are trimmed to [0.01, 0.99] before weighting. Standard error: 2,000-replicate bootstrap.
5. Doubly robust (augmented IPW, AIPW): same propensity model as IPW, and an outcome regression for the comparison group (OLS on the covariates, fit only on comparison units). ATT estimate combines both. Standard error: influence-function formula, with bootstrap as the cross-check.
6. Difference-in-differences (DiD): change in earnings from 1975 to 1978 (re78 minus re75) for treated versus comparison group, using the same comparison group. Identification assumption: parallel trends. The check for this assumption is in Section 3e.

Synthetic control: NOT APPLIED. Reason: it requires a panel with many units observed over many periods to build a weighted control unit. The NSW-PSID and NSW-CPS data are individual-level cross sections with only two pre-program years (1974 and 1975) and one post-program year (1978). It cannot be built in a way that matches the other designs, so it is reported as not applicable rather than being forced in.

### 3d. Primary metrics and pass/fail rules (decided now)

Bias definition. Every estimator uses the same 185 NSW treated people as the benchmark does, so the treated mean cancels:

bias = (observational ATT) minus (benchmark ATT) = (NSW experimental control mean) minus (observational counterfactual mean for the treated)

So bias is measured on the counterfactual only. This removes the shared treated-group noise from the comparison. Bootstrap: 2,000 replicates, resampling the NSW treated, NSW control, and comparison group separately, with the full estimator re-run in each replicate. The 95% bootstrap confidence interval for the bias is reported in dollars and as a percentage of the benchmark mean.

Precision guard. Let H_b be the half-width of the benchmark's 95% confidence interval (Section 3a, primary) and H_o the half-width of the bias confidence interval. A cell is "informative" only if H_o is at most 2 times H_b. A cell that is not informative is marked "uninformative" and never counts as a pass.

Primary pass rule (benchmark test): a cell passes if it is informative AND the 95% bootstrap confidence interval for the bias contains zero.

Secondary metrics, reported and not used for pass/fail: the bias in dollars and percent; the observational 95% confidence interval for the ATT itself (to show whether it contains the benchmark point estimate); root mean squared error across the robustness grid (Section 3f).

Summary claim allowed: "Observational methods recovered the experimental answer" only if at least 4 of the 6 estimators pass (informative and covering zero) for at least one comparison group, AND those same estimators pass the placebo tests in 3e. Uninformative cells do not count toward the 4. Otherwise the summary says that the methods did not recover the benchmark.

No estimator will be dropped, added, or relabelled after results are in, except through the Deviations section.

### 3e. Placebo and falsification tests

- Placebo outcome for cross-sectional estimators (1-4 above): re75 (1975 earnings), treated as a pre-program outcome as the Dehejia-Wahba data classify it. The source page does not give program dates, so the timing is an assumption stated here, not verified. Because randomization is assumed to hold in the NSW sample, the experimental placebo difference is expected to be near zero. The placebo benchmark is the NSW experimental difference in re75. The same bias definition and bootstrap as 3d apply, with re75 removed from the covariate set. Pass if the placebo bias interval contains zero and the cell is informative (same guard as 3d).
- Pre-trend check for DiD: change from 1974 to 1975 (re75 minus re74). Benchmark: the same change in the NSW experimental sample. Pass under the same bias, bootstrap, and guard rules.

A method that fails the placebo test is marked "fails placebo" even if it passes the benchmark test. This is important: a method can match the benchmark by chance while being biased in a way the placebo exposes. Any partial overlap of 1975 earnings with program participation for some enrollees would contaminate the placebo; this is listed as a limit in Section 8.

### 3f. Robustness grid (fixed, all reported)

For each estimator with a tuning choice, the following variations are run and all results are reported, whichever direction they go:
- Propensity trimming: [0.01, 0.99] (primary), [0.05, 0.95].
- Matching caliper: 0.1 SD (secondary), 0.2 SD (primary), no caliper.
- Covariate set: primary set; primary set without re74 (test of how much the pre-program earnings control matters).
- Standard error method: Abadie-Imbens and bootstrap for matching; influence function and bootstrap for AIPW.

No variation is chosen as the "best" after seeing results. The primary specification is the one listed in 3c.

### 3g. Hold-out role of the benchmark

The benchmark is used to score the estimators, not to tune them. Every tuning choice above was fixed in this document before any estimate was computed. The benchmark is never used to pick a caliper, trimming level, or covariate set.

---

## 4. Analysis components

At least four distinct components are required. This design has six, plus one gated.

**Component A: Benchmark construction and randomization check.** Experimental ATT with Neyman standard error and 95% CI; permutation p-value (10,000 permutations); covariate balance table (SMD) between NSW treated and NSW control. Output: `results/benchmark.csv`, `results/balance_nsw.csv`.

**Component B: Observational estimates against the benchmark.** Six estimators (3c) for each of two comparison groups (PSID, CPS). Pass/fail against the benchmark (3d). Output: `results/observational_vs_benchmark.csv`, with bias in dollars and percent and the pass/fail flag.

**Component C: Diagnostics for failure.** Covariate balance before and after weighting or matching (SMD for each of the 8 covariates, plus SMD of the propensity score); effective sample size for IPW and AIPW (Kish formula); largest weight; number of matches and dropped units for the matching estimator; overlap plot of the propensity score distribution, treated versus comparison. For each failed cell, the write-up names the covariate with the largest remaining SMD. Output: `results/balance_after_adjustment.csv`, `results/overlap_*.png`.

**Component D: Placebo and pre-trend tests.** Section 3e. Output: `results/placebo_re75.csv`, `results/pretrend_did.csv`.

**Component E: Sensitivity analysis for hidden confounding.** Three analyses, all fixed:
1. Rosenbaum bounds (sensitivity parameter Gamma, from 1 to 3 in steps of 0.25) for the matched sample, using the Wilcoxon signed-rank test on matched differences (R package `sensitivitymv` or `rbounds`). Reports the Gamma at which the conclusion (effect different from zero) would no longer hold at the 0.05 level.
2. Robustness value (RV) from the `sensemakr` package for the regression adjustment estimator. The RV is the strength of a single hidden confounder, in partial R-squared terms, needed to reduce the estimate to zero and, separately, to move it to the benchmark value.
3. Benchmarking of the hidden confounder: the partial R-squared of each observed covariate (re75, education, re74, and others) with the outcome and with treatment. These are compared with the partial R-squared values needed to move the estimate to the benchmark. If a hidden confounder would need to be stronger than any observed covariate to explain the gap, that is reported as "the gap is not explained by a confounder as strong as the observed ones." If not, it is reported as "a confounder of observed strength could explain the gap."

Output: `results/sensitivity_gamma.csv`, `results/sensemakr_rv.csv`.

**Component F: Robustness grid.** Section 3f. Output: `results/robustness_grid.csv`, and a figure of estimates against benchmark for each specification: `results/robustness_plot.png`.

**Component G (gated on license check, Section 2b): Second randomized dataset with constructed confounding.** Real Hillstrom data only. The fixed details are below, written before any run.

- Arms used: "Womens E-Mail" (treated, the assigned e-mail) and "No E-Mail" (control). The "Mens E-Mail" arm is dropped.
- Benchmark: the full randomized difference in outcome rates, Womens arm minus No-E-Mail arm, with a Neyman 95% CI. Because assignment is randomized over all 64,000 rows, this is the benchmark for the constructed sample too.
- Primary outcome: `visit` (visited the site after the e-mail, binary). Secondary outcome: `conversion` (binary, about 0.9% of rows: 578 of 64,000). Conversion has too few events for a precise comparison, so it is reported but its results are not used for pass/fail alone.
- Constructed confounding (fixed formula). Standardize two covariates over the full 64,000 rows: z_recency = (recency - mean)/SD and z_history = (history - mean)/SD. Each row is kept in the observational sample with an independent Bernoulli probability:
  - Womens-arm rows: keep probability = plogis(beta * (z_recency + z_history)), with no intercept.
  - No-E-Mail rows: keep probability = plogis(-beta * (z_recency + z_history)), with no intercept.
  - Two fixed strengths: beta = 0.5 (mild) and beta = 1.0 (strong). Each strength is one run, with seed 20261029 for the thinning. The kept share and the covariate balance after thinning are reported.
- Selection depends only on pre-treatment covariates, so the randomized benchmark still applies to the thinned sample in principle. Any difference between observational and benchmark estimates is therefore an estimator failure under the constructed selection, not a change in the true effect.
- Estimators (same definitions as 3c, binary outcome): adjustment set = recency, history, mens, womens, newbie, channel, zip_code (categorical). Naive difference; regression adjustment (comparison-group-only linear probability model, predicted counterfactual); 1:1 matching with the same caliper rule (0.2 SD of the logit propensity score), risk difference on the matched sample; IPW; AIPW. The DiD estimator and the placebo tests are not applicable because the file has no pre-period outcome. The only falsification test for this component is the covariate balance after adjustment.
- Pass/fail uses the same bias definition, bootstrap, and precision guard as Section 3d.
- Label on every output: "constructed confounding on real randomized data". All outcomes are the real observed values. Output: `results/hillstrom_*.csv`.

---

## 5. Failure criteria

A result counts as a failure for an estimator if any of these hold (all decided now):

1. Benchmark failure: the bias confidence interval (Section 3d) does not contain zero, in an informative cell. An uninformative cell is neither a pass nor a failure; it is reported as uninformative.
2. Placebo failure: the placebo bias confidence interval does not contain zero, in an informative cell (Section 3e).
3. Overlap failure: effective sample size below 50, or any weight above 10 times the mean weight, for IPW or AIPW. Reported as a diagnostic failure, not a benchmark failure.
4. Balance failure: any |SMD| above 0.1 for a covariate after adjustment, in the primary specification.
5. Robustness failure: the estimator passes the benchmark in the primary specification but fails in more than half of the robustness grid cells. This is marked "fragile".

Project-level failure: if the omnibus randomization test (Section 3a) gives a permutation p-value below 0.01, the project stops and the problem is written up before any further step. The current value (p = 0.034, Section 3a) does not trigger this stop.

Null and negative results are expected to be possible. If every observational estimator fails the benchmark test, the README states that plainly, with the bias numbers, and does not present a "best method".

No claims will be made beyond these numbers. Statements about "how well" or "why" a method fails are limited to what the diagnostics in Component C and E show.

---

## 6. Repository layout

```
causal-estimators-vs-randomized-benchmark/
  DESIGN.md                 this file
  README.md                 plain-English write-up (question, data, method, results, limits, what failed)
  data/
    README.md               sources, terms, citations, and how to re-download
    download_lalonde.R      downloads the NSW and PSID/CPS files; checks row counts
    download_hillstrom.R    downloads the Hillstrom file (gated, Section 2b)
    raw/                    created by the download scripts; listed in .gitignore, not committed
  src/
    01_benchmark.R          Component A
    02_observational.R      Component B (six estimators, two comparison groups)
    03_diagnostics.R        Component C
    04_placebo.R            Component D
    05_sensitivity.R        Component E
    06_robustness.R         Component F
    07_hillstrom.R          Component G (gated)
    08_figures.R            figures from results/ (R, ggplot2, following the repo's plot style)
  results/                  CSV tables and PNG figures produced by src/
```

Language: R for all analysis and final figures. Each script is run from the command line in order, and the README lists the exact commands and the R package versions used (recorded with `sessionInfo()` into `results/session_info.txt`).

Reproduction: `Rscript data/download_lalonde.R` then `Rscript src/01_benchmark.R` through `Rscript src/08_figures.R`. The random seed is set to a fixed value at the top of each script that uses randomness (permutations, bootstrap, and the Hillstrom selection step).

---

## 7. Why this is hard to reproduce in a day

The benchmark can be computed in minutes, but the project takes a long time because each piece has to be checked against the others. The observational estimators need propensity models, overlap checks, and two different standard error methods (Abadie-Imbens and bootstrap) that disagree in small samples, and the comparison groups differ by a factor of about six in size. Each estimator has to be run on two comparison groups, passed through the benchmark test and the placebo test, and then through a robustness grid of about a dozen variations. The sensitivity analysis needs the Rosenbaum bounds and the partial R-squared calculations, and both need to be read carefully to say what they do and do not show. A second randomized dataset, with a constructed confounding rule, adds a separate data download, license check, and selection step. The write-up also has to explain every step in plain English for a reader without a coding background, which means each diagnostic must be interpreted, not just reported. Finally, the protocol is fixed in advance, so the work includes deciding the rules before results exist and then reporting results that may fail those rules. Doing all of that honestly, with checks on each estimator and each comparison group, takes more than a day.

---

## 8. Limits to state in the write-up

- The NSW benchmark is one experiment with 185 treated people. Its confidence interval is wide, and a "pass" means the observational estimate is inside that interval, not that it matches exactly.
- The Dehejia-Wahba subsample is a subset of the full NSW experiment (it keeps people with 1974 earnings). Results apply to this subsample.
- The comparison groups come from 1970s survey data. Results may not carry over to other programs or periods.
- Bootstrap standard errors for matching with replacement are approximate.
- Sensitivity analyses bound the effect of unobserved confounders relative to observed ones. They do not rule out every hidden confounder.
- The Hillstrom component (if run) tests a constructed confounding rule chosen by the author. It shows how the methods behave under that rule, not under real-world selection.
- The NSW experimental sample has imbalance on nodegree (SMD -0.30) and hispanic (SMD -0.17) by chance (omnibus permutation p = 0.034). Adjusted estimators must handle this, and the unadjusted benchmark is the primary one.
- The placebo outcome re75 is assumed pre-program. Dehejia-Wahba classify it that way, but the source page gives no program dates. Any overlap with enrollment would weaken the placebo test.

---

## 9. Deviations

**Deviation 1 (2026-10-29, before any estimator was run): Component G (Hillstrom e-mail data) dropped.**
- Reason: the license or terms check in Section 2b could not be completed with a confirmed answer. The MineThatData page (http://blog.minethatdata.com/2008/03/minethatdata-e-mail-analytics-and-data.html) describes the file as a 2008 data mining challenge and does not state a license or reuse terms in the text read. The only "terms" text on the page concerns blog subscriber e-mail addresses, not the dataset. Written confirmation was not obtained.
- Effect: Section 2b, Section 4 Component G, the Hillstrom parts of Section 6 (`download_hillstrom.R`, `07_hillstrom.R`), and the Hillstrom limit in Section 8 do not apply. `data/download_hillstrom.R` was not written and nothing from that file is committed.
- Data handling: the Hillstrom CSV was fetched once during the terms check and deleted without being read or counted. No Hillstrom rows were analysed.
- Consequence for the benchmark: none. Components A to F use only the NSW and PSID/CPS data, which have confirmed terms (CC BY-NC, stated on the source page).
- Re-entry: Component G can be re-added as a new deviation if written permission from the source is obtained.

**Deviation 2 (2026-10-29): data preparation script added.** `src/00_prep_data.R` was added to read the raw files, check that the Stata and text copies of the NSW sample agree, and write `data/processed/` and `results/data_summary.csv`. It is not in the Section 6 layout list. It uses no outcome values to build any sample and computes no effect estimate. The layout otherwise follows Section 6, plus `DATA.md` at the repository root and `data/README.md`.

**Deviation 3 (2026-10-29, disclosure of an outcome count, after the design was locked): one descriptive tally of the outcome by arm.**
- What happened: while checking the NSW experimental file during the data preparation step, one cross-tabulation of `re78 == 0` by `treat` was printed once. Result: 45 of 185 treated and 92 of 260 control people had zero 1978 earnings.
- What it is not: it is not an effect estimate, p-value, or comparison of estimates. No estimator was written or run. It did not touch any Section 3 to 5 choice.
- Protocol impact: none. No tuning choice, estimator, metric, pass/fail rule, or sample definition was changed in response to it. The design text above is unchanged apart from Deviations 1 and 2.
- Going forward: no further `re78` statistic is computed by arm before the protocol is run. The pre-treatment top-coding and zero counts that were added to `src/00_prep_data.R` are not split by arm.
