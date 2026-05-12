# CNED EEG Analysis Pipeline (R)

**Author:** Dryden Arseneau  
**Affiliations:** Schabrun Lab · Seminowicz Lab  
**Dataset:** CNED – Zhao et al., Sci Data (2025)  
https://doi.org/10.1038/s41597-025-05900-1

---

## 1. Overview

This R pipeline takes the outputs of the MATLAB preprocessing/spectral pipeline and the Python source localization pipeline, merges them with behavioural and demographic data, and fits a comprehensive set of statistical models. Two parallel analytical streams — **channel-space** (scalp EEG) and **source-space** (sLORETA ROI time courses) — run side by side, enabling a formal comparison of whether localizing alpha oscillations to anatomical sources improves the prediction of pain perception.

The primary statistical tool is the **generalized additive mixed model (GAMM)**, which is well-suited to trial-level pain data because it: (a) accommodates nonlinear predictor–outcome relationships without imposing a functional form, (b) handles deeply nested random effects (trials within subjects within experiments), and (c) supports multiple simultaneous predictors with adaptive smoothness selection via `select = TRUE`.

Classical inferential tests (paired t-tests, repeated-measures ANOVA, Rayleigh test) complement the GAMMs by answering questions the GAMMs do not address directly — effect existence, spatial distribution across ROIs, and phase clustering.

---

## 2. Pipeline Architecture

The pipeline has two merge streams feeding into parallel model-fitting streams, plus a cross-stream comparison and classical tests layer.

```
MATLAB spectral outputs               Python source outputs
(chan_trial_*.csv per subject)        (sub-XXX_source_trial.csv per subject)
            ↓                                       ↓
  merge_behavioural.R                    merge_source_spectral.R
  merge_participants_into_behavioural.R  merge_source_ga.R
  merge_spectral_behaviour.R                       ↓
            ↓                              source_pain_master.csv
  alpha_pain_master.csv                  source_ga_master.csv
  behavioural_demo_master.csv            source_ga_fooof_master.csv
            ↓                                       ↓
  run_gamm_alpha_metrics_v2.R          run_gamm_source.R
  gamm_outputs_v2/                     gamm_outputs_source/<ROI>/
            ↓                                       ↓
            └──────────────────┬────────────────────┘
                               ↓
                  compare_channel_source_gamm.R
                  run_classical_tests.R
```

---

## 3. File Reference

### Merge Layer

| File | Input | Output | Role |
|---|---|---|---|
| `merge_behavioural.R` | Per-experiment behavioural CSVs | `behavioural_master.csv` | Standardizes and concatenates raw behavioural data across all experiments |
| `merge_participants_into_behavioural.R` | `behavioural_master.csv`, `participants.tsv`, `cap_size.csv` | `behavioural_demo_master.csv` | Adds age, sex, and EEG cap size; assigns `global_subjid` and `subjid_uid` cross-experiment identifiers |
| `merge_spectral_behaviour.R` | MATLAB `chan_trial_*.csv` files, `behavioural_demo_master.csv` | `alpha_pain_master.csv` | Merges channel-space spectral metrics with behavioural/demographic master; upserts on `(experiment_id, subjid, trial)` |
| `merge_source_spectral.R` | Python `sub-XXX_source_trial.csv` files, `behavioural_demo_master.csv`, `sub-XXX_source_ga_fooof.csv` | `source_pain_master.csv` | Same merge pattern for source-space data; broadcasts GA FOOOF metrics to every trial row so FOOOF-controlled GAMMs have access to aperiodic parameters |
| `merge_source_ga.R` | Python `sub-XXX_source_ga.csv` and `sub-XXX_source_ga_fooof.csv` files | `source_ga_master.csv`, `source_ga_fooof_master.csv` | Collates grand-average source metrics (TVI_alpha, ITC, FOOOF) across subjects into group-level tables |

### Model-Fitting Layer

| File | Input | Output | Role |
|---|---|---|---|
| `run_gamm_alpha_metrics_v2.R` | `alpha_pain_master.csv` | `gamm_outputs_v2/` | Fits the channel-space GAMM model family; 14 models per dataset ranging from baseline covariates through FOOOF-controlled alpha metrics |
| `run_gamm_source.R` | `source_pain_master.csv` | `gamm_outputs_source/<ROI>/` | Fits the same GAMM model family in source space, iterating over every ROI; produces per-ROI comparison tables and a cross-ROI summary |

### Comparison and Classical Tests Layer

| File | Input | Output | Role |
|---|---|---|---|
| `compare_channel_source_gamm.R` | `gamm_outputs_v2/`, `gamm_outputs_source/` | `gamm_comparison/` | Three-layer comparison: AIC/deviance fit quality, smooth-term significance, and partial-effect concordance (shape correlation) between channel and source estimates |
| `run_classical_tests.R` | `source_pain_master.csv`, `alpha_pain_master.csv`, fitted model RDS files | `classical_tests/` | Five targeted inferential tests: ERD significance, channel vs source R², rmANOVA on metrics across ROIs, ANOVA on deviance across ROIs, Rayleigh phase uniformity |

---

## 4. Key Data Files

| File | Rows | Description |
|---|---|---|
| `behavioural_demo_master.csv` | 1 per (experiment, subject, trial) | Pain ratings, laser power, trial index, age, sex, cap size, cross-experiment IDs |
| `alpha_pain_master.csv` | 1 per (experiment, subject, trial) | Channel-space alpha metrics + behavioural data |
| `source_pain_master.csv` | 1 per (experiment, subject, trial, ROI) | Source-space trial metrics + behavioural data + FOOOF broadcast |
| `source_ga_master.csv` | 1 per (experiment, subject, ROI) | GA source metrics: TVI_alpha, ITC, pre/post spectral, LEP |
| `source_ga_fooof_master.csv` | 1 per (experiment, subject, ROI) | GA FOOOF aperiodic parameters |

All merge scripts use an **anti-join upsert pattern**: existing rows for the incoming subjects are removed and replaced, while rows for subjects not in the current batch are preserved. This means the master files can be grown incrementally as new subjects are processed without duplication.

---

## 5. GAMM Model Family

Both `run_gamm_alpha_metrics_v2.R` and `run_gamm_source.R` fit the same conceptual model family. All continuous predictors are z-scored. All models include laser power, trial index, and age as smooth covariates, plus sex and cap size as parametric terms, with crossed random effects for subject and experiment.

| Model | Predictor added | Scientific question |
|---|---|---|
| m00 | Baseline (covariates only) | Covariate-only fit |
| m01 | `s(BI_pre)` | Does sub-band balance index predict pain? |
| m02 | `s(LR_pre)` | Does log ratio predict pain? |
| m03 | `s(CoG_pre)` | Does peak alpha frequency (CoG) predict pain? |
| m04 | `s(psi_cog)` | Does the BI × CoG interaction predict pain? |
| m05 | `te(pow_slow, pow_fast)` | Does the raw 2-D slow×fast space predict pain? |
| m06 | `sin_phase + cos_phase` | Does slow-alpha phase at stimulus onset predict pain? |
| m07 | `s(delta_ERD)` | Does sub-band ERD asymmetry predict pain? |
| m08 | `s(n2p2_amp)` | Does LEP amplitude predict pain? |
| m09 | `s(BI_pre) + s(delta_ERD)` | Do pre and post-stimulus alpha carry independent variance? |
| m10 | `s(BI_pre) + s(n2p2_amp)` | Does combining alpha state with LEP improve prediction? |
| m11 | `s(BI_pre) + s(delta_ERD) + s(n2p2_amp)` | Full combined model |
| m12 | FOOOF baseline | Baseline with aperiodic terms controlled |
| m13 | FOOOF + `s(BI_pre)` | Does BI_pre predict pain independently of the aperiodic slope? |
| m14 | FOOOF + `s(LR_pre)` | Does LR_pre predict pain independently of the aperiodic slope? |

**Implementation notes:**
- All `bam()` calls use `method = "fREML"`, `discrete = TRUE`, `nthreads = 4`, `select = TRUE`
- `select = TRUE` adds a shrinkage penalty to every smooth, allowing terms to be shrunk to zero if unsupported by the data — this is the correct mechanism for adaptive term selection in mgcv
- k-values are computed via `safe_k()` from the actual number of unique values per predictor to avoid overfitting with small N
- FOOOF columns are broadcast from GA to trial level in `merge_source_spectral.R`, allowing FOOOF-controlled trial-level models

---

## 6. Classical Tests

`run_classical_tests.R` runs five targeted tests that address questions outside the GAMM framework.

### Test 1 — ERD Significance (paired t-test)
**Question:** Does alpha power actually change after the stimulus?  
**Method:** Per-subject mean slow/fast/total alpha power in pre- vs post-stimulus windows; paired t-test across subjects; FDR correction across all ROI × band combinations.  
**Output:** `test1_erd_ttest.csv`, `test1_erd_plot.png`

### Test 2 — Channel vs Source GAMM Fit (paired t-test)
**Question:** Does source localization produce a better-fitting GAMM than channel-space metrics?  
**Method:** Per-subject R² (from `cor(fitted, observed)²`) computed from each domain's trial-level fitted CSVs; paired across subjects for each matched model concept; best source ROI selected as comparator.  
**Output:** `test2_fit_ttest.csv`, `test2_fit_plot.png`

### Test 3 — Metrics Across ROIs (repeated-measures ANOVA)
**Question:** Are BI_pre, CoG_pre, and delta_ERD spatially uniform or concentrated in specific ROIs?  
**Method:** Subject-level means per ROI as the within-subject factor; afex::aov_ez with Greenhouse-Geisser correction (falls back to base `aov()` if afex unavailable); FDR-corrected pairwise post-hoc.  
**Output:** `test3_rmanova_<metric>.csv`, `test3_posthoc_<metric>.csv`, `test3_rmanova_plot.png`

### Test 4 — Deviance Explained Across ROIs (one-way ANOVA)
**Question:** Does any source ROI significantly outperform the others in GAMM fit quality?  
**Method:** Per-subject R² per ROI from source fitted CSVs; one-way ANOVA with ROI as factor, separately per model concept; FDR-corrected pairwise post-hoc.  
**Output:** `test4_roi_anova_<model>.csv`, `test4_posthoc_<model>.csv`, `test4_roi_anova_plot.png`

### Test 5 — Phase Uniformity (Rayleigh test)
**Question:** Is slow-alpha phase at stimulus onset clustered (non-uniform), as required for phase to be a meaningful pain predictor?  
**Method:** Per-subject Rayleigh test using the Mardia-Jupp approximation (self-contained implementation, no `circular` package required); applied to both `slow_phase` (pre-stim) and `slow_phase_post`; FDR correction across subject × ROI cells.  
**Output:** `test5_rayleigh_<phase_col>.csv`, `test5_rayleigh_summary_<phase_col>.csv`, `test5_rayleigh_plot_<phase_col>.png`

---

## 7. Channel vs Source Comparison

`compare_channel_source_gamm.R` provides a three-layer formal comparison.

**Layer 1 — Fit quality:** AIC, BIC, and deviance explained side-by-side for matched model concepts. For each concept, the source comparator is the best ROI (lowest AIC), making the comparison fair. Positive ΔAIC (channel − source) means source fits better.

**Layer 2 — Effect significance:** Full smooth-term tables (EDF, F-statistic, p-value) extracted from every fitted model in both domains. Allows you to ask: is a term significant in channel space, source space, both, or neither?

**Layer 3 — Smooth concordance:** Partial effect curves from matched channel and source models are interpolated onto the same x-grid and correlated. A high positive r (> 0.7) means both domains agree on the direction and shape of the effect — strong evidence the metric is capturing a genuine biological signal. The concordance heatmap (ROI × metric) visualises where in the brain the channel estimates are best reproduced.

**Output files:**
- `fit_quality_comparison.csv`
- `effect_significance_comparison.csv`
- `smooth_concordance.csv`
- `fit_quality_plot.png` — deviance explained bar chart (channel vs source)
- `delta_aic_plot.png` — signed ΔAIC per model concept
- `smooth_overlay_<metric>.png` — channel vs all source ROI smooth curves on one figure
- `concordance_heatmap.png` — ROI × metric correlation heatmap
- `comparison_report.txt` — plain-text narrative summary with domain recommendation

---

## 8. Run Sequence

Run scripts in this order. Steps 1–3 only need to be run once per dataset; steps 4–8 can be re-run as new subjects are added.

```
1.  merge_behavioural.R                  → behavioural_master.csv
2.  merge_participants_into_behavioural.R → behavioural_demo_master.csv
3.  merge_spectral_behaviour.R           → alpha_pain_master.csv        (channel stream)
4.  merge_source_spectral.R              → source_pain_master.csv       (source stream)
5.  merge_source_ga.R                    → source_ga_master.csv + fooof_master.csv
6.  run_gamm_alpha_metrics_v2.R          → gamm_outputs_v2/
7.  run_gamm_source.R                    → gamm_outputs_source/
8.  compare_channel_source_gamm.R        → gamm_comparison/
9.  run_classical_tests.R               → classical_tests/
```

**Important sequencing constraint:** `merge_source_spectral.R` (step 4) reads the FOOOF GA files (`sub-XXX_source_ga_fooof.csv`) that are produced by the Python source pipeline. Run the Python pipeline with `fooof.enabled = true` first, then run the merge. If FOOOF files are absent, the merge proceeds without FOOOF columns and a warning is printed — models m12–m14 in `run_gamm_source.R` will run but the FOOOF predictors will be all-NA.

---

## 9. Dependencies

```r
# Core
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(mgcv)       # GAMMs (bam, gam, vis.gam)
library(ggplot2)

# Classical tests
library(ggpubr)     # stat_compare_means (significance brackets)
library(broom)      # tidy() on aov objects
library(afex)       # aov_ez with Greenhouse-Geisser (optional but recommended)
```

`afex` is optional — `run_classical_tests.R` falls back to base R `aov()` if it is not installed, but without Greenhouse-Geisser sphericity correction. For any publication-quality rmANOVA, `afex` should be installed.

The `%||%` null-coalescing operator is used in `compare_channel_source_gamm.R`. It is available natively in R ≥ 4.4; a fallback definition is included for older installations.

---

## 10. Output Directory Structure

```
R/analysis/behavioural_analysis/
├── behavioural_master.csv
├── behavioural_demo_master.csv
├── alpha_pain_master.csv
├── source_pain_master.csv
├── source_ga_master.csv
├── source_ga_fooof_master.csv
│
├── gamm_outputs_v2/
│   ├── model_comparison_v2.csv
│   ├── trial_level_fitted_values_v2.csv
│   ├── <model_name>.rds
│   ├── <model_name>_summary.txt
│   └── <model_name>_diagnostics.png
│
├── gamm_outputs_source/
│   ├── model_comparison_all_rois.csv
│   ├── best_model_per_roi.csv
│   └── <ROI>/
│       ├── model_comparison_<ROI>.csv
│       ├── model_input_<ROI>.csv
│       ├── trial_level_fitted_<ROI>.csv
│       ├── subject_ga_fitted_<ROI>.csv
│       ├── <model_name>.rds
│       ├── <model_name>_summary.txt
│       ├── <model_name>_diagnostics.png
│       ├── <model_name>_smooth_<term>.png
│       └── obs_vs_fitted_best_<ROI>.png
│
├── gamm_comparison/
│   ├── fit_quality_comparison.csv
│   ├── effect_significance_comparison.csv
│   ├── smooth_concordance.csv
│   ├── fit_quality_plot.png
│   ├── delta_aic_plot.png
│   ├── smooth_overlay_<metric>.png
│   ├── concordance_heatmap.png
│   └── comparison_report.txt
│
└── classical_tests/
    ├── test1_erd_ttest.csv / test1_erd_plot.png
    ├── test2_fit_ttest.csv / test2_fit_plot.png
    ├── test3_rmanova_<metric>.csv + test3_posthoc_<metric>.csv / test3_rmanova_plot.png
    ├── test4_roi_anova_<model>.csv + test4_posthoc_<model>.csv / test4_roi_anova_plot.png
    ├── test5_rayleigh_<col>.csv + test5_rayleigh_summary_<col>.csv / test5_rayleigh_plot_<col>.png
    └── classical_tests_report.txt
```