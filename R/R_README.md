# CNED EEG Analysis Pipeline (R)

**Author:** Dryden Arseneau  
**Affiliation:** Schabrun Lab | Seminowicz Lab  
**Dataset:** CNED - Zhao et al., Sci Data (2025)  
https://doi.ord/10.1038/s41597-025-05900-1

---

## 1. Overview

This R pipeline takes the outputs of the MATLAB preprocessing/spectral pipeline and the Python source localization pipeline, ingests them into a mornalized DuckDB database, and fits a comprehensive set of statistical models. Two parallel analytical streams - **channel-space** (scalp EEG, pooled across electrodes) and **source-space** (sLORETA ROI projections onto fsaverage) - run side by side, enabling a formal comparison of whether source localization improves the prediction of trial-level pain perception.

The primary statistical tool ts the **Generalized Additive Mixed Models (GAMM)** implemented via `mgcv::bam()` with fREML estimation and double-penalty selection (`select = TRUE`). GAMMs are well-suited to trial-level pain data because they: (a) accomodate nonlinear predictor-outcome relationships without imposing a functional form; (b) handle deeply nested random effects (trials within subjects within experiments); and (c) support adaptive smoothness selection that shrinks unsupported terms to zero.

The model-fitting stage is **combo-capable**: a single engine (`gamm_comvo_core.R`) generates and fits all metric combinations of a specified size, from single-metric screens to pairwise or higher-order interaction models, without modifying the fitting or storage logic. All GAMM results are written incrementally to `results.duckdb` as models complete, so a crash mid-run does not lose earlier fits.

Classical inferential tests (`run_classical_tests.R`) and a channel-vs-source comparison (`compare_channel_source_gamm.R`) complement the GAMMs. These scripts currently read the old CSV master files and are flagges for migration to DuckDB in a future pass.

---

## 2. Architecture

### Two-Database Design

| Database | Role | Tables |
|---|---|---|
| `data.duckdb` | Input data - read-only during GAMM-stage | `subjects`, `trials`, `source_metrics`, `spectral_metrics`, `subject_split` |
| `results.duckdb` | Model results - written by GAMM and validate stages | `gamm_fitted`, `power_analysis` |

`data.duckdb` is populated and maintained by the merge stage. `results.duckdb` is written by `source_core.R`, `channel_core.R`, and `power_analysis.R`. The two databases are never written to simultaneously by the same stage.

### Metric Storage: Long/Tidy Format

`source_metrics` and `spectral_metrics` use a long/tidy schema (`global_subjid`, `trial`, `roi`/`channel`, `metric_name`, `metric_value`) rather than wide columns. New metrics from future pipeline revisions appear automatically as new `metric_name` values - no schema migration, no `ALTER TABLE`, no changes to the GAMM script. `db_helpers.R` pivots to wide format transparently before handing data to the fitting loop.

All `metric_name` values are **lowercase** throughout, regardless of source (spectral files already use lowercase; mixed-case source columns such as `BI_pre`, `CoG_pre`, `delta_ERD` are lowercased at ingestion in `read_helpers.R`).

### Train/Holdout Split

`subject_split` records each subject's assignment (`train` or `holdout`), startified by cap size using a quota=based greedy algorithm that converges to exactly the target fraction (default 75/25). Assignments are made incrementally as new subejcts are ingested - a subject's group is set once and never changes by later runs. The GAMM stage reads training subjects only; held-out subjects are reserved for cross-validation (scripts not yet written).

### Pipeline Flow

```
Per-subject MATLAB/Python outputs
(sub-XXX_spectral_chan_by_trial.csv, sub-XXX_source_trial.csv)
Experiment behavioural CSVs, participants.tsv, cap_size.csv
                    │
                    ▼
             R/merge/merge_core.R
                    │
                    ▼
             data.duckdb
      (subjects, trials, source_metrics,
       spectral_metrics, subject_split)
                    │
          ┌─────────┴──────────┐
          ▼                    ▼
  R/gamm/source_core.R   R/gamm/channel_core.R
  (per source-space ROI) (pooled across channels)
          │                    │
          └─────────┬──────────┘
                    ▼
             results.duckdb
             (gamm_fitted)
                    │
                    ▼
    R/validate/power_analysis.R
    (all significant models, per cap_size + pooled)
                    │
                    ▼
             results.duckdb
             (power_analysis)
                    │
          ┌─────────┴──────────┐
          ▼                    ▼
compare_channel_source_gamm.R  run_classical_tests.R
(pending DuckDB migration)     (pending DuckDB migration)
```

---

## 3. File Reference

### `R/merge/`

| File | Role |
|---|---|
| `merge_defaults.R` | All configurable paths and parameters: `data_root`, `behav_dir`, `cap_size_file`, `data_db_path`, experiment lookup tables, split fraction, and seed. Edit here when paths change or a new experiment is added. |
| `merge_core.R` | Main orchestrator. Discovers and reads all per-subject pipline CSVs, assembles `subjects`/`trials`/`source_metrics`/`spectral_metrics` tables in `data.duckdb`, then calls split assignment. Idempotent - safe to re-run as new subjects clear the pipeline. |
| `herlpers/read_helpers.R` | Per-file-type readers: `read_behavioural_file()`, `read_participants_file()`, `read_cap_size_file()`, `read_spectral_file()`, `read_source_file()`. Each returns a tidy data frame in a known shape; metric-name lowercasing happens here. |
| `helpers/db_write_helpers.R` | Schema initialization (`init_data_schema()`) and upsert functions (`upsert_subjects()`, `upsert_trials()`, `upsert_metrics_long()`). Delete-then-insert per subejct - re-ingesting a subject's data on re-run is always safe. |
| `helpers/split_helpers.R` | `assign_new_subjects_to_split()` - quota-based, cap-size-stratified, incremental train/holdout assignment. Existing assignments are never touched. Reports running proportions per cap_size after each call. |

**Replaces:** `merge_behavioural.R`, `merge_participants_into_behavioural.R`, `merge_spectral_behaviour.R`, `merge_source_spectral.R` (trial-level portion). `merge_source_ga.R` (GA master tables) is deferred - not yet ported.

---

### `R/gamm/`

| File | Role |
|---|---|
| `gamm_defaults.R` | GAMM stage configuration: database paths, output directory, QC thresholds (`min_fooof_r2`, `min_trials_grp`), combo sizes (`COMBO_SIZES`), tensor pairs, k ceiling. Flip `COMBO_SIZES` from `1L` to `c(1L, 2L)` here when ready for pairwise interaction screening — no other changes needed. |
| `source_core.R` | Source-space GAMM combo screen. Loops over every ROI present in `source_metrics` (training subjects only). At startup, prints the full train/holdout split breakdown per cap_size as a runtime confirmation that the filter is active. Writes `.rds` per model and one row per fit to `gamm_fitted` incrementally. |
| `channel_core.R` | Channel-space GAMM combo screen. Structurally identical to `source_core.R` but reads `spectral_metrics` and includes `s(channel, bs="re")` in the base formula as a pooled sensor-space random effect — one model per metric-combo across all electrodes simultaneously, not a per-electrode loop. Serves as a secondary validation layer against the source-space primary analysis. |
| `helpers/gamm_combo_core.R` | The shared fitting engine. `generate_metric_combos()` enumerates combinations. `build_combo_formula()` constructs the bam() formula for any metric vector (additive or tensor). `fit_gamm_combo()` fits, saves `.rds`, and returns a result row. `init_results_db()` / `write_result_row()` / `close_results_db()` handle incremental DuckDB writes. |
| `helpers/db_helpers.R` | Database connection and read helpers for the GAMM stage. `connect_data_db()` opens `data.duckdb` read-only (write attempts fail loudly — intentional). `load_source_metrics_wide()` and `load_channel_metrics_wide()` filter to training subjects, pivot long→wide, and join subjects + trials in one query. `load_holdout_source_metrics_wide()` is the held-out counterpart, reserved for CV scripts. `list_available_rois()` / `list_available_channels()` return only groups with training-set coverage. |
 
**Replaces:** `run_gamm_source.R`, `run_gamm_alpha_metrics_v2.R`, `run_gamm_alpha_metrics.R`, `gamm_combo_core.R` (old flat version), `run_gamm_source_screen.R`.

---

### `R/validate/`
 
| File | Role |
|---|---|
| `power_analysis.R` | Validates the train/holdout split across **all** significant single-metric source models in `gamm_fitted` (not just the best one). For each: refits as `lmer()` to extract variance components, runs `simr::powerCurve()` pooled and per cap_size, bootstraps holdout evaluation metrics. Reports the binding constraint (metric × cap_size requiring the most subjects), a pass/fail verdict against current training N, and a suggested `split_train_frac` adjustment if the current split is insufficient. Results written to `power_analysis` table in `results.duckdb`. |
 
**Run frequency:** twice — once on interim data to make the split-fraction decision, once near full N as a confirmatory check. Do not run repeatedly and adjust the split based on results (optional stopping).
 
---
 
### `R/analysis/` (pending DuckDB migration)
 
These scripts remain functional against the old CSV master files. They will be updated to read from `data.duckdb` in a future pass.
 
| File | Status |
|---|---|
| `compare_channel_source_gamm.R` | Reads `gamm_outputs_v2/` and `gamm_outputs_source/` CSV outputs; pending migration to query `gamm_fitted` in `results.duckdb` |
| `run_classical_tests.R` | Reads `source_pain_master.csv` and `alpha_pain_master.csv`; pending migration to read from `data.duckdb` |
 
---

## 4. `gamm_fitted` Schema

All GAMM results land in the `gamm_fitted` table in `results.duckdb`. Primary key: `(level, group_value, model_name)`.

| Column | Type | Description |
|---|---|---|
| `fitted_at` | VARCHAR | Timestamp of this fit |
| `level` | VARCHAR | `'source'` or `'channel'` — which stage fit this model |
| `group_value` | VARCHAR | ROI name (e.g. `postcentral-lh`) for source; `'all_channels'` for channel |
| `model_name` | VARCHAR | Auto-generated from metric combo (e.g. `m_1metric_cog_pre`) |
| `metrics` | VARCHAR | Pipe-delimited metric names (e.g. `cog_pre_z\|bi_pre_z` for a 2-metric combo) |
| `n_metrics` | INTEGER | Number of metrics in this combo (0 = baseline) |
| `edf` | VARCHAR | Pipe-delimited effective degrees of freedom per metric smooth |
| `f_stat` | VARCHAR | Pipe-delimited F statistics per metric smooth |
| `p_value` | VARCHAR | Pipe-delimited p-values per metric smooth |
| `aic` | DOUBLE | Model AIC |
| `bic` | DOUBLE | Model BIC |
| `dev_expl` | DOUBLE | Deviance explained (proportion) |
| `r_sq` | DOUBLE | Adjusted R² |
| `n_obs` | INTEGER | Number of observations used |
| `rds_path` | VARCHAR | Absolute path to the saved `.rds` model object |

**Useful Queries:**

```sql
-- ALL significant source-space single-metric results, ranked by r_sq
SELECT group_value, metrics, r_sq, p_value, aic
FROM gamm_fitted
WHERE level = `source` AND n_metrics = 1
    AND CAST(p_value AS DOUBLE) < 0.05
ORDER BY r_sq DESC:

-- Compare baseline vs single-metric AIC per ROI
SELECT as.group_value AS roi, b.metrics, b.aic - a.aic AS delta_aic
FROM gamm_fitted a
JOIN gamm_fitted b ON a.group_value = b.group_value AND a.level = b.level
WHERE a.model_name = 'm00_baseline' AND b.n_metrics = 1
    AND a.level = 'source'
ORDER BY delta_aic;

-- SHow all 2-metric combos that beat the single-metric baseline AIC
SELECT group_value, metrics, aic, r_sq
FROM gamm_fitted
WHERE level = 'source' AND n_metrics = 2 AND aic < (
    SELECT MIN(aic) FROM gamm_fitted gf2
    WHERE gf2.level = 'source' AND gf2.n_metrics = 1
    AND gf2.group_value = gamm_fitted.group_value
)
ORDER BY aic;
```

---

## 5. GAMM Model Structure

All models share the same base formula. Continuous predictors are z-scored. `select = TRUE` adds a shrinkage penalty to every smooth, allowing unsupported terms to be shrunk to zero.

**Base Formula (all models):**
```
pain_rating ~ s(laser_power_z, k=k_laser) + s(trial_index_z, k=k_trial)
            + age_z + [sex] + [cap_size]
            + s(global_subjid, bs="re") + [s(experiment_id, bs="re")]
```

Fixed effects (`age`, `cap_size`) and random effects (`experiment_id`) are included conditionally - only when the factor has more than one level in the current ROI/channel subset after QC filtering, to avoid "contrasts not defined" errors on small or homogenous subsets. `age_z` is always linear (not a smooth) - appropriate at current interim N and consistent across channel and source GAMM families.

**Per-metric smooth term (added to base):**
```
+ s(metric_z, k = safe_k(data[[metric]], k_max = 10))
```

**Tensor interaction term (for designated pairs, e.g., pow_slow * pow_fast):**
```
+ te(metric_a_z, metric_b_z, k = c(k_a, k_b))
```

`safe_k()` caps smooth complexity at `min(k_max, n_unique_values - 1)` with a floor of 3, preventing over-parameterization when a metric has a few unique values at small N.

---

## 6. Multiple Comparisons Stance

The GAMM screening is explicitly **exploratory**. No correction is applied at the screening stage. The train/holdout split serves as the de facto correction: a metric x ROI that reaches significance in screening must replicate in the held-out sample (via the cross-validation scripts, not yet written) before being treated as a confirmed finding.

This decision is recorded in `gamm_defaults.R` as `apply_fdr_correction = FALSE` so it is visible in configuration, not buried in code.

---

## 7. Run Sequence
Scripts must be run from their stage directory (`R/merge/` or `R/gamm/`) so relative `source()` calls resolve correctly. Steps 1–2 are run once to initialize the database; step 3 is re-run each time new subjects clear the pipeline.
 
```
Step 1 (once):   R/merge/merge_core.R
                   → populates data.duckdb
                   → assigns initial train/holdout split
 
Step 2 (re-run): R/merge/merge_core.R
                   → re-ingests updated subjects (upsert, safe to repeat)
                   → assigns split for any new subjects only
 
Step 3:          R/gamm/source_core.R
                   → screens all metrics × all source-space ROIs
                   → writes to results.duckdb (gamm_fitted)
 
Step 4:          R/gamm/channel_core.R
                   → screens all metrics, pooled channel-space
                   → writes to results.duckdb (gamm_fitted)
 
Step 5 (once):   R/validate/power_analysis.R
                   → runs on all significant gamm_fitted hits
                   → reports binding constraint + split adequacy verdict
                   → writes to results.duckdb (power_analysis)
                   → [adjust split_train_frac in merge_defaults.R if needed,
                      then re-run merge_core.R split step only]
 
Step 6 (once):   R/validate/power_analysis.R  [confirmatory, near full N]
 
--- pending ---
 
Step 7:          Cross-validation scripts (not yet written)
                   → LOSO internal validation on training set
                   → Stratified holdout external validation
                   → Reads held-out subjects via
                     load_holdout_source_metrics_wide()
 
Step 8:          compare_channel_source_gamm.R  (pending DuckDB migration)
Step 9:          run_classical_tests.R           (pending DuckDB migration)
```
 
---
 
## 8. Dependencies
 
```r
# Core
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(mgcv)       # GAMMs — bam(), select=TRUE, fREML
library(ggplot2)
 
# Database
library(DBI)
library(duckdb)
 
# Power analysis
library(lme4)       # lmer() refit for simr
library(simr)       # powerCurve(), powerSim()
 
# Classical tests (R/analysis/ scripts, pending migration)
library(ggpubr)
library(broom)
library(afex)       # optional; falls back to base aov() if absent
```
 
`afex` is optional for `run_classical_tests.R` — without it, repeated-measures ANOVA runs without Greenhouse-Geisser sphericity correction. For publication-quality output, install `afex`.
 
`simr` requires `lme4`. The power analysis script refits `bam()` models as `lmer()` for simulation purposes only — see the `power_analysis.R` header for the methodological rationale.
 
---
 
## 9. Output Directory Structure
 
```
db/
├── data.duckdb          ← subjects, trials, source_metrics,
│                           spectral_metrics, subject_split
└── results.duckdb       ← gamm_fitted, power_analysis
 
R/gamm_outputs/
├── <roi>/               ← one folder per source-space ROI
│   ├── model_input_<roi>.csv
│   └── <model_name>.rds
├── all_channels/
│   └── <model_name>.rds
└── power_analysis/
    └── power_curve_<metric>_<roi>.png
 
R/analysis/              ← legacy; pending DuckDB migration
├── behavioural_analysis/
│   └── (old merge + GAMM scripts, still functional against CSV masters)
└── experiment/
    └── (per-experiment behavioural CSVs)
```
 
---
 
## 10. Key Design Decisions (and why)
 
**Long/tidy metric tables over wide columns.** New metrics appear as new rows — no schema change, no `ALTER TABLE`, no downstream script edits. The pivot to wide happens once in `db_helpers.R` transparently.
 
**Lowercase metric names throughout.** Spectral pipeline files already use lowercase (`bi_pre`, `cog_pre`); source pipeline files use mixed case (`BI_pre`, `CoG_pre`). Lowercasing at ingestion in `read_helpers.R` ensures `bi_pre` and `BI_pre` become the same `metric_name` value, not two separate metrics that silently fragment the long-format table.
 
**`global_subjid` = `subjid_uid` string (e.g. `E02_S019`), not a row-number integer.** Row-number IDs shift across re-runs when new subjects are added, silently desyncing every downstream table that references them. A deterministic string derived from `(experiment_id, subjid)` is stable across all re-runs.
 
**Quota-based stratified split, not independent random draws.** For the smallest cap_size group (n=124 32-channel subjects), independent draws at p=0.75 have a 5th–95th percentile range of roughly 68–81% — a ±5–10pp deviation that meaningfully affects the holdout group's stability. Quota-based assignment converges to exactly 75/25 per cap_size group regardless of N.
 
**Incremental split assignment, never reassigned.** A subject added today gets a `split_group` immediately, so screening can proceed on interim data while the dataset grows. Because assignments are never overwritten, the split is reproducible even as new subjects arrive in later processing batches.
 
**`data.duckdb` opened read-only during the GAMM stage.** `connect_data_db()` enforces `read_only = TRUE` — an accidental write from the GAMM stage fails loudly rather than silently corrupting the upstream data that the merge stage owns.
 
**Multiple comparisons: no correction at the screening stage.** The train/holdout split IS the correction — a finding must replicate in held-out data. This is methodologically preferable to correcting on a non-independent pilot set and is explicitly documented in `gamm_defaults.R` so the decision is traceable.