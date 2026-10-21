# CNED EEG Analysis Pipeline (R)

**Author:** Dryden Arseneau  
**Affiliation:** Schabrun Lab | Seminowicz Lab  
**Dataset:** CNED - Zhao et al., Sci Data (2025)  
https://doi.org/10.1038/s41597-025-05900-1

**Feature set:** v3.1 (spectral_core.m 3.1.0, source_core.py 4.1.0, preproc_core.m 2.2.1)  
**R stage:** V3 (Oct 2026): non-discrete fits by default, decomposition-keyed results, robustness suite, channel-vs-source comparison and classical tests on DuckDB

---

## 1. Overview

This R pipeline takes the outputs of the MATLAB preprocessing/spectral pipeline and the Python source-localization pipeline, checks that every subject's EEG trials line up with their behavioural trials, loads everything into a normalized DuckDB database, fits a screen of statistical models, and then tests every hit: held-out replication, sex differences, power, robustness, and channel vs source.

Two parallel analytical streams run side by side:

- **Channel space** — scalp electroencephalography (EEG) features, one model set per electrode.
- **Source space** — standardized low-resolution electromagnetic tomography (sLORETA) projections onto fsaverage, one model set per region of interest (ROI).

Together they support a formal comparison of whether source localization improves the prediction of trial-level pain ratings (`compare_levels.R`).

The primary statistical tool is the **generalized additive mixed model (GAMM)**, fitted with `mgcv::bam()` using fast restricted maximum likelihood (fREML) and double-penalty selection (`select = TRUE`). GAMMs suit trial-level pain data because they:

- accommodate nonlinear predictor–outcome relationships without imposing a functional form;
- handle nested random effects (trials within subjects within experiments);
- shrink unsupported smooth terms towards zero.

The model-fitting stage is **combo-capable**. One engine (`gamm_combo_core.R`) generates and fits every metric combination of a chosen size, from single-metric screens to pairwise or higher-order models, without changes to the fitting or storage logic. Results are written to `results.duckdb` one model at a time, so a crash mid-run loses at most the model in flight.

---

## 2. Architecture

### Two-Database Design

| Database | Role | Tables |
|---|---|---|
| `data.duckdb` | Input data; read-only during the GAMM, validate and robustness stages | `subjects`, `trials`, `epochs`, `alignment_qc`, `spectral_metrics`, `source_metrics`, `spectral_subject_metrics`, `source_subject_metrics`, `ingest_log`, `subject_split` |
| `results.duckdb` | Model results; written by the GAMM, validate and robustness stages | `gamm_fitted`, `gamm_sex_diff`, `holdout_replication`, `power_analysis`, `robustness` |

The merge stage builds `data.duckdb`. `results.duckdb` is written by the scripts in `R/gamm/`, `R/validate/` and `R/robustness/`. No stage writes to both databases.

### `data.duckdb` Tables

| Table | Grain | Contents |
|---|---|---|
| `subjects` | subject | `global_subjid`, `subjid`, `experiment_id`, `experiment_name`, `age`, `sex`, `cap_size` |
| `trials` | subject × trial | `trial`, `trial_index`, `laser_power`, `pain_rating` (from behaviour) |
| `epochs` | subject × epoch | `trial` (= epoch index), `event_code`, `condition`, `intensity` (low/high/unknown) from `07_epoch/sub-XXX_epoch_metadata.csv` |
| `alignment_qc` | subject | trial-alignment check result: counts, intensity check, `status`, `reason`, `ingested` (Section 2, Trial Alignment) |
| `spectral_metrics` | subject × trial × channel × metric | trial-level channel-space features (long) |
| `source_metrics` | subject × trial × ROI × metric | trial-level source-space features (long) |
| `spectral_subject_metrics` | subject × channel × metric | subject-level summaries: trial-variability index (TVI), inter-trial coherence (ITC), fitting oscillations and one-over-f (FOOOF) parameters, trial-mean window features |
| `source_subject_metrics` | subject × ROI × metric | the same summaries per ROI |
| `ingest_log` | subject × file | which file produced how many rows, under which feature set (rebuilt each run; append-only history if `rebuild_derived_tables = FALSE`) |
| `subject_split` | subject | `train` / `holdout` assignment (persistent; Section 2, Train/Holdout Split) |

### Metric Storage: Long/Tidy Format

All four metric tables use a long schema (`global_subjid`, [`trial`], `roi`/`channel`, `metric_name`, `metric_value`) rather than wide columns. A new feature from a later pipeline revision appears as a new `metric_name` value, with no schema migration. `db_helpers.R` pivots to wide format before handing data to the fitting loop.

Rules applied at ingestion (`ingest_feature_csv()` in `db_write_helpers.R`):

- every `metric_name` is **lowercased**;
- `NaN`, `Inf` and empty values are **not stored** (a missing row means "no value");
- the per-subject comma-separated values (CSV) file is read and unpivoted inside DuckDB in one SQL statement, so the long tables are never held in R memory (`spectral_metrics` alone is of the order of 10⁸ rows at full N: subjects × trials × 64 channels × 58 metrics).

With feature set v3.1, the trial-level tables hold the same 58 metrics per unit in both spaces (spectral and source parity is tested in `tests/parity_test.py`):

| Group | Metrics |
|---|---|
| Alpha families (10) × windows `whole_`, `pre_`, `post_`, `delta_` | `pow_slow_alpha`, `pow_fast_alpha`, `pow_alpha_total`, `paf_cog_hz` (peak alpha frequency, centre of gravity), `rel_slow_alpha`, `rel_fast_alpha`, `sf_ratio`, `sf_logratio`, `sf_balance`, `psi_cog` |
| Event-related desynchronization (ERD) | `erd_slow`, `erd_fast`, `erd_pow_alpha_total`, `erd_asym` |
| Quality control (QC) flag | `p5_flag` |
| Phase | `slow_phase`, `sin_phase`, `cos_phase` and their `_post` versions |
| Laser-evoked potential (LEP) | `n2_amp`, `n2_lat_ms`, `p2_amp`, `p2_lat_ms`, `n2p2_amp`, `n2_mean`, `p2_mean` |

### Trial Alignment (merge stage)

Behaviour is joined to EEG by `trial` = epoch index. That is only correct if the k-th epoch is the k-th behavioural trial; one extra or missing epoch shifts every later rating onto the wrong trial with no error. Before any EEG row is written, `assess_alignment()` (`helpers/alignment_helpers.R`) checks each subject:

1. Behavioural trials must be numbered 1..n.
2. The EEG trial count must equal the behavioural trial count. The reference count is the epoch metadata if present, otherwise the spectral file, otherwise the source file; the spectral and source files must also agree with it.
3. For experiments whose epochs carry low/high labels (two laser levels), every `high` epoch must have the subject's higher `laser_power` and every `low` epoch the lower one. Laser order is randomized, so a one-trial shift breaks about half the trials and is caught.

| `status` | Meaning | EEG rows written? |
|---|---|---|
| `ok` | all checks passed (or the intensity check did not apply) | yes |
| `count_mismatch` | epoch / spectral / source / behaviour counts differ | no (`exclude_misaligned = TRUE`) |
| `intensity_mismatch` | low/high labels disagree with `laser_power` | no (`exclude_misaligned = TRUE`) |
| `behaviour_noncontiguous` | behavioural trial numbers are not 1..n | no |
| `rating_out_of_range` | a rating outside `rating_range` (data-entry error) | no |
| `no_eeg` | behaviour only (not yet processed) | — |
| `no_behaviour` | EEG outputs but no behavioural rows | no |
| `read_error` / `ingest_error` | a file could not be read / inserted | no (transaction rolled back) |

`intensity_check` records whether step 3 ran: `pass`, `fail`, `labels_unknown` (epoch metadata written before preproc_core.m 2.2.0 has no low/high labels; re-running preprocessing fixes this), `single_level`, `no_metadata`, `not_run`. `alignment_qc` doubles as the missing-data audit: every subject with behaviour or EEG appears there with a reason.

### Train/Holdout Split

`subject_split` records each subject's assignment (`train` or `holdout`), stratified by cap size using a quota-based greedy algorithm that converges to the target fraction (default 75/25).

- Assignments are made incrementally as new subjects are ingested. A subject's group is set once and never changed by later runs.
- Every subject with behaviour is assigned, including subjects excluded by the alignment check. The split belongs to the person, so a subject who is fixed and re-ingested later keeps their original group.
- `subject_split` is the only table the merge rebuild never drops.
- The GAMM screen reads training subjects only; held-out subjects are reserved for `holdout_replicate.R` and the held-out parts of the robustness tests.

### Pipeline Flow

```
Per-subject MATLAB / Python outputs
  preproc/…/sub-XXX_epoch_metadata.csv
  spec/…/sub-XXX_spectral_chan_by_trial.csv, sub-XXX_spectral_chan_summary.csv
  source/…/sub-XXX_source_trial.csv,          sub-XXX_source_ga.csv
Experiment behavioural CSVs, participants.tsv, cap_size.csv
                    │
                    ▼
             R/merge/merge_core.R
        (alignment check → ingest → split)
                    │
                    ▼
               data.duckdb
                    │
          ┌─────────┴──────────┐
          ▼                    ▼
  R/gamm/source_core.R   R/gamm/channel_core.R
  (per ROI)              (per channel)
          └─────────┬──────────┘
                    ▼
       results.duckdb (gamm_fitted)
                    │
   ┌────────────┬───┴────────┬──────────────┬─────────────────┐
   ▼            ▼            ▼              ▼                 ▼
sex_diff_   plot_gamm_   holdout_       compare_levels.R   ../validate/
core.R      hits.R       replicate.R    (channel vs        power_analysis.R
                             │           source)
                             ▼
                 ../robustness/robustness_core.R
                 (replicated hits: LOEO, k-fold,
                  permutation, split-half, FOOOF)

classical_tests.R reads data.duckdb only (ERD t-tests, rmANOVA, Rayleigh).
```

---

## 3. File Reference

### `R/merge/`

| File | Role |
|---|---|
| `merge_defaults.R` | Paths (`data_root`, `behav_dir`, `cap_size_file`, `data_db_path`), `feature_set`, `rebuild_derived_tables`, the per-experiment file layout (`files`), `exclude_misaligned`, `rating_range`, the experiment lookup, split seed and fraction. Edit here when paths change or an experiment is added. |
| `merge_core.R` | Orchestrator. (1) drops and recreates the derived tables; (2) reads behaviour, participants and cap size; (3) discovers each experiment's pipeline outputs; (4) per subject: alignment check, then all of that subject's EEG files ingested in one transaction; (5) writes `alignment_qc` and prints a status summary; (6) assigns the split to new subjects. Stops if EEG files exist but no subject was ingested. |
| `helpers/read_helpers.R` | `read_behavioural_file()`, `read_participants_file()`, `read_cap_size_file()`, `read_epoch_metadata()`, and `discover_pipeline_files()` (subject ID from the file name; errors if one subject has two files of the same kind). |
| `helpers/db_write_helpers.R` | Schema (`init_data_schema()`), rebuild (`drop_derived_tables()`), delete-then-insert helpers, `ingest_feature_csv()` (SQL unpivot of one CSV) and `log_ingest()`. |
| `helpers/alignment_helpers.R` | `assess_alignment()` (see Trial Alignment). |
| `helpers/split_helpers.R` | `assign_new_subjects_to_split()`: quota-based, cap-size-stratified, incremental. Stops if a subject has no cap size. |

### `R/gamm/`

All scripts are run from `R/gamm/`.

| File | Role |
|---|---|
| `gamm_defaults.R` | Database paths, output directory, `feature_set`, metric selection (`metric_include` / `metric_exclude`), `metric_decomposition` / `primary_part`, unit selection (`roi_include` / `channel_include`), QC (`exclude_p5_flag`, `min_fooof_r2`, `min_trials_grp`), `combo_sizes`, `tensor_pairs`, `k_max`, `nthreads`, `bam_discrete`, follow-up and held-out settings, figure settings, `robustness`, `compare` and `classical` blocks, `apply_fdr_correction`. |
| `source_core.R` | Entry point for the source-space screen: `run_gamm_screen("source")`. |
| `channel_core.R` | Entry point for the channel-space screen: `run_gamm_screen("channel")`. One model set per channel; channels are not pooled into one model. |
| `helpers/gamm_screen.R` | The per-unit loop shared by both entry points: load the screened metrics for training subjects, QC, z-score, build the base formula, fit every combo, write each result. |
| `helpers/gamm_combo_core.R` | Fitting engine: `generate_metric_combos()`, `build_combo_formula()` (additive or tensor), `fit_gamm_combo()` (fit, save `.rds`, return a result row), `model_out_dir()`, `init_results_db()` / `write_result_row()` / `close_results_db()`. |
| `sex_diff_core.R` | Follow-up on every significant single-metric model: refits it on the same rows with sex as an ordered factor, so a **difference curve** (other − reference sex) tests whether the curve shape differs. Writes `gamm_sex_diff` (Section 4). |
| `holdout_replicate.R` | **The replication step.** For every significant single-metric model: the training-fit metric model and a nested baseline (same training rows, metric term removed) both predict the held-out subjects at population level. Change in prediction error at trial and subject level, 95% CI from resampling held-out subjects, Benjamini–Hochberg (BH) across hits. Writes `holdout_replication` and a forest plot per level. |
| `plot_gamm_hits.R` | Figures for every significant single-metric model (Section 5, Figures). Reads only `results.duckdb` and the `.rds` files. |
| `compare_levels.R` | **Channel vs source.** Per metric: share of units with a screen hit at each level, held-out replication counts; for metrics with a hit at both levels, a paired held-out comparison of the best ROI and the best channel (chosen by screen p, training data only) and the correlation of their fitted curves. Writes `compare/level_summary.csv`, `compare/level_paired_holdout.csv`, `figures/level_comparison.png`. Replaces `R/analysis/behavioural_analysis/compare_channel_source_gamm.R`. |
| `classical_tests.R` | Classical tests on the v3.1 source tables: ERD (paired t-tests, post vs pre power per ROI × band), metric differences across ROIs (repeated-measures ANOVA, Greenhouse–Geisser via `afex` if installed, BH pairwise), Rayleigh tests of phase uniformity per subject × ROI. Writes CSVs, figures and a text report to `classical$out_dir`. Replaces `R/analysis/behavioural_analysis/run_classical_tests.R`. |
| `helpers/hit_helpers.R` | Shared by the follow-up scripts: `load_hits()`, `ensure_derived_table()`, `fig_unit_dir()`, `metric_label()` / `unit_label()`, population-level predictions (`marginal_lpmatrix()`, `predict_from_lp()`), held-out prediction (`predict_population()`, `prepare_holdout()`, `refit_baseline()`, `boot_subjects()`), the figure theme, `fmt_p()`. |
| `helpers/db_helpers.R` | `connect_data_db()` (read-only), `list_available_units()` / `list_available_metrics()`, `select_metrics()` / `select_units()`, `load_metrics_wide()` (trial level; wrappers `load_source_metrics_wide()`, `load_channel_metrics_wide()`, `load_holdout_source_metrics_wide()`), `load_subject_metrics_wide()` (subject level). |
| `check_discrete_between_p.R` | 30-second check of whether `bam(discrete = TRUE)` gives valid p-values for subject-level smooths on your mgcv version (Section 5, `bam_discrete`). |

**Metrics screened.** `metric_include` selects the 10 alpha families × 4 windows plus the 4 ERD metrics: **44 metrics** per unit. Not screened:

- `p5_flag` — a QC flag. With `exclude_p5_flag = TRUE` it drops trials instead.
- Phase columns — circular, so a smooth on a linear axis is wrong (tested for uniformity in `classical_tests.R` instead).
- LEP columns — a different signal and question.

**ROIs.** `source_metrics` holds 18 ROIs: S1, S2, ACC (anterior cingulate cortex), Insula, dlPFC (dorsolateral prefrontal cortex) and M1, each bilateral plus `_lh` and `_rh`. Screening all three versions tests each region three times with partly redundant data. Restrict with `roi_include` if that is not wanted.

### `R/validate/`

| File | Role |
|---|---|
| `power_analysis.R` | For every significant single-metric source model of the current `feature_set` and decomposition: refits it as `lmer()` on the screen model's own rows (same QC and scaling; both parts of a split metric, linear; sex and cap size as in the screen), runs `simr::powerCurve()` for the **primary part** pooled and per cap size (after `simr::extend()` to reach the projected N), at the observed effect **and at half of it** (`effect_scales`: screen hits are selected for looking large, so their estimates are inflated), and bootstraps holdout prediction error. Per-cap-size requirements are converted to total training N. Reports the binding constraint and a verdict at the smallest effect scale. Writes `power_analysis` in `results.duckdb`. Run from `R/gamm/`. |

**Run frequency:** twice — once on interim data to make the split-fraction decision, once near full N as a confirmatory check. Do not run repeatedly and adjust the split based on results (optional stopping).

### `R/robustness/`

| File | Role |
|---|---|
| `robustness_core.R` | Robustness tests for hits (by default only those that **replicated** in `holdout_replication`): leave-one-experiment-out (LOEO), grouped k-fold cross-validation (CV) within the training set, subject-level permutation test on the held-out set, split-half stability of the fitted curve, FOOOF-fit sensitivity. Writes `robustness` (Section 4) and `<metric>_robustness.png` per hit. Run from `R/gamm/`: `source("../robustness/robustness_core.R")`. |
| `helpers/robustness_helpers.R` | `refit_formula()` / `refit_pair()` (screen model and nested baseline refit on a subset, terms rebuilt for that subset), `score_predictions()`, `permute_primary()`, `stratified_folds()`, the `robustness` table. |

The drafts `run_robustness.R`, `run_cv_internal.R`, `run_cv_external.R`, `cv_helpers.R` and `cv_split.R` are replaced by these two files (external validation is `holdout_replicate.R`; leave-one-subject-out is replaced by grouped k-fold, which needs k instead of N refits per hit).

### Not yet migrated

| Location | Status |
|---|---|
| `R/analysis/behavioural_analysis/` (old merge / GAMM scripts, `gamm_combo_core.R` copy) | Legacy CSV-based scripts, superseded by `R/merge/` and `R/gamm/`. Kept for reference only; do not source them alongside the current stage. |

---

## 4. `results.duckdb` Schemas

### `gamm_fitted`

Primary key: `(feature_set, decomposition, level, group_value, model_name)`. `init_results_db()` refuses to open a `gamm_fitted` table without the `feature_set`, `decomposition` and `bam_discrete` columns and asks for `results.duckdb` to be archived (those rows came from discrete fits, whose between-subject p-values are invalid).

| Column | Type | Description |
|---|---|---|
| `fitted_at` | VARCHAR | Timestamp of this fit |
| `feature_set` | VARCHAR | Pipeline feature set the metrics came from (e.g. `v3.1`) |
| `level` | VARCHAR | `'source'` or `'channel'` |
| `group_value` | VARCHAR | ROI (e.g. `Insula_lh`) or channel label (e.g. `Cz`) |
| `model_name` | VARCHAR | `m00_baseline` or generated from the combo (e.g. `m_1metric_pre_paf_cog_hz`) |
| `metrics` | VARCHAR | Pipe-delimited z-scored metric names (e.g. `pre_paf_cog_hz_z\|erd_asym_z`) |
| `n_metrics` | INTEGER | Metrics in this combo (0 = baseline) |
| `decomposition` | VARCHAR | `between_within` or `pooled` |
| `primary_part` | VARCHAR | `between`, `within` or `pooled`: which term `edf` / `f_stat` / `p_value` describe |
| `edf` | VARCHAR | Pipe-delimited effective degrees of freedom per metric smooth |
| `f_stat` | VARCHAR | Pipe-delimited F statistics per metric smooth |
| `p_value` | VARCHAR | Pipe-delimited p-values per metric smooth (primary part) |
| `edf_other`, `f_other`, `p_other` | VARCHAR | The same for the other part of a split metric |
| `aic` | DOUBLE | Akaike information criterion |
| `bic` | DOUBLE | Bayesian information criterion |
| `dev_expl` | DOUBLE | Deviance explained (proportion) |
| `r_sq` | DOUBLE | Adjusted R² |
| `n_obs` | INTEGER | Observations used |
| `n_subjects` | INTEGER | Subjects contributing at least one complete row |
| `bam_discrete` | BOOLEAN | Whether the fit used `discrete = TRUE` (p-values of subject-level terms invalid if TRUE) |
| `rds_path` | VARCHAR | Path to the saved model object |

Smooth-term statistics are matched to metrics by the variables inside each term label (`s(x)`, `te(a,b)`), not by substring. In a tensor pair both metrics report the shared `te()` term.

**Useful queries:**

```sql
-- Significant source-space single-metric results, ranked by p
SELECT group_value, metrics, p_value, p_other, r_sq
FROM gamm_fitted
WHERE feature_set = 'v3.1' AND decomposition = 'between_within'
  AND level = 'source' AND n_metrics = 1
  AND TRY_CAST(p_value AS DOUBLE) < 0.05
ORDER BY TRY_CAST(p_value AS DOUBLE);

-- Hits that replicated, with their robustness headline numbers
SELECT h.level, h.group_value, h.metric, h.subj_pct, r.test, r.variant, r.stat, r.p_value
FROM holdout_replication h
JOIN robustness r USING (feature_set, decomposition, level, group_value, metric)
WHERE h.replicated AND r.variant IN ('all', 'subject_level')
ORDER BY h.level, h.group_value, h.metric, r.test;
```

**Useful `data.duckdb` queries:**

```sql
-- Who was excluded, and why
SELECT experiment_name, status, COUNT(*) AS n
FROM alignment_qc GROUP BY ALL ORDER BY ALL;

SELECT global_subjid, status, reason
FROM alignment_qc WHERE NOT ingested AND status <> 'no_eeg';
```

### `gamm_sex_diff`

One row per significant single-metric model (primary key `(feature_set, decomposition, level, group_value, metric)`). Derived from `gamm_fitted` and the `.rds` files; an older layout is dropped and rebuilt.

| Column | Description |
|---|---|
| `status` | `ok`, or why the refit was skipped: `no_sex_term`, `not_two_sexes`, `too_few_subjects` (fewer than `sex_min_subjects` per sex), `fit_failed`, `rds_missing` |
| `term_var` | the term compared between sexes (the primary part, e.g. `x_between_z`) |
| `sex_reference`, `sex_other` | e.g. `F`, `M`; the difference curve is other − reference |
| `n_subj_ref`, `n_subj_other` | subjects of each sex in the model |
| `term_type` | `smooth` (difference smooth) or `linear` (metric entered linearly in the screen; slope difference instead) |
| `diff_edf`, `diff_f`, `diff_p` | statistics of the difference curve: does the shape differ by sex? |
| `diff_p_fdr` | BH-adjusted `diff_p` across all hits of the run |
| `sex_offset`, `sex_offset_p` | overall difference in level (other − reference) |
| `aic_screen`, `aic_sex`, `delta_aic` | AIC of the screen model, of the sex model, and the change (negative = sex model better; same rows, so comparable) |
| `rds_path` | `<out_dir>/<level>/<unit>/[pooled/]sexdiff_<metric>.rds` |

### `holdout_replication`

One row per significant single-metric model (primary key `(feature_set, decomposition, level, group_value, metric)`). Derived; an older layout is rebuilt.

| Column | Description |
|---|---|
| `status` | `ok`, or why the hit was not evaluated |
| `n_holdout_subj`, `n_holdout_trials` | held-out data used |
| `trial_mse_base`, `trial_mse_metric`, `trial_dmse` (+ `_lo`, `_hi`), `trial_pct`, `trial_p` | trial-level prediction error: baseline, metric model, change (negative = metric better) with 95% CI, % change, one-sided bootstrap p |
| `subj_*` | the same for each subject's mean rating (who is more pain-sensitive), plus `subj_r_base` / `subj_r_metric` (correlation of observed and predicted subject means) |
| `primary_level` | `subject` or `trial` (`holdout_primary`) |
| `primary_p_fdr` | BH-adjusted p of the primary level across all hits |
| `replicated` | primary error at least `holdout_min_pct` % lower than baseline **and** `primary_p_fdr` < `holdout_alpha` |

### `power_analysis`

One row per metric × ROI × group (`pooled`, `cap_size_<n>`) × effect scale × power target, plus one `holdout_bootstrap` row per metric × ROI. An older layout is renamed to `power_analysis_pre_v3`, not dropped.

| Column | Description |
|---|---|
| `feature_set`, `decomposition`, `metric`, `roi`, `group_label` | identity |
| `effect_scale` | multiplier on the metric's estimated coefficient (1 = observed) |
| `power_target`, `achieved`, `max_power` | target, whether the curve reached it within the projected N, highest power on the curve |
| `n_current` | subjects in the group now |
| `n_required` | subjects in the group needed for the target |
| `n_required_total` | the same as total training N (`n_required` / the group's share of training subjects) |
| `holdout_r`, `holdout_r_lo`, `holdout_r_hi`, `holdout_rmse`, `holdout_stable` | holdout bootstrap rows only |

### `robustness`

One row per hit × test × variant (primary key `(feature_set, decomposition, level, group_value, metric, test, variant)`). Derived; an older layout is rebuilt.

| `test` | `variant` | `stat` / `stat_lo` / `stat_hi` | `p_value` |
|---|---|---|---|
| `loeo` | `exp_<id>`, `all` | primary change in error (%) for subjects of the left-out experiment, 95% CI | one-sided bootstrap p |
| `kfold` | `fold_01`…, `all` | the same for left-out training subjects | one-sided bootstrap p |
| `permutation` | `subject_level` | observed primary change in error (%); 2.5 / 97.5% of the null | (1 + #null ≤ observed) / (n_perm + 1) |
| `split_half` | `all` | median correlation of the two half-sample metric curves; 2.5 / 97.5% over repetitions | — |
| `fooof_sensitivity` | `fooof_r2_<thr>` | primary change in error (%) after dropping subjects with `fooof_r2` < thr | one-sided bootstrap p |

Prediction tests also fill `trial_*` and `subj_*` (as in `holdout_replication`, plus `_pct_lo` / `_pct_hi`). `status` is `ok` or the reason a variant was skipped; `note` explains special cases (e.g. a cap size absent from a LOEO refit).

---

## 5. GAMM Model Structure

All models share the same base formula. Continuous predictors are z-scored across all training rows for the unit. `select = TRUE` adds a shrinkage penalty to every smooth so unsupported terms can be shrunk to zero.

**Base formula (all models):**
```
pain_rating ~ s(laser_power_z, k = k_laser) + s(trial_index_z, k = k_trial)
            + age_z + [sex] + [cap_size]
            + s(global_subjid, bs = "re") + [s(experiment_id, bs = "re")]
```

Bracketed terms are included only when the factor has more than one level in the current unit after QC, to avoid "contrasts not defined" errors on small or homogeneous subsets. `age_z` is always linear (a between-subject covariate).

**Per-metric smooth term (added to base):**
```
+ s(metric_z, k = safe_k(data[[metric]], k_max = 10))
```

**Tensor interaction term (pairs listed in `tensor_pairs`, e.g. `pre_pow_slow_alpha_z` × `pre_pow_fast_alpha_z`):**
```
+ te(metric_a_z, metric_b_z, k = c(k_a, k_b))
```

`safe_k()` caps smooth complexity at `min(k_max, n_unique_values - 1)` with a floor of 3. A covariate or metric with fewer than 3 unique values (e.g. `laser_power` in a unit where only two-level experiments remain, or a 0/1 metric) enters **linearly** instead; its row in `gamm_fitted` then reports edf = 1, F = t² and the parametric p-value.

**`bam_discrete` (default `FALSE`).** `discrete = TRUE` is several times faster, but on mgcv 1.9-1 its smooth-term p-values are invalid for terms that are constant within a subject: the `_between_z` part and the subject random effect. In synthetic tests a true between-subject effect gave p = 0.999 with `discrete = TRUE` and p = 7 × 10⁻⁷ without it (same edf, same F), and null between-subject terms also gave p ≈ 1, at 66, 135 and 360 subjects (3–9 experiments). Estimates and predictions agree; only the tests differ. Because the between-subject p-value defines a hit, discrete fits make the screen return no valid hits. Run `check_discrete_between_p.R` before considering `TRUE` on another mgcv version; `run_gamm_screen()` warns if `TRUE` is combined with the between/within split, and every `gamm_fitted` row records the setting.

**Sex-difference model (`sex_diff_core.R`).** The screen model with `sex` replaced by an ordered factor `sex_o` and one term added:
```
... + sex_o + s(metric_z) + s(metric_z, by = sex_o)
```
With an ordered factor, `s(metric_z)` is the reference sex's curve and `s(metric_z, by = sex_o)` is the difference curve, so its p-value tests directly whether the shapes differ. (An unordered `by = sex` would give two curves each tested against zero, which does not answer that question.) The refit uses the screen model's own data frame, so both models see exactly the same trials.

### Figures (`plot_gamm_hits.R`)

| File | Content |
|---|---|
| `<fig_dir>/<level>/<unit>/<metric>_effect.png` | Predicted pain rating across the metric, in its original units; 95% confidence band; rug of the data; edf, F, p and N. With the split: A = subject-mean curve (the tested part, blue), B = within-person curve (grey), same y-scale |
| `<fig_dir>/<level>/<unit>/<metric>_sex.png` | A: predicted curve per sex (direct-labelled with n). B: difference curve (Male − Female) with 95% band; x-ranges where the band excludes 0 are shaded. Written only for hits with an `ok` row in `gamm_sex_diff`. |
| `<fig_dir>/<level>/<unit>/<metric>_robustness.png` | Primary change in held-out error per LOEO experiment, k-fold, FOOOF variant (from `robustness_core.R`); permutation p and split-half r in the subtitle |

Pooled-decomposition figures go to a `pooled/` sub-folder of each unit.

"Predicted" means:
- the y-axis is in pain-rating units, unlike mgcv's `plot.gam()`, which shows the centred partial effect;
- laser power, trial index and age are held at their mean;
- the prediction is averaged over sex and cap size in proportion to the sample (sex is held fixed in the sex figure);
- subject and experiment random effects are excluded, so a curve is the expected rating for a typical subject.

The x-axis covers the 1st–99th percentile of the metric, and original units come from the mean/SD stored on each model (`mod$zscale`).

Colours are a colour-blind-safe pair (`col_sex` in `gamm_defaults.R`), validated for contrast on white. Labels use plain characters only, so they render the same on any machine.

**Held-out prediction (`holdout_replicate.R`, `robustness_core.R`, `compare_levels.R`).**
- Predictions for held-out people leave out the subject random effect, since they are new people. The experiment effect is kept for experiments seen in training.
- Held-out metrics are z-scored with the **training** mean/SD.
- The baseline is refit on the metric model's own training rows, so the two models differ only by the metric.

**Between / within split (default, `metric_decomposition = "between_within"`).**

Each metric enters the screen as two terms:
```
... + s(<metric>_between_z) + s(<metric>_within_z)
```
- `_between_z` is the subject's mean of the z-scored metric over their trials: who has high or low values.
- `_within_z` is the trial's deviation from that mean: this trial compared with the person's own average.

Why: `s(global_subjid, bs = "re")` gives each person their own intercept, which absorbs differences between people. A single pooled `s(<metric>_z)` is therefore estimated almost entirely from within-person variation and misses trait-like biomarkers.

In a synthetic test, a metric whose subject mean truly predicted pain sensitivity behaved as follows:
- pooled model: p = 0.67, so never a hit;
- split model (non-discrete fit): between p < 10⁻⁶;
- held-out prediction: 10–19% lower error in held-out people's mean ratings; LOEO, k-fold and permutation tests all in the same direction; split-half curve correlation ≈ 1.

`primary_part` (default `"between"`) decides which term's p-value goes into `edf` / `f_stat` / `p_value`, and therefore what counts as a hit. The other term is in `edf_other` / `f_other` / `p_other`. `"between"` answers the research question (does a person's typical value predict their pain sensitivity?). Within-person effects are still estimated and stored, but are not followed up unless `primary_part = "within"`.

Both decompositions can be stored side by side: every results table is keyed on `decomposition`, pooled `.rds` files and figures go to `pooled/` sub-folders, and every follow-up script reads the decomposition set in `gamm_defaults.R`.

Two cautions for between-person terms:
- **In-sample AIC barely changes.** The subject intercepts already fit each person's mean, so "AIC change vs baseline" is about 0 even for a real effect. Held-out prediction is the test.
- **N is the number of subjects, not trials.** The between term has as many data points as there are subjects.

---

## 6. Multiple Comparisons Stance

The GAMM screen is explicitly **exploratory**. No correction is applied at the screening stage. The train/holdout split is the de facto correction: a metric × ROI that reaches significance in screening must replicate in the held-out sample before it is treated as a confirmed finding.

This decision is recorded in `gamm_defaults.R` as `apply_fdr_correction = FALSE` (false discovery rate), so it is visible in configuration.

At 44 metrics × 18 ROIs, the source screen alone is 792 single-metric tests, so about 40 "hits" are expected at p < 0.05 under the null. Held-out replication applies BH across all hits; the robustness tests then probe the replicated ones.

---

## 7. Run Sequence

Scripts are run from their stage directory (`R/merge/` or `R/gamm/`) so relative `source()` calls resolve.

```
Step 0 (once, when moving to R stage V3):
                 archive results.duckdb (e.g. rename to results_pre_v3.duckdb)
                   → old gamm_fitted rows came from discrete fits, whose
                     between-subject p-values are invalid; init_results_db()
                     refuses to append to them
                 data.duckdb does NOT need archiving: merge_core.R rebuilds
                 every derived table and keeps subject_split

Step 1 (every time outputs change):   R/merge/merge_core.R
                   → rebuilds data.duckdb from the files on disk
                   → alignment check per subject → alignment_qc
                   → assigns split to new subjects only

Step 2:          R/gamm/source_core.R
                   → 44 metrics × each ROI → gamm_fitted

Step 3:          R/gamm/channel_core.R
                   → 44 metrics × each channel → gamm_fitted

Step 4:          R/gamm/sex_diff_core.R
                   → refits every significant model with a sex difference
                     curve → gamm_sex_diff

Step 5:          R/gamm/plot_gamm_hits.R
                   → <metric>_effect.png (+ <metric>_sex.png) per hit

Step 6:          R/gamm/holdout_replicate.R
                   → held-out prediction for every hit → holdout_replication
                     + holdout_replication_<level>.png

Step 7:          R/gamm/compare_levels.R
                   → channel vs source summary + paired held-out comparison

Step 8:          R/robustness/robustness_core.R   (run from R/gamm/)
                   → LOEO, k-fold, permutation, split-half, FOOOF sensitivity
                     for replicated hits → robustness + <metric>_robustness.png

Step 9 (twice):  R/validate/power_analysis.R      (run from R/gamm/)
                   → all significant source hits → power_analysis

Any time:        R/gamm/classical_tests.R
                   → ERD t-tests, rmANOVA across ROIs, Rayleigh tests
```

---

## 8. Dependencies

```r
# Core
library(readr); library(dplyr); library(tidyr); library(stringr)
library(purrr); library(tibble)
library(mgcv)       # GAMMs — bam(), select = TRUE, fREML
library(ggplot2)    # figures
library(patchwork)  # two-panel figures

# Database
library(DBI)
library(duckdb)

# Power analysis
library(lme4)       # lmer() refit for simr
library(simr)       # powerCurve(), extend(), fixef<-

# Optional
library(afex)       # rmANOVA with Greenhouse–Geisser (classical_tests.R falls back to aov())
```

---

## 9. Output Directory Structure

```
db/
├── data.duckdb          ← Section 2 tables
└── results.duckdb       ← gamm_fitted, gamm_sex_diff, holdout_replication,
                           power_analysis, robustness

R/gamm_outputs/
├── source/<roi>/<model_name>.rds, sexdiff_<metric>.rds   (pooled/ for pooled models)
├── channel/<channel>/<model_name>.rds, sexdiff_<metric>.rds
├── figures/<level>/<unit>/<metric>_effect.png, <metric>_sex.png, <metric>_robustness.png
├── figures/holdout_replication_<level>.png
├── figures/level_comparison.png
├── compare/level_summary.csv, level_paired_holdout.csv
├── classical_tests/test1_*, test2_*, test3_*, classical_tests_report.txt
└── power_analysis/power_curve_<metric>_<roi>.png
```

---

## 10. Key Design Decisions (and why)

**Check alignment before ingesting, exclude on failure.** A shifted trial index produces plausible-looking but wrong data with no error. The check is exact where it can be (two-level laser experiments), and every exclusion is recorded with a reason rather than dropped silently.

**Rebuild derived tables each run, keep the split.** The database always mirrors the files on disk: a subject whose files were deleted or re-run cannot leave stale rows behind. The split is the only state that must outlive a rebuild.

**One transaction per subject.** A file that fails half-way through ingestion rolls back that subject entirely; the database never holds a partial subject.

**Long/tidy metric tables over wide columns.** New metrics appear as new rows; no schema change, no downstream edits.

**`global_subjid` is a string (e.g. `E02_S019`), not a row number.** Row numbers shift when subjects are added; a string derived from `(experiment_id, subjid)` is stable.

**Quota-based stratified split, not independent random draws.** For the smallest cap-size group, independent draws at p = 0.75 have a 5th–95th percentile range of roughly 68–81%. Quota-based assignment converges to exactly 75/25 per cap size.

**Explicit metric selection.** The screen fits only what `metric_include` names, so a new column in the feature CSVs (for example a QC flag) can never silently become a predictor.

**`data.duckdb` opened read-only during the GAMM stage.** `connect_data_db()` enforces `read_only = TRUE`, so an accidental write fails loudly.

**Non-discrete fits.** Correct p-values for subject-level terms matter more than speed: they define the hits.

**Prediction is the criterion; every check uses the same machinery.** Held-out replication, LOEO, k-fold and the channel-vs-source comparison all compare a metric model against its nested baseline on people the fit never saw, with subject-level resampling. Permutations shuffle whole subjects, never trials, so within-subject correlation is preserved.

**Multiple comparisons: no correction at the screening stage.** The holdout replication is the correction, and the decision is recorded in configuration.