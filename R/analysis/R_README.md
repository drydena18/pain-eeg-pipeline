# R Analysis Pipeline

This page documents the R workflow: how raw behavioural data is merged, how demographics and EEG spectral features are incorporated, and how the GAMM models are structured.

---

## Overview

The R pipeline has three sequential merge steps followed by two GAMM workflows.

```
merge_behavioural.R
    ↓  behavioural_master.csv
merge_participants_into_behavioural.R
    ↓  behavioural_demo_master.csv
merge_spectral_behaviour.R
    ↓  alpha_pain_master.csv
        ├── run_gamm_alpha_metrics.R        (v1: metric-by-metric)
        └── run_gamm_alpha_metrics_v2.R     (v2: interaction + aperiodic)
```

---

## Experiment Lookup Table

All R scripts share a fixed experiment lookup table:

| experiment_name | experiment_id |
|-----------------|---------------|
| 26ByBiosemi | 1 |
| 29ByANT | 2 |
| 39ByBP | 3 |
| 30ByANT | 4 |
| 65ByANT | 5 |
| 95ByBP | 6 |
| 142ByBiosemi | 7 |
| 223ByBP | 8 |
| 29ByBP | 9 |

This lookup is reproduced in each script as a `tibble::tribble`. Experiment IDs are used to construct globally unique subject identifiers (`subjid_uid = sprintf("E%02d_S%03d", experiment_id, subjid)`).

---

## Step 1: `merge_behavioural.R`

**Input:** Individual per-experiment behavioural CSVs in `R/analysis/experiment/`

**Output:** `behavioural_master.csv`

Each CSV is matched to an experiment name by checking whether the experiment name appears in the filename. The script standardizes column names via flexible `rename(any_of(...))`, tolerating column name variation across experiments.

Required output columns:

| Column | Description |
|--------|-------------|
| `experiment_name` | Experiment identifier string |
| `experiment_id` | Integer (from lookup) |
| `subjid` | Integer subject ID within experiment |
| `subjid_uid` | Globally unique: `E01_S001` format |
| `global_subjid` | Integer row number across all subjects (used as random effect) |
| `trial` | Trial number as recorded in the raw CSV |
| `trial_index` | Within-subject sequential trial counter (1-based) |
| `laser_power` | Stimulus intensity |
| `pain_rating` | Subjective pain rating |

`trial_index` is computed as `row_number()` within each subject after sorting by trial, which handles missing or non-sequential trial numbers.

---

## Step 2: `merge_participants_into_behavioural.R`

**Inputs:**
- `behavioural_master.csv`
- One `participants.tsv` per experiment (from `da-analysis/<experiment_name>/`)
- `cap_size.csv`

**Output:** `behavioural_demo_master.csv`

### Participants TSV

Each `participants.tsv` is read from `<PROJ_ROOT>/<experiment_name>/participants.tsv`. The script handles column name variation in participant TSVs (`participant_id`, `participant`, `age`, `Age`, `sex`, `Sex`, `gender`).

Sex values are harmonized: `"m"`, `"male"` → `"M"`; `"f"`, `"female"` → `"F"`.

Subject IDs are parsed by stripping the `sub-` prefix and converting to integer.

### Cap Size CSV

`cap_size.csv` must contain: `experiment_name`, `experiment_id`, `cap_size`. One row per experiment.

Cap size is stored as a factor.

### Merge Strategy

`participants_master` is joined to `behavioural_master` on `(experiment_name, experiment_id, subjid)`. `cap_df` is joined on `(experiment_name, experiment_id)` — one cap size per experiment. Both are left joins; missing matches produce `NA`.

---

## Step 3: `merge_spectral_behaviour.R`

**Inputs:**
- `behavioural_demo_master.csv`
- All `*_ga_by_trial.csv` files found recursively under `da-analysis/`

**Output:** `alpha_pain_master.csv` (upserted incrementally)

### Spectral File Discovery

Scans `spectral_root` recursively for files matching `_ga_by_trial\.csv$`. The experiment name is extracted by parsing the file path: `da-analysis/<experiment_name>/preproc/...`.

### Per-Subject Merge

For each spectral file:
1. Parse `experiment_name` from the file path
2. Look up `experiment_id` from the fixed table
3. Join spectral features to the corresponding rows in `behavioural_demo_master` on `(experiment_name, experiment_id, subjid, subjid_uid, trial)`

Column names are normalized via `rename(any_of(...))` to tolerate naming drift.

### Incremental Upsert

If `alpha_pain_master.csv` already exists:
1. Read the existing file
2. Remove any rows matching `(experiment_id, subjid, trial)` from the incoming new data
3. Bind the new rows
4. Re-sort by `(experiment_id, subjid, trial)`

This allows re-running the merge for new subjects or corrected data without duplicating or corrupting existing rows.

---

## Step 4a: `run_gamm_alpha_metrics.R` (v1)

**Input:** `alpha_pain_master.csv`

**Output:** `gamm_outputs/`

### QC Filtering

Rows are removed if any of the following are missing: `pain_rating`, `laser_power`, `trial_index`, `global_subjid`, `experiment_id`, `age`, `sex`, `cap_size`. If FOOOF columns are present, rows with `fooof_r2 < 0.80` are also removed.

### Scaling

All continuous predictors are z-scored within the modelling dataset:
- `age_z`, `laser_power_z`, `trial_index_z`
- Any available spectral metric: `paf_cog_hz_z`, `pow_slow_z`, `pow_fast_z`, `sf_ratio_z`, `sf_logratio_z`, `sf_balance_z`, `slow_frac_z`, `fooof_offset_z`, `fooof_exponent_z`

### Model Structure

All models use `mgcv::bam` with `method = "fREML"` and `discrete = TRUE` for computational efficiency with large datasets.

**Baseline formula:**
```r
pain_rating ~
    s(laser_power_z, k = 10) +
    s(trial_index_z, k = 10) +
    s(age_z, k = 10) +
    sex +
    s(global_subjid, bs = "re") +
    s(experiment_id, bs = "re")
```

Random effects for `global_subjid` and `experiment_id` account for between-subject and between-experiment variance. The `cap_size` factor is added in `_cap` variants.

### Model Set (v1)

| Model | Added term |
|-------|-----------|
| `m00_baseline` | Baseline only |
| `m01_baseline_cap` | + `cap_size` |
| `m02_paf` | + `s(paf_cog_hz_z)` |
| `m03_pow_slow` | + `s(pow_slow_z)` |
| `m04_pow_fast` | + `s(pow_fast_z)` |
| `m05_sf_ratio` | + `s(sf_ratio_z)` |
| `m06_sf_logratio` | + `s(sf_logratio_z)` |
| `m07_sf_balance` | + `s(sf_balance_z)` |
| `m08_slow_frac` | + `s(slow_frac_z)` |
| `m09_sf_logratio_cap` | + `s(sf_logratio_z)` + `cap_size` |
| `m10_sf_balance_cap` | + `s(sf_balance_z)` + `cap_size` |
| `m11_slow_frac_cap` | + `s(slow_frac_z)` + `cap_size` |

### Outputs (v1)

- `mXX_<name>_summary.txt` — full `summary(mod)` output
- `mXX_<name>_diagnostics.png` — `gam.check()` diagnostic plots
- `mXX_<name>_observed_vs_fitted.png` — scatter of observed vs. fitted pain ratings
- `mXX_<name>.rds` — serialized model object
- `model_comparison.csv` — AIC, BIC, log-likelihood, deviance explained, N for all models
- `trial_level_fitted_values.csv` — full dataset with fitted values and residuals for all models
- `subject_level_ga_fitted_values.csv` — subject-level grand averages of fitted values
- `alpha_pain_master_model_input.csv` — filtered and scaled modelling dataset

---

## Step 4b: `run_gamm_alpha_metrics_v2.R` (v2)

Extends v1 with three additions:

1. **Derived interaction metrics** — `sf_logratio`, `sf_balance`, `slow_frac` without vs. with aperiodic controls
2. **Raw slow×fast tensor interaction** — `te(pow_slow_z, pow_fast_z)` bivariate tensor product smooth
3. **Aperiodic-controlled comparator models** — all alpha metrics tested against a baseline that includes `s(fooof_offset_z)` and `s(fooof_exponent_z)` as covariates

### Aperiodic Baseline Formula (v2)

```r
pain_rating ~
    s(laser_power_z, k = 10) +
    s(trial_index_z, k = 10) +
    s(age_z, k = 10) +
    sex +
    cap_size +
    s(fooof_offset_z, k = 10) +
    s(fooof_exponent_z, k = 10) +
    s(global_subjid, bs = "re") +
    s(experiment_id, bs = "re")
```

### Model Set (v2)

| Model | Description |
|-------|-------------|
| `m00_baseline` | Baseline (no aperiodic) |
| `m01_baseline_aperiodic` | Baseline + aperiodic offset + exponent |
| `m02_sf_logratio` | + `s(sf_logratio_z)` (no aperiodic) |
| `m03_sf_balance` | + `s(sf_balance_z)` (no aperiodic) |
| `m04_slow_frac` | + `s(slow_frac_z)` (no aperiodic) |
| `m05_tensor_slow_fast` | + `te(pow_slow_z, pow_fast_z)` |
| `m06_sf_logratio_aperiodic` | Aperiodic baseline + `s(sf_logratio_z)` |
| `m07_sf_balance_aperiodic` | Aperiodic baseline + `s(sf_balance_z)` |
| `m08_slow_frac_aperiodic` | Aperiodic baseline + `s(slow_frac_z)` |
| `m09_tensor_slow_fast_aperiodic` | Aperiodic baseline + `te(pow_slow_z, pow_fast_z)` |
| `m10_raw_main_effects_aperiodic` | Aperiodic baseline + `s(pow_slow_z)` + `s(pow_fast_z)` |

### Additional v2 Outputs

- `m05_tensor_slow_fast_surface.png` — `vis.gam` contour plot of the slow×fast interaction surface
- `m09_tensor_slow_fast_aperiodic_surface.png` — same for the aperiodic-controlled tensor model
- Smooth effect plots for all scalar alpha metrics (both raw and aperiodic-controlled variants)

---

## Random Effects Design

Both workflows use two crossed random effects:

- `s(global_subjid, bs = "re")` — accounts for between-subject differences in pain sensitivity
- `s(experiment_id, bs = "re")` — accounts for between-experiment differences (equipment, protocol, population)

`global_subjid` is an integer assigned globally across all experiments, ensuring subjects from different experiments never share the same ID. It is constructed in `merge_behavioural.R` as `row_number()` over a `distinct(experiment_id, subjid)` table.

---

## Interpreting Model Comparisons

`model_comparison.csv` (and `_v2.csv`) ranks all models by AIC. Lower AIC indicates better fit penalized for complexity.

To assess whether an alpha metric improves fit over baseline:
1. Compare `mXX_<metric>` AIC against `m00_baseline` AIC
2. Check `dev_expl` — larger values indicate more variance explained
3. Inspect the smooth effect plot — a flat smooth near zero suggests no effect

For aperiodic-controlled models, the key comparison is `m06`/`m07`/`m08` vs. `m01_baseline_aperiodic` — this tests whether alpha interaction metrics explain pain variance *beyond* what the aperiodic component already explains.