# =========================================================
# run_gamm_alpha_metrics.R
# ---------------------------------------------------------
# Pooled trial-level GAMM workflow: baseline + V1 full-epoch
# alpha metrics + V2 pre-stimulus interaction metrics.
# V 3.0.0
#
# V2.0.0 changes:
#   - Column names updated to match dynamic CSV writer output.
#   - Added V2 pre-stim metric model candidates.
#   - Phase (phase_slow_rad) handled via sin/cos decomposition.
#   - fitted_trial_df and subject_ga_df use any_of().
#   - Global QC filter requires only identity + core covariates.
#
# V3.0.0 changes:
#   - safe_k(): k is computed from actual unique values in the
#     fit data, clamped to [k_min, k_max]. Prevents mgcv error
#     "fewer unique covariate combinations than specified maximum
#     degrees of freedom" when laser_power or other covariates
#     have fewer unique values than the hardcoded k=10.
#   - build_baseline_formula(): formula is constructed
#     dynamically from fit_data after per-metric NA removal, so
#     adaptive k is always computed on the real sample.
#     Also guards against single-level sex and cap_size factors,
#     which cause "contrasts not defined for 0 degrees of freedom"
#     when only one sex or one cap size is present in the data.
#   - select=TRUE added to all bam() calls. This is NOT the same
#     as k=-1 (which just uses mgcv's default k). select=TRUE
#     applies a double penalty that can shrink individual smooth
#     terms toward zero if the data do not support them, giving
#     automatic smooth term selection via fREML without requiring
#     manual model comparison for each term.
#   - model_specs uses tibble::tibble() with explicit vectors
#     instead of tribble(), which fails on older tibble versions.
#   - age_z is a linear fixed effect (not a smooth) because age
#     is constant within subject and cannot support a nonlinear
#     smooth with small N; restore s(age_z) once N >> 50.
# =========================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(mgcv)
library(tibble)
library(ggplot2)

# =========================================================
# USER SETTINGS
# =========================================================
data_file       <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/alpha_pain_master.csv"
out_dir         <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/gamms_v1"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_fooof_r2    <- 0.80
max_fooof_error <- Inf
use_discrete    <- TRUE
nthreads_to_use <- 4

# =========================================================
# HELPERS
# =========================================================
clean_names_local <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2",  "r2") %>%
    str_replace_all("\\.",   "_")
}

zscore      <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))

# safe_k: returns the largest k that won't exceed the number of
# unique values in x, clamped to [k_min, k_max].
# k is the upper bound on smooth complexity; the actual effective
# degrees of freedom are estimated from data via fREML penalisation.
safe_k <- function(x, k_max = 10, k_min = 3) {
  n_unique <- length(unique(na.omit(x)))
  max(k_min, min(k_max, n_unique - 1L))
}

# build_baseline_formula: constructs the baseline formula from the
# actual fit_data so that adaptive k and factor-level guards are
# applied to the real sample, not the full df_model.
build_baseline_formula <- function(data, include_cap = FALSE) {
  k_laser <- safe_k(data$laser_power_z)
  k_trial <- safe_k(data$trial_index_z)

  # Fixed effects: only include factors with >1 level
  fixed_terms <- "age_z"
  if (nlevels(data$sex) > 1)
    fixed_terms <- c(fixed_terms, "sex")
  if (include_cap && nlevels(data$cap_size) > 1)
    fixed_terms <- c(fixed_terms, "cap_size")

  rhs <- paste(
    c(
      sprintf("s(laser_power_z, k = %d)", k_laser),
      sprintf("s(trial_index_z, k = %d)", k_trial),
      fixed_terms,
      "s(global_subjid, bs = 're')",
      "s(experiment_id,  bs = 're')"
    ),
    collapse = " + "
  )

  as.formula(paste("pain_rating ~", rhs))
}

# =========================================================
# LOAD DATA
# =========================================================
df <- read_csv(data_file, show_col_types = FALSE)
names(df) <- clean_names_local(names(df))

required_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
  "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating"
)

missing_required <- setdiff(required_cols, names(df))
if (length(missing_required) > 0) {
  stop("Missing required columns in alpha_pain_master.csv: ",
       paste(missing_required, collapse = ", "))
}

df <- df %>%
  mutate(
    experiment_name = safe_factor(experiment_name),
    experiment_id   = safe_factor(experiment_id),
    subjid          = as.integer(subjid),
    subjid_uid      = safe_factor(subjid_uid),
    global_subjid   = safe_factor(global_subjid),
    age             = suppressWarnings(as.numeric(age)),
    sex             = safe_factor(sex),
    cap_size        = safe_factor(cap_size),
    trial           = as.integer(trial),
    trial_index     = as.integer(trial_index),
    laser_power     = suppressWarnings(as.numeric(laser_power)),
    pain_rating     = suppressWarnings(as.numeric(pain_rating))
  )

# =========================================================
# PHASE DECOMPOSITION
# Circular predictors (radians) cannot enter a GAMM directly.
# Decompose into sin and cos components (bounded [-1, 1]) as
# linear additive terms — a standard approach for phase data.
# =========================================================
if ("phase_slow_rad" %in% names(df)) {
  df <- df %>%
    mutate(
      phase_sin = sin(phase_slow_rad),
      phase_cos = cos(phase_slow_rad)
    )
}

# =========================================================
# QC FILTERING
# Only required identity + covariate columns are filtered
# globally. Per-metric NAs are handled inside fit_metric_gamm.
# =========================================================
df_model <- df %>%
  filter(
    !is.na(pain_rating),
    !is.na(laser_power),
    !is.na(trial_index),
    !is.na(global_subjid),
    !is.na(experiment_id),
    !is.na(age),
    !is.na(sex),
    !is.na(cap_size)
  )

if ("fooof_r2" %in% names(df_model)) {
  df_model <- df_model %>%
    filter(is.na(fooof_r2) | fooof_r2 >= min_fooof_r2)
}
if ("fooof_error" %in% names(df_model) && is.finite(max_fooof_error)) {
  df_model <- df_model %>%
    filter(is.na(fooof_error) | fooof_error <= max_fooof_error)
}

# =========================================================
# SCALE CONTINUOUS VARIABLES
# Z-scores only columns present in df_model.
# =========================================================
metric_candidates <- c(
  # V1 full-epoch metrics (canonical names from spec_write_ga_trial_csv)
  "paf_cog_hz",
  "pow_slow_alpha", "pow_fast_alpha", "pow_alpha_total",
  "rel_slow_alpha", "rel_fast_alpha",
  "sf_ratio", "sf_logratio", "sf_balance", "slow_alpha_frac",
  # V2 pre-stimulus interaction metrics
  "bi_pre", "lr_pre", "cog_pre", "psi_cog", "delta_erd",
  # FOOOF aperiodic
  "fooof_offset", "fooof_exponent"
)

df_model <- df_model %>%
  mutate(
    age_z         = zscore(age),
    laser_power_z = zscore(laser_power),
    trial_index_z = zscore(trial_index)
  )

for (metric in metric_candidates) {
  if (metric %in% names(df_model)) {
    df_model[[paste0(metric, "_z")]] <- zscore(df_model[[metric]])
  }
}

if ("phase_sin" %in% names(df_model)) {
  df_model$phase_sin_z <- zscore(df_model$phase_sin)
  df_model$phase_cos_z <- zscore(df_model$phase_cos)
}

write_csv(df_model, file.path(out_dir, "alpha_pain_master_model_input.csv"))

message("Model input rows: ", nrow(df_model))
message("Subjects: ",         n_distinct(df_model$subjid_uid))
message("Experiments: ",      n_distinct(df_model$experiment_id))
message("Unique laser_power values: ",
        length(unique(na.omit(df_model$laser_power))))

# =========================================================
# MODEL FITTER
# ---------------------------------------------------------
# Per-metric NA rows are dropped BEFORE formula construction
# so that safe_k and factor-level guards operate on the
# actual sample being fit, not the full df_model.
#
# select=TRUE enables automatic smooth term selection: each
# smooth gets a second penalty that can shrink it to zero.
# This is the correct mechanism for adaptive complexity —
# NOT k=-1, which merely uses mgcv's default k.
# =========================================================
fit_metric_gamm <- function(data, metric_z, model_name,
                             include_cap = FALSE, extra_terms = NULL) {

  # Guard: cap models need at least two cap levels
  if (include_cap && nlevels(data$cap_size) < 2) {
    warning("Skipping ", model_name,
            ": cap_size has only 1 level in current data.")
    return(NULL)
  }

  # Guard: metric column must exist and not be all-NA
  if (!is.na(metric_z)) {
    if (!metric_z %in% names(data)) {
      warning("Skipping ", model_name, ": column not found -> ", metric_z)
      return(NULL)
    }
    if (all(is.na(data[[metric_z]]))) {
      warning("Skipping ", model_name, ": all-NA -> ", metric_z)
      return(NULL)
    }
  }

  # Drop per-metric NAs BEFORE building formula
  fit_data <- if (!is.na(metric_z)) {
    data %>% filter(!is.na(.data[[metric_z]]))
  } else {
    data
  }

  if (nrow(fit_data) == 0) {
    warning("Skipping ", model_name, ": no rows after NA removal.")
    return(NULL)
  }

  # Build baseline formula from actual fit_data (adaptive k, level guards)
  base_formula <- build_baseline_formula(fit_data, include_cap)

  # Add metric smooth with adaptive k
  if (!is.na(metric_z)) {
    k_metric      <- safe_k(fit_data[[metric_z]])
    formula_to_fit <- update(
      base_formula,
      as.formula(sprintf(". ~ . + s(%s, k = %d)", metric_z, k_metric))
    )
  } else {
    formula_to_fit <- base_formula
  }

  # Extra linear terms (e.g., phase sin + cos)
  if (!is.null(extra_terms)) {
    for (term in extra_terms) {
      formula_to_fit <- update(formula_to_fit,
                               as.formula(paste(". ~ . +", term)))
    }
  }

  message("Fitting: ", model_name, "  (n=", nrow(fit_data), ")")

  tryCatch(
    bam(
      formula  = formula_to_fit,
      data     = fit_data,
      method   = "fREML",
      select   = TRUE,           # automatic smooth term selection
      discrete = use_discrete,
      nthreads = nthreads_to_use
    ),
    error = function(e) {
      warning("Model failed [", model_name, "]: ", e$message)
      NULL
    }
  )
}

# =========================================================
# MODEL SPECIFICATIONS
# NA in metric_z → baseline model (no spectral predictor).
# include_cap adds cap_size as a fixed effect (skipped
# automatically when cap_size has only one level).
# =========================================================
model_specs <- tibble::tibble(
  model_name = c(
    # Baselines
    "m00_baseline",        "m01_baseline_cap",
    # V1 full-epoch metrics
    "m02_paf_cog",         "m03_pow_slow",         "m04_pow_fast",
    "m05_sf_ratio",        "m06_sf_logratio",       "m07_sf_balance",
    "m08_slow_frac",
    # V1 cap-controlled
    "m09_sf_logratio_cap", "m10_sf_balance_cap",    "m11_slow_frac_cap",
    # V2 pre-stimulus interaction metrics
    "m12_bi_pre",          "m13_lr_pre",            "m14_cog_pre",
    "m15_psi_cog",         "m16_delta_erd",
    # V2 cap-controlled
    "m17_bi_pre_cap",      "m18_psi_cog_cap"
  ),
  metric_z = c(
    NA_character_,         NA_character_,
    "paf_cog_hz_z",        "pow_slow_alpha_z",      "pow_fast_alpha_z",
    "sf_ratio_z",          "sf_logratio_z",          "sf_balance_z",
    "slow_alpha_frac_z",
    "sf_logratio_z",       "sf_balance_z",           "slow_alpha_frac_z",
    "bi_pre_z",            "lr_pre_z",               "cog_pre_z",
    "psi_cog_z",           "delta_erd_z",
    "bi_pre_z",            "psi_cog_z"
  ),
  include_cap = c(
    FALSE, TRUE,
    FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE,
    TRUE,  TRUE,  TRUE,
    FALSE, FALSE, FALSE, FALSE, FALSE,
    TRUE,  TRUE
  )
)

# =========================================================
# FIT MODELS
# =========================================================
model_list <- pmap(
  model_specs,
  function(model_name, metric_z, include_cap) {
    fit_metric_gamm(
      data        = df_model,
      metric_z    = metric_z,
      model_name  = model_name,
      include_cap = include_cap
    )
  }
)
names(model_list) <- model_specs$model_name

# Phase model: sin(phase) + cos(phase) as linear additive terms
if (all(c("phase_sin_z", "phase_cos_z") %in% names(df_model))) {
  phase_data <- df_model %>% filter(!is.na(phase_sin_z), !is.na(phase_cos_z))
  model_list[["m19_phase_sincos"]] <- fit_metric_gamm(
    data        = phase_data,
    metric_z    = NA_character_,
    model_name  = "m19_phase_sincos",
    include_cap = FALSE,
    extra_terms = c("phase_sin_z", "phase_cos_z")
  )
}

# Drop failed / skipped models
valid_idx  <- !map_lgl(model_list, is.null)
model_list <- model_list[valid_idx]

if (length(model_list) == 0) stop("No GAMMs were successfully fit.")
message("Successfully fit models: ", paste(names(model_list), collapse = ", "))

# =========================================================
# SAVE SUMMARIES
# =========================================================
iwalk(model_list, function(mod, nm) {
  capture.output(
    summary(mod),
    file = file.path(out_dir, paste0(nm, "_summary.txt"))
  )
})

# =========================================================
# MODEL COMPARISON TABLE
# =========================================================
extract_model_info <- function(mod, nm) {
  tibble(
    model_name = nm,
    AIC        = AIC(mod),
    BIC        = BIC(mod),
    logLik     = as.numeric(logLik(mod)),
    dev_expl   = summary(mod)$dev.expl,
    n          = nobs(mod)
  )
}

model_comparison <- imap_dfr(model_list, extract_model_info) %>%
  arrange(AIC)

write_csv(model_comparison, file.path(out_dir, "model_comparison.csv"))
message("Model comparison saved.")

# =========================================================
# FITTED VALUES
# Spectral columns are passed through via any_of() so new
# metrics appear automatically as they enter the pipeline.
# Fitted/residual values are aligned by original row index
# because bam() may have been fit on a metric-NA-dropped
# subset; unfitted rows receive NA.
# =========================================================
spectral_passthrough <- c(
  "paf_cog_hz",
  "pow_slow_alpha", "pow_fast_alpha", "pow_alpha_total",
  "rel_slow_alpha", "rel_fast_alpha",
  "sf_ratio", "sf_logratio", "sf_balance", "slow_alpha_frac",
  "bi_pre", "lr_pre", "cog_pre", "psi_cog",
  "erd_slow", "erd_fast", "delta_erd", "p5_flag",
  "phase_slow_rad",
  "fooof_offset", "fooof_exponent", "fooof_r2", "fooof_error",
  "fooof_alpha_cf", "fooof_alpha_pw", "fooof_alpha_bw"
)

fitted_trial_df <- df_model %>%
  select(
    experiment_name, experiment_id, subjid, subjid_uid, global_subjid,
    age, sex, cap_size, trial, trial_index, laser_power, pain_rating,
    any_of(spectral_passthrough)
  )

for (nm in names(model_list)) {
  fit_vals   <- rep(NA_real_, nrow(fitted_trial_df))
  resid_vals <- rep(NA_real_, nrow(fitted_trial_df))
  mod_rows   <- as.integer(rownames(model_list[[nm]]$model))
  common     <- intersect(mod_rows, seq_len(nrow(fitted_trial_df)))
  if (length(common) > 0) {
    fit_vals[common]   <- fitted(model_list[[nm]])
    resid_vals[common] <- residuals(model_list[[nm]])
  }
  fitted_trial_df[[paste0(nm, "_fitted")]] <- fit_vals
  fitted_trial_df[[paste0(nm, "_resid")]]  <- resid_vals
}

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values.csv"))

# Subject-level grand-average summary
subject_ga_df <- fitted_trial_df %>%
  group_by(
    experiment_name, experiment_id, subjid, subjid_uid,
    global_subjid, age, sex, cap_size
  ) %>%
  summarise(
    n_trials         = n(),
    mean_pain_rating = mean(pain_rating,   na.rm = TRUE),
    mean_laser_power = mean(laser_power,   na.rm = TRUE),
    across(any_of(spectral_passthrough),
           ~ mean(.x, na.rm = TRUE),
           .names = "mean_{.col}"),
    across(ends_with("_fitted"),
           ~ mean(.x, na.rm = TRUE),
           .names = "ga_{.col}"),
    .groups = "drop"
  )

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values.csv"))

# =========================================================
# DIAGNOSTICS + OBSERVED VS FITTED PLOTS
# =========================================================
iwalk(model_list, function(mod, nm) {
  png(file.path(out_dir, paste0(nm, "_diagnostics.png")),
      width = 1800, height = 1400, res = 180)
  gam.check(mod)
  dev.off()
})

for (nm in names(model_list)) {
  fitted_col <- paste0(nm, "_fitted")
  if (!fitted_col %in% names(fitted_trial_df)) next
  p <- ggplot(
    fitted_trial_df %>% filter(!is.na(.data[[fitted_col]])),
    aes(x = .data[[fitted_col]], y = pain_rating)
  ) +
    geom_point(alpha = 0.15) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste("Observed vs Fitted:", nm),
      x     = "Fitted pain rating",
      y     = "Observed pain rating"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    file.path(out_dir, paste0(nm, "_observed_vs_fitted.png")),
    plot   = p,
    width  = 7,
    height = 5,
    dpi    = 300
  )
}

# =========================================================
# SAVE MODEL OBJECTS
# =========================================================
iwalk(model_list, function(mod, nm) {
  saveRDS(mod, file.path(out_dir, paste0(nm, ".rds")))
})

message("Done. Outputs saved to: ", out_dir)