# =========================================================
# run_gamm_alpha_metrics.R
# ---------------------------------------------------------
# Pooled trial-level GAMM workflow: baseline + V1 full-epoch
# alpha metrics + V2 pre-stimulus interaction metrics.
# V 2.0.0
#
# V2.0.0 changes:
#   - Column names updated to match dynamic CSV writer output:
#       pow_slow_alpha   (was pow_slow)
#       pow_fast_alpha   (was pow_fast)
#       pow_alpha_total  (was pow_alpha)
#       rel_slow_alpha   (was rel_slow)
#       rel_fast_alpha   (was rel_fast)
#       slow_alpha_frac  (was slow_frac)
#   - Added V2 pre-stim metric model candidates:
#       bi_pre, lr_pre, cog_pre, psi_cog, delta_erd
#   - Phase (phase_slow_rad) handled via sin/cos decomposition.
#   - fitted_trial_df and subject_ga_df use any_of() so new
#     columns are included automatically as they appear.
#   - Global QC filter requires only identity + core covariates;
#     per-metric NAs are handled gracefully by fit_metric_gamm.
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
data_file        <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/alpha_pain_master.csv"
out_dir          <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_fooof_r2     <- 0.80
max_fooof_error  <- Inf
use_discrete     <- TRUE
nthreads_to_use  <- 4

# =========================================================
# HELPERS
# =========================================================
clean_names_local <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2", "r2") %>%
    str_replace_all("\\.", "_")
}

zscore     <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))

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
# Circular predictors cannot enter a GAMM as raw radians.
# Decompose into sin and cos components (each bounded [-1,1])
# which are used as linear additive terms.
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
# globally.  Spectral metric NAs are handled per-model.
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
# Only z-scores columns that are present in df_model.
# =========================================================
metric_candidates <- c(
  # V1 full-epoch metrics (new canonical names)
  "paf_cog_hz",
  "pow_slow_alpha",
  "pow_fast_alpha",
  "pow_alpha_total",
  "rel_slow_alpha",
  "rel_fast_alpha",
  "sf_ratio",
  "sf_logratio",
  "sf_balance",
  "slow_alpha_frac",
  # V2 pre-stim interaction metrics
  "bi_pre",
  "lr_pre",
  "cog_pre",
  "psi_cog",
  "delta_erd",
  # FOOOF aperiodic
  "fooof_offset",
  "fooof_exponent"
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

# Phase components are already bounded; z-score for comparability
if ("phase_sin" %in% names(df_model)) {
  df_model$phase_sin_z <- zscore(df_model$phase_sin)
  df_model$phase_cos_z <- zscore(df_model$phase_cos)
}

write_csv(df_model, file.path(out_dir, "alpha_pain_master_model_input.csv"))

message("Model input rows: ", nrow(df_model))
message("Subjects: ",         n_distinct(df_model$subjid_uid))

# =========================================================
# BASE FORMULAS
# =========================================================
baseline_formula <- pain_rating ~
  s(laser_power_z, k = 10) +
  s(trial_index_z, k = 10) +
  s(age_z, k = 10) +
  sex +
  s(global_subjid, bs = "re") +
  s(experiment_id, bs = "re")

baseline_cap_formula <- pain_rating ~
  s(laser_power_z, k = 10) +
  s(trial_index_z, k = 10) +
  s(age_z, k = 10) +
  sex +
  cap_size +
  s(global_subjid, bs = "re") +
  s(experiment_id, bs = "re")

# =========================================================
# MODEL FITTER
# Returns NULL (with warning) when the predictor column is
# absent or is all-NA after the global filter.
# =========================================================
fit_metric_gamm <- function(data, metric_z, model_name,
                             include_cap = FALSE, extra_terms = NULL) {
  base_formula <- if (include_cap) baseline_cap_formula else baseline_formula

  if (!is.na(metric_z)) {
    if (!metric_z %in% names(data)) {
      warning("Skipping ", model_name, ": column not found -> ", metric_z)
      return(NULL)
    }
    if (all(is.na(data[[metric_z]]))) {
      warning("Skipping ", model_name, ": all-NA -> ", metric_z)
      return(NULL)
    }
    formula_to_fit <- update(
      base_formula,
      paste(". ~ . + s(", metric_z, ", k = 10)")
    )
  } else {
    formula_to_fit <- base_formula
  }

  # Extra terms (e.g. phase sin + cos as linear additive)
  if (!is.null(extra_terms)) {
    for (term in extra_terms) {
      formula_to_fit <- update(formula_to_fit,
                                paste(". ~ . +", term))
    }
  }

  # Drop rows NA in the specific predictor before fitting
  fit_data <- if (!is.na(metric_z)) {
    data %>% filter(!is.na(.data[[metric_z]]))
  } else {
    data
  }

  if (nrow(fit_data) == 0) {
    warning("Skipping ", model_name, ": no rows after NA removal.")
    return(NULL)
  }

  message("Fitting: ", model_name, "  (n=", nrow(fit_data), ")")

  tryCatch(
    bam(
      formula  = formula_to_fit,
      data     = fit_data,
      method   = "fREML",
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
# NA in metric_z = baseline model (no spectral predictor).
# include_cap adds cap_size as a fixed effect.
# =========================================================
model_specs <- tribble(
  ~model_name,                    ~metric_z,              ~include_cap,
  # --- Baselines ---
  "m00_baseline",                 NA_character_,           FALSE,
  "m01_baseline_cap",             NA_character_,           TRUE,
  # --- V1 full-epoch metrics ---
  "m02_paf_cog",                  "paf_cog_hz_z",          FALSE,
  "m03_pow_slow",                 "pow_slow_alpha_z",      FALSE,
  "m04_pow_fast",                 "pow_fast_alpha_z",      FALSE,
  "m05_sf_ratio",                 "sf_ratio_z",            FALSE,
  "m06_sf_logratio",              "sf_logratio_z",         FALSE,
  "m07_sf_balance",               "sf_balance_z",          FALSE,
  "m08_slow_frac",                "slow_alpha_frac_z",     FALSE,
  "m09_sf_logratio_cap",          "sf_logratio_z",         TRUE,
  "m10_sf_balance_cap",           "sf_balance_z",          TRUE,
  "m11_slow_frac_cap",            "slow_alpha_frac_z",     TRUE,
  # --- V2 pre-stim interaction metrics ---
  "m12_bi_pre",                   "bi_pre_z",              FALSE,
  "m13_lr_pre",                   "lr_pre_z",              FALSE,
  "m14_cog_pre",                  "cog_pre_z",             FALSE,
  "m15_psi_cog",                  "psi_cog_z",             FALSE,
  "m16_delta_erd",                "delta_erd_z",           FALSE,
  "m17_bi_pre_cap",               "bi_pre_z",              TRUE,
  "m18_psi_cog_cap",              "psi_cog_z",             TRUE
)

# Phase model handled separately (two linear terms, no smooth z column)
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
  model_list[["m19_phase_sincos"]] <- fit_metric_gamm(
    data        = df_model %>% filter(!is.na(phase_sin_z), !is.na(phase_cos_z)),
    metric_z    = NA_character_,
    model_name  = "m19_phase_sincos",
    include_cap = FALSE,
    extra_terms = c("phase_sin_z", "phase_cos_z")
  )
}

# Drop failed models
valid_idx  <- !map_lgl(model_list, is.null)
model_list <- model_list[valid_idx]

if (length(model_list) == 0) stop("No GAMMs were successfully fit.")
message("Successfully fit models: ", paste(names(model_list), collapse = ", "))

# =========================================================
# SAVE SUMMARIES
# =========================================================
iwalk(model_list, function(mod, nm) {
  capture.output(summary(mod),
                 file = file.path(out_dir, paste0(nm, "_summary.txt")))
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
# FITTED VALUES  (all spectral columns included via any_of)
# =========================================================
spectral_passthrough <- c(
  "paf_cog_hz", "pow_slow_alpha", "pow_fast_alpha", "pow_alpha_total",
  "rel_slow_alpha", "rel_fast_alpha", "sf_ratio", "sf_logratio",
  "sf_balance", "slow_alpha_frac",
  "bi_pre", "lr_pre", "cog_pre", "psi_cog", "erd_slow", "erd_fast",
  "delta_erd", "p5_flag", "phase_slow_rad",
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
  # Fitted values may be on a subset of rows; align by row name
  fit_vals  <- rep(NA_real_, nrow(fitted_trial_df))
  resid_vals <- rep(NA_real_, nrow(fitted_trial_df))
  mod_rows  <- as.integer(rownames(model_list[[nm]]$model))
  # bam() preserves original row indices in model$model
  common    <- intersect(mod_rows, seq_len(nrow(fitted_trial_df)))
  if (length(common) > 0) {
    fit_vals[common]   <- fitted(model_list[[nm]])
    resid_vals[common] <- residuals(model_list[[nm]])
  }
  fitted_trial_df[[paste0(nm, "_fitted")]] <- fit_vals
  fitted_trial_df[[paste0(nm, "_resid")]]  <- resid_vals
}

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values.csv"))

# Subject-level GA summary
subject_ga_df <- fitted_trial_df %>%
  group_by(experiment_name, experiment_id, subjid, subjid_uid,
           global_subjid, age, sex, cap_size) %>%
  summarise(
    n_trials         = n(),
    mean_pain_rating = mean(pain_rating, na.rm = TRUE),
    mean_laser_power = mean(laser_power, na.rm = TRUE),
    across(any_of(spectral_passthrough), ~ mean(.x, na.rm = TRUE),
           .names = "mean_{.col}"),
    across(ends_with("_fitted"),         ~ mean(.x, na.rm = TRUE),
           .names = "ga_{.col}"),
    .groups = "drop"
  )

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values.csv"))

# =========================================================
# DIAGNOSTICS + OBSERVED VS FITTED
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
    labs(title = paste("Observed vs Fitted:", nm),
         x = "Fitted pain rating", y = "Observed pain rating") +
    theme_minimal()

  ggsave(file.path(out_dir, paste0(nm, "_observed_vs_fitted.png")),
         plot = p, width = 7, height = 5, dpi = 300)
}

# =========================================================
# SAVE MODEL OBJECTS
# =========================================================
iwalk(model_list, function(mod, nm) {
  saveRDS(mod, file.path(out_dir, paste0(nm, ".rds")))
})

message("Done. Outputs saved to: ", out_dir)