# =========================================================
# run_gamm_alpha_metrics_v2.R
# ---------------------------------------------------------
# V2 pooled GAMM workflow: raw slow-fast tensor interaction,
# aperiodic-controlled comparators, and V2 pre-stim
# interaction metrics (BI_pre, CoG_pre, psi_cog, DELTA_ERD,
# phase sin/cos).
# V 2.0.0
#
# V2.0.0 changes:
#   - Fixed typos: an.factor -> as.factor, use_dicrete ->
#     use_discrete, paf_cof_hz_z -> paf_cog_hz_z.
#   - Updated column names throughout (pow_slow_alpha, etc.)
#   - Required cols and QC filter use new canonical names.
#   - plot_smooth: term identification uses x$term[1] instead
#     of x$term (safe for both univariate and tensor smooths).
#   - Added V2 interaction metric models (m11-m17) covering
#     bi_pre, cog_pre, psi_cog, delta_erd, phase sin/cos,
#     and aperiodic-controlled versions of key V2 metrics.
#   - Refactored to tribble + pmap model spec table (matches
#     v1 structure; explicit bam() calls removed).
#   - fitted_trial_df and subject_ga_df use any_of() so new
#     spectral columns flow through automatically.
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
data_file       <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/alpha_pain_master.csv"
out_dir         <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_v2"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_fooof_r2    <- 0.80
max_fooof_error <- Inf
use_discrete    <- TRUE          # was: use_dicrete = TRUE (typo + wrong operator)
nthreads_to_use <- 4

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

zscore      <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))   # was: an.factor (typo)

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

# plot_smooth: find the smooth index by matching x$term[1].
# Using x$term[1] (not x$term) is safe for both s() and te()
# smooths — te() sets x$term to a length-2+ vector, and
# x$label to something like "te(a,b)"; neither equals a
# single term string reliably across mgcv versions.
plot_smooth <- function(model, term, file_name) {
  smooth_idx <- which(
    sapply(model$smooth, function(x) x$term[1]) == term
  )
  if (length(smooth_idx) == 0) {
    warning("plot_smooth: term '", term, "' not found in model.")
    return(invisible(NULL))
  }
  png(file.path(out_dir, file_name), width = 1600, height = 1200, res = 200)
  plot(model,
       select    = smooth_idx[1],
       shade     = TRUE,
       shade.col = "lightblue",
       main      = paste("Effect of", term, "on Pain Rating"),
       xlab      = term,
       ylab      = "Partial Effect on Pain")
  abline(h = 0, lty = 2)
  dev.off()
}

# =========================================================
# LOAD DATA
# =========================================================
df <- read_csv(data_file, show_col_types = FALSE)
names(df) <- clean_names_local(names(df))

# V1 columns required for all V2 raw-power / tensor models
required_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
  "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating",
  "pow_slow_alpha", "pow_fast_alpha",          # updated from pow_slow, pow_fast
  "sf_logratio", "sf_balance",
  "slow_alpha_frac"                            # updated from slow_frac
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

# Phase decomposition for circular predictor
if ("phase_slow_rad" %in% names(df)) {
  df <- df %>%
    mutate(
      phase_sin = sin(phase_slow_rad),
      phase_cos = cos(phase_slow_rad)
    )
}

# =========================================================
# QC FILTERING
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
    !is.na(cap_size),
    !is.na(pow_slow_alpha),          # required for tensor/raw-power models
    !is.na(pow_fast_alpha),
    !is.na(sf_logratio),
    !is.na(sf_balance),
    !is.na(slow_alpha_frac)
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
# =========================================================
df_model <- df_model %>%
  mutate(
    age_z           = zscore(age),
    laser_power_z   = zscore(laser_power),
    trial_index_z   = zscore(trial_index),
    pow_slow_alpha_z = zscore(pow_slow_alpha),   # updated name
    pow_fast_alpha_z = zscore(pow_fast_alpha),   # updated name
    sf_logratio_z   = zscore(sf_logratio),
    sf_balance_z    = zscore(sf_balance),
    slow_alpha_frac_z = zscore(slow_alpha_frac)  # updated name
  )

# Optional columns — z-score only when present
optional_z <- c("paf_cog_hz", "fooof_offset", "fooof_exponent",
                 "bi_pre", "lr_pre", "cog_pre", "psi_cog", "delta_erd",
                 "phase_sin", "phase_cos")

for (col in optional_z) {
  if (col %in% names(df_model) && !all(is.na(df_model[[col]]))) {
    df_model[[paste0(col, "_z")]] <- zscore(df_model[[col]])
  }
}

# Fix typo from v1: paf_cof_hz_z -> paf_cog_hz_z (already handled above)

write_csv(df_model, file.path(out_dir, "alpha_pain_master_model_input_v2.csv"))
message("Model input rows: ", nrow(df_model),
        "  Subjects: ", n_distinct(df_model$subjid_uid))

# =========================================================
# BASE FORMULAS
# =========================================================
base_formula <- pain_rating ~
  s(laser_power_z, k = 10) +
  s(trial_index_z, k = 10) +
  s(age_z, k = 10) +
  sex +
  cap_size +
  s(global_subjid, bs = "re") +
  s(experiment_id, bs = "re")

base_formula_aperiodic <- pain_rating ~
  s(laser_power_z, k = 10) +
  s(trial_index_z, k = 10) +
  s(age_z, k = 10) +
  sex +
  cap_size +
  s(fooof_offset_z, k = 10) +
  s(fooof_exponent_z, k = 10) +
  s(global_subjid, bs = "re") +
  s(experiment_id, bs = "re")

# =========================================================
# MODEL FITTER
# Gracefully skips when predictor is absent or all-NA.
# Uses aperiodic base when use_aperiodic = TRUE and FOOOF
# columns are available; falls back to standard base if not.
# =========================================================
fit_v2_gamm <- function(data, metric_z, model_name,
                         use_aperiodic = FALSE, extra_terms = NULL) {

  has_aperiodic <- use_aperiodic &&
    all(c("fooof_offset_z", "fooof_exponent_z") %in% names(data)) &&
    !all(is.na(data$fooof_offset_z))

  base <- if (has_aperiodic) base_formula_aperiodic else base_formula

  if (!is.na(metric_z)) {
    if (!metric_z %in% names(data)) {
      warning("Skipping ", model_name, ": column not found -> ", metric_z)
      return(NULL)
    }
    if (all(is.na(data[[metric_z]]))) {
      warning("Skipping ", model_name, ": all-NA -> ", metric_z)
      return(NULL)
    }
    formula_to_fit <- update(base,
                              paste(". ~ . + s(", metric_z, ", k = 10)"))
    data <- data %>% filter(!is.na(.data[[metric_z]]))
  } else {
    formula_to_fit <- base
  }

  if (!is.null(extra_terms)) {
    for (term in extra_terms) {
      formula_to_fit <- update(formula_to_fit, paste(". ~ . +", term))
    }
  }

  if (nrow(data) == 0) {
    warning("Skipping ", model_name, ": no rows after NA removal.")
    return(NULL)
  }

  message("Fitting: ", model_name, "  (n=", nrow(data), ")")
  tryCatch(
    bam(formula  = formula_to_fit,
        data     = data,
        method   = "fREML",
        discrete = use_discrete,
        nthreads = nthreads_to_use),
    error = function(e) {
      warning("Model failed [", model_name, "]: ", e$message)
      NULL
    }
  )
}

# =========================================================
# MODEL SPECIFICATIONS
# Columns: model_name | metric_z | use_aperiodic
# =========================================================
model_specs <- tribble(
  ~model_name,                         ~metric_z,               ~use_aperiodic,
  # --- Baselines ---
  "m00_baseline",                      NA_character_,            FALSE,
  "m01_baseline_aperiodic",            NA_character_,            TRUE,
  # --- V1 full-epoch: standard ---
  "m02_sf_logratio",                   "sf_logratio_z",          FALSE,
  "m03_sf_balance",                    "sf_balance_z",           FALSE,
  "m04_slow_frac",                     "slow_alpha_frac_z",      FALSE,
  # --- V1 full-epoch: aperiodic-controlled ---
  "m05_sf_logratio_aperiodic",         "sf_logratio_z",          TRUE,
  "m06_sf_balance_aperiodic",          "sf_balance_z",           TRUE,
  "m07_slow_frac_aperiodic",           "slow_alpha_frac_z",      TRUE,
  # --- V2 pre-stim interaction metrics: standard ---
  "m08_bi_pre",                        "bi_pre_z",               FALSE,
  "m09_cog_pre",                       "cog_pre_z",              FALSE,
  "m10_psi_cog",                       "psi_cog_z",              FALSE,
  "m11_delta_erd",                     "delta_erd_z",            FALSE,
  # --- V2 pre-stim interaction metrics: aperiodic-controlled ---
  "m12_bi_pre_aperiodic",              "bi_pre_z",               TRUE,
  "m13_psi_cog_aperiodic",             "psi_cog_z",              TRUE,
  "m14_delta_erd_aperiodic",           "delta_erd_z",            TRUE
)

model_list <- pmap(
  model_specs,
  function(model_name, metric_z, use_aperiodic) {
    fit_v2_gamm(
      data          = df_model,
      metric_z      = metric_z,
      model_name    = model_name,
      use_aperiodic = use_aperiodic
    )
  }
)
names(model_list) <- model_specs$model_name

# Tensor (slow x fast) models — added outside tribble because
# te() syntax can't be expressed as a single metric_z string
tensor_data <- df_model  # already filtered for pow_slow_alpha/pow_fast_alpha NAs

model_list[["m15_tensor_slow_fast"]] <- tryCatch(
  { message("Fitting: m15_tensor_slow_fast  (n=", nrow(tensor_data), ")")
    bam(update(base_formula,
               . ~ . + te(pow_slow_alpha_z, pow_fast_alpha_z, k = c(10, 10))),
        data = tensor_data, method = "fREML",
        discrete = use_discrete, nthreads = nthreads_to_use) },
  error = function(e) { warning("m15 failed: ", e$message); NULL }
)

model_list[["m16_tensor_slow_fast_aperiodic"]] <- tryCatch(
  { message("Fitting: m16_tensor_slow_fast_aperiodic  (n=", nrow(tensor_data), ")")
    bam(update(base_formula_aperiodic,
               . ~ . + te(pow_slow_alpha_z, pow_fast_alpha_z, k = c(10, 10))),
        data = tensor_data %>% filter(!is.na(fooof_offset_z)),
        method = "fREML", discrete = use_discrete, nthreads = nthreads_to_use) },
  error = function(e) { warning("m16 failed: ", e$message); NULL }
)

model_list[["m17_raw_main_effects_aperiodic"]] <- tryCatch(
  { message("Fitting: m17_raw_main_effects_aperiodic")
    bam(update(base_formula_aperiodic,
               . ~ . + s(pow_slow_alpha_z, k = 10) + s(pow_fast_alpha_z, k = 10)),
        data = tensor_data %>% filter(!is.na(fooof_offset_z)),
        method = "fREML", discrete = use_discrete, nthreads = nthreads_to_use) },
  error = function(e) { warning("m17 failed: ", e$message); NULL }
)

# Phase model (sin + cos as linear additive terms)
if (all(c("phase_sin_z", "phase_cos_z") %in% names(df_model))) {
  model_list[["m18_phase_sincos"]] <- fit_v2_gamm(
    data        = df_model %>% filter(!is.na(phase_sin_z), !is.na(phase_cos_z)),
    metric_z    = NA_character_,
    model_name  = "m18_phase_sincos",
    use_aperiodic = FALSE,
    extra_terms = c("phase_sin_z", "phase_cos_z")
  )
  model_list[["m19_phase_sincos_aperiodic"]] <- fit_v2_gamm(
    data        = df_model %>% filter(!is.na(phase_sin_z), !is.na(phase_cos_z),
                                       !is.na(fooof_offset_z)),
    metric_z    = NA_character_,
    model_name  = "m19_phase_sincos_aperiodic",
    use_aperiodic = TRUE,
    extra_terms = c("phase_sin_z", "phase_cos_z")
  )
}

# Drop failed models
valid_idx  <- !map_lgl(model_list, is.null)
model_list <- model_list[valid_idx]

if (length(model_list) == 0) stop("No GAMMs were successfully fit.")
message("Successfully fit: ", paste(names(model_list), collapse = ", "))

# =========================================================
# SMOOTH EFFECT PLOTS  (univariate predictors only)
# =========================================================
smooth_plot_specs <- tribble(
  ~model_name,                   ~term,
  "m02_sf_logratio",             "sf_logratio_z",
  "m03_sf_balance",              "sf_balance_z",
  "m04_slow_frac",               "slow_alpha_frac_z",
  "m05_sf_logratio_aperiodic",   "sf_logratio_z",
  "m06_sf_balance_aperiodic",    "sf_balance_z",
  "m07_slow_frac_aperiodic",     "slow_alpha_frac_z",
  "m08_bi_pre",                  "bi_pre_z",
  "m09_cog_pre",                 "cog_pre_z",
  "m10_psi_cog",                 "psi_cog_z",
  "m11_delta_erd",               "delta_erd_z",
  "m12_bi_pre_aperiodic",        "bi_pre_z",
  "m13_psi_cog_aperiodic",       "psi_cog_z",
  "m14_delta_erd_aperiodic",     "delta_erd_z"
)

pwalk(smooth_plot_specs, function(model_name, term) {
  if (!model_name %in% names(model_list)) return(invisible(NULL))
  plot_smooth(
    model     = model_list[[model_name]],
    term      = term,
    file_name = paste0(model_name, "_", term, "_effect.png")
  )
})

# =========================================================
# INTERACTION SURFACE PLOTS (tensor models)
# =========================================================
tensor_surface_specs <- list(
  list(nm = "m15_tensor_slow_fast",
       v  = c("pow_slow_alpha_z", "pow_fast_alpha_z"),
       title = "Slow x Fast Alpha Tensor GAMM"),
  list(nm = "m16_tensor_slow_fast_aperiodic",
       v  = c("pow_slow_alpha_z", "pow_fast_alpha_z"),
       title = "Slow x Fast Alpha Tensor GAMM (Aperiodic Controlled)")
)

for (spec in tensor_surface_specs) {
  if (!spec$nm %in% names(model_list)) next
  png(file.path(out_dir, paste0(spec$nm, "_surface.png")),
      width = 1800, height = 1400, res = 200)
  vis.gam(model_list[[spec$nm]],
          view      = spec$v,
          plot.type = "contour",
          too.far   = 0.05,
          main      = spec$title)
  dev.off()
}

# =========================================================
# SAVE SUMMARIES + MODEL COMPARISON
# =========================================================
iwalk(model_list, function(mod, nm) {
  capture.output(summary(mod),
                 file = file.path(out_dir, paste0(nm, "_summary.txt")))
})

model_comparison <- imap_dfr(model_list, extract_model_info) %>%
  arrange(AIC)

write_csv(model_comparison, file.path(out_dir, "model_comparison_v2.csv"))
message("Model comparison saved.")

# =========================================================
# FITTED VALUES
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

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values_v2.csv"))

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

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values_v2.csv"))

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

message("v2 GAMM workflow complete. Outputs saved to: ", out_dir)