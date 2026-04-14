# =========================================================
# run_gamm_alpha_metrics_v2.R
# ---------------------------------------------------------
# V2 pooled GAMM workflow: raw slow-fast tensor interaction,
# aperiodic-controlled comparators, and V2 pre-stimulus
# interaction metrics (BI_pre, CoG_pre, psi_cog, DELTA_ERD,
# phase sin/cos).
# V 3.0.0
#
# V2.0.0 changes:
#   - Fixed typos: an.factor, use_dicrete, paf_cof_hz_z.
#   - Updated canonical column names throughout.
#   - plot_smooth: term identification uses x$term[1].
#   - Added V2 interaction metric models (m08-m19).
#   - Refactored to tibble model spec table + pmap.
#   - fitted_trial_df / subject_ga_df use any_of().
#
# V3.0.0 changes (mirrors v1 V3.0.0 fixes):
#   - safe_k(): adaptive k from actual unique values in fit
#     data, clamped to [k_min, k_max]. Fixes "fewer unique
#     covariate combinations" errors for discrete covariates
#     (laser power, fooof columns, etc.).
#   - build_v2_formula(): builds formula dynamically from
#     fit_data after per-metric NA removal. Guards against
#     single-level sex and cap_size (present in ALL v2
#     baselines). Conditionally includes fooof smooths only
#     when fooof columns are non-NA in fit_data.
#   - fit_v2_gamm(): NA removal (metric + fooof) happens
#     BEFORE formula construction so safe_k and level guards
#     operate on the real sample.
#   - fit_tensor_gamm(): dedicated helper for te() models;
#     adapts k per marginal dimension via safe_k.
#   - extra_smooths parameter on fit_v2_gamm() handles m17
#     (two additional s() terms) without separate bam() calls.
#   - select=TRUE in all bam() calls for automatic smooth
#     term selection via double penalty + fREML.
#   - tribble() replaced with tibble::tibble() in both
#     model_specs and smooth_plot_specs.
#   - age_z treated as linear fixed effect (not smooth) —
#     between-subject variable, not estimable as nonlinear
#     smooth at small N. Restore s(age_z) once N >> 50.
#   - Observed-vs-fitted plots: geom_jitter + opaque white
#     background + improved alpha for small N.
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
out_dir         <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/figures_v2"

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

# safe_k: upper bound on smooth complexity from unique values.
# k sets the maximum degrees of freedom; actual EDF is
# estimated from data via fREML penalisation.
safe_k <- function(x, k_max = 10, k_min = 3) {
  n_unique <- length(unique(na.omit(x)))
  max(k_min, min(k_max, n_unique - 1L))
}

# build_v2_formula: constructs the baseline formula from the
# actual fit_data. All v2 models include cap_size, so the
# level guard applies to every model (not just cap variants).
# Fooof smooths are included only when columns are present
# and non-NA in fit_data.
build_v2_formula <- function(data, use_aperiodic = FALSE) {
  k_laser <- safe_k(data$laser_power_z)
  k_trial <- safe_k(data$trial_index_z)

  smooth_terms <- c(
    sprintf("s(laser_power_z, k = %d)", k_laser),
    sprintf("s(trial_index_z, k = %d)", k_trial)
  )

  # Fooof smooths: only when aperiodic requested AND data present
  if (use_aperiodic) {
    if ("fooof_offset_z" %in% names(data) &&
        !all(is.na(data$fooof_offset_z))) {
      smooth_terms <- c(smooth_terms,
        sprintf("s(fooof_offset_z,   k = %d)", safe_k(data$fooof_offset_z)))
    }
    if ("fooof_exponent_z" %in% names(data) &&
        !all(is.na(data$fooof_exponent_z))) {
      smooth_terms <- c(smooth_terms,
        sprintf("s(fooof_exponent_z, k = %d)", safe_k(data$fooof_exponent_z)))
    }
  }

  # Fixed effects: guard against single-level factors
  fixed_terms <- "age_z"
  if (nlevels(data$sex)      > 1) fixed_terms <- c(fixed_terms, "sex")
  if (nlevels(data$cap_size) > 1) fixed_terms <- c(fixed_terms, "cap_size")

  rhs <- paste(
    c(smooth_terms, fixed_terms,
      "s(global_subjid, bs = 're')",
      "s(experiment_id,  bs = 're')"),
    collapse = " + "
  )

  as.formula(paste("pain_rating ~", rhs))
}

# extract_model_info: AIC/BIC/logLik/dev_expl summary row.
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

# plot_smooth: saves a partial effect plot for a single
# univariate smooth term. Uses x$term[1] for identification
# so it is safe for both s() and te() smooths.
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

required_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
  "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating",
  "pow_slow_alpha", "pow_fast_alpha",
  "sf_logratio", "sf_balance", "slow_alpha_frac"
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

# Phase decomposition: circular predictor → sin + cos components
if ("phase_slow_rad" %in% names(df)) {
  df <- df %>%
    mutate(
      phase_sin = sin(phase_slow_rad),
      phase_cos = cos(phase_slow_rad)
    )
}

# =========================================================
# QC FILTERING
# V2 requires the core V1 sub-band metrics globally because
# tensor and raw-power models need them on every row.
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
    !is.na(pow_slow_alpha),
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
    age_z             = zscore(age),
    laser_power_z     = zscore(laser_power),
    trial_index_z     = zscore(trial_index),
    pow_slow_alpha_z  = zscore(pow_slow_alpha),
    pow_fast_alpha_z  = zscore(pow_fast_alpha),
    sf_logratio_z     = zscore(sf_logratio),
    sf_balance_z      = zscore(sf_balance),
    slow_alpha_frac_z = zscore(slow_alpha_frac)
  )

# Optional columns: z-score only when present and not all-NA
optional_z <- c(
  "paf_cog_hz",
  "fooof_offset", "fooof_exponent",
  "bi_pre", "lr_pre", "cog_pre", "psi_cog", "delta_erd",
  "phase_sin", "phase_cos"
)

for (col in optional_z) {
  if (col %in% names(df_model) && !all(is.na(df_model[[col]]))) {
    df_model[[paste0(col, "_z")]] <- zscore(df_model[[col]])
  }
}

write_csv(df_model, file.path(out_dir, "alpha_pain_master_model_input_v2.csv"))

message("Model input rows: ", nrow(df_model))
message("Subjects: ",         n_distinct(df_model$subjid_uid))
message("Experiments: ",      n_distinct(df_model$experiment_id))
message("Unique laser_power values: ",
        length(unique(na.omit(df_model$laser_power))))

# =========================================================
# MODEL FITTERS
# ---------------------------------------------------------
# fit_v2_gamm: standard and aperiodic-controlled scalar
#   metric models. Per-metric AND fooof NAs are dropped
#   BEFORE formula construction so safe_k and level guards
#   operate on the real sample. extra_smooths adds additional
#   s() terms with adaptive k (used for m17).
#
# fit_tensor_gamm: te() interaction surface models. Cannot
#   use the metric_z path because te() requires two variable
#   names; k is adapted per marginal dimension.
# =========================================================
fit_v2_gamm <- function(data, metric_z, model_name,
                         use_aperiodic = FALSE,
                         extra_terms   = NULL,
                         extra_smooths = NULL) {

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

  # For aperiodic models, also drop fooof NA rows
  if (use_aperiodic) {
    if ("fooof_offset_z" %in% names(fit_data))
      fit_data <- fit_data %>% filter(!is.na(fooof_offset_z))
    if ("fooof_exponent_z" %in% names(fit_data))
      fit_data <- fit_data %>% filter(!is.na(fooof_exponent_z))
  }

  if (nrow(fit_data) == 0) {
    warning("Skipping ", model_name, ": no rows after NA removal.")
    return(NULL)
  }

  # Build baseline formula from actual fit_data
  base_formula <- build_v2_formula(fit_data, use_aperiodic)

  # Add metric smooth with adaptive k
  formula_to_fit <- if (!is.na(metric_z)) {
    k_metric <- safe_k(fit_data[[metric_z]])
    update(base_formula,
           as.formula(sprintf(". ~ . + s(%s, k = %d)", metric_z, k_metric)))
  } else {
    base_formula
  }

  # Extra smooths with adaptive k (e.g., two raw power terms in m17)
  if (!is.null(extra_smooths)) {
    for (col in extra_smooths) {
      k_extra <- safe_k(fit_data[[col]])
      formula_to_fit <- update(
        formula_to_fit,
        as.formula(sprintf(". ~ . + s(%s, k = %d)", col, k_extra))
      )
    }
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
    bam(formula  = formula_to_fit,
        data     = fit_data,
        method   = "fREML",
        select   = TRUE,
        discrete = use_discrete,
        nthreads = nthreads_to_use),
    error = function(e) {
      warning("Model failed [", model_name, "]: ", e$message)
      NULL
    }
  )
}

fit_tensor_gamm <- function(data, model_name, use_aperiodic = FALSE) {
  fit_data <- data

  if (use_aperiodic) {
    if ("fooof_offset_z" %in% names(fit_data))
      fit_data <- fit_data %>% filter(!is.na(fooof_offset_z))
    if ("fooof_exponent_z" %in% names(fit_data))
      fit_data <- fit_data %>% filter(!is.na(fooof_exponent_z))
  }

  if (nrow(fit_data) == 0) {
    warning("Skipping ", model_name, ": no rows after NA removal.")
    return(NULL)
  }

  base_formula <- build_v2_formula(fit_data, use_aperiodic)
  k_slow <- safe_k(fit_data$pow_slow_alpha_z)
  k_fast <- safe_k(fit_data$pow_fast_alpha_z)

  formula_to_fit <- update(
    base_formula,
    as.formula(sprintf(
      ". ~ . + te(pow_slow_alpha_z, pow_fast_alpha_z, k = c(%d, %d))",
      k_slow, k_fast
    ))
  )

  message("Fitting: ", model_name, "  (n=", nrow(fit_data), ")")

  tryCatch(
    bam(formula  = formula_to_fit,
        data     = fit_data,
        method   = "fREML",
        select   = TRUE,
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
# Tensor and phase models are handled outside pmap below.
# =========================================================
model_specs <- tibble::tibble(
  model_name = c(
    # Baselines
    "m00_baseline",               "m01_baseline_aperiodic",
    # V1 sub-band metrics: standard
    "m02_sf_logratio",            "m03_sf_balance",          "m04_slow_frac",
    # V1 sub-band metrics: aperiodic-controlled
    "m05_sf_logratio_aperiodic",  "m06_sf_balance_aperiodic","m07_slow_frac_aperiodic",
    # V2 pre-stimulus interaction metrics: standard
    "m08_bi_pre",                 "m09_cog_pre",
    "m10_psi_cog",                "m11_delta_erd",
    # V2 pre-stimulus interaction metrics: aperiodic-controlled
    "m12_bi_pre_aperiodic",       "m13_psi_cog_aperiodic",   "m14_delta_erd_aperiodic"
  ),
  metric_z = c(
    NA_character_,       NA_character_,
    "sf_logratio_z",     "sf_balance_z",      "slow_alpha_frac_z",
    "sf_logratio_z",     "sf_balance_z",      "slow_alpha_frac_z",
    "bi_pre_z",          "cog_pre_z",
    "psi_cog_z",         "delta_erd_z",
    "bi_pre_z",          "psi_cog_z",         "delta_erd_z"
  ),
  use_aperiodic = c(
    FALSE, TRUE,
    FALSE, FALSE, FALSE,
    TRUE,  TRUE,  TRUE,
    FALSE, FALSE, FALSE, FALSE,
    TRUE,  TRUE,  TRUE
  )
)

# =========================================================
# FIT MODELS
# =========================================================
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

# Tensor interaction models (te() requires dedicated helper)
model_list[["m15_tensor_slow_fast"]] <-
  fit_tensor_gamm(df_model, "m15_tensor_slow_fast", use_aperiodic = FALSE)

model_list[["m16_tensor_slow_fast_aperiodic"]] <-
  fit_tensor_gamm(df_model, "m16_tensor_slow_fast_aperiodic", use_aperiodic = TRUE)

# Raw main effects + aperiodic: two s() terms added via extra_smooths
model_list[["m17_raw_main_effects_aperiodic"]] <- fit_v2_gamm(
  data          = df_model,
  metric_z      = NA_character_,
  model_name    = "m17_raw_main_effects_aperiodic",
  use_aperiodic = TRUE,
  extra_smooths = c("pow_slow_alpha_z", "pow_fast_alpha_z")
)

# Phase sin + cos as linear additive terms
if (all(c("phase_sin_z", "phase_cos_z") %in% names(df_model))) {
  phase_data <- df_model %>% filter(!is.na(phase_sin_z), !is.na(phase_cos_z))

  model_list[["m18_phase_sincos"]] <- fit_v2_gamm(
    data          = phase_data,
    metric_z      = NA_character_,
    model_name    = "m18_phase_sincos",
    use_aperiodic = FALSE,
    extra_terms   = c("phase_sin_z", "phase_cos_z")
  )

  model_list[["m19_phase_sincos_aperiodic"]] <- fit_v2_gamm(
    data          = phase_data,
    metric_z      = NA_character_,
    model_name    = "m19_phase_sincos_aperiodic",
    use_aperiodic = TRUE,
    extra_terms   = c("phase_sin_z", "phase_cos_z")
  )
}

# Drop failed / skipped models
valid_idx  <- !map_lgl(model_list, is.null)
model_list <- model_list[valid_idx]

if (length(model_list) == 0) stop("No GAMMs were successfully fit.")
message("Successfully fit: ", paste(names(model_list), collapse = ", "))

# =========================================================
# SMOOTH EFFECT PLOTS (univariate predictors only)
# =========================================================
smooth_plot_specs <- tibble::tibble(
  model_name = c(
    "m02_sf_logratio",            "m03_sf_balance",          "m04_slow_frac",
    "m05_sf_logratio_aperiodic",  "m06_sf_balance_aperiodic","m07_slow_frac_aperiodic",
    "m08_bi_pre",                 "m09_cog_pre",
    "m10_psi_cog",                "m11_delta_erd",
    "m12_bi_pre_aperiodic",       "m13_psi_cog_aperiodic",   "m14_delta_erd_aperiodic"
  ),
  term = c(
    "sf_logratio_z",     "sf_balance_z",      "slow_alpha_frac_z",
    "sf_logratio_z",     "sf_balance_z",      "slow_alpha_frac_z",
    "bi_pre_z",          "cog_pre_z",
    "psi_cog_z",         "delta_erd_z",
    "bi_pre_z",          "psi_cog_z",         "delta_erd_z"
  )
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
  list(nm    = "m15_tensor_slow_fast",
       v     = c("pow_slow_alpha_z", "pow_fast_alpha_z"),
       title = "Slow x Fast Alpha Tensor GAMM"),
  list(nm    = "m16_tensor_slow_fast_aperiodic",
       v     = c("pow_slow_alpha_z", "pow_fast_alpha_z"),
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
  capture.output(
    summary(mod),
    file = file.path(out_dir, paste0(nm, "_summary.txt"))
  )
})

model_comparison <- imap_dfr(model_list, extract_model_info) %>%
  arrange(AIC)

write_csv(model_comparison, file.path(out_dir, "model_comparison_v2.csv"))
message("Model comparison saved.")

# =========================================================
# FITTED VALUES
# Spectral columns pass through via any_of() automatically.
# Fitted/residual values are aligned by original row index;
# rows not in a model's fit sample receive NA.
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

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values_v2.csv"))

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

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values_v2.csv"))

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
    geom_jitter(width = 0, height = 0.15, alpha = 0.5, size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, colour = "steelblue") +
    labs(
      title = paste("Observed vs Fitted:", nm),
      x     = "Fitted pain rating",
      y     = "Observed pain rating (jittered)"
    ) +
    theme_bw(base_size = 13)

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

message("v2 GAMM workflow complete. Outputs saved to: ", out_dir)