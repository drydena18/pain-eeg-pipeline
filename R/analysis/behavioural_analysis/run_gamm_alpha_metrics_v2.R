# =========================================================
# run_gamm_alpha_metrics_v2.R
# ---------------------------------------------------------
# v2 pooled GAMM workflow — pre-stimulus alpha interaction metrics
#
# Model family:
#   m00  baseline (no spectral predictor)
#   m01  baseline + aperiodic controls
#   m02  LR_pre  (pre-stim log ratio; replaces full-epoch sf_logratio)
#   m03  BI_pre  (pre-stim balance index; replaces full-epoch sf_balance)
#   m04  CoG_pre (pre-stim centre of gravity)
#   m05  delta_erd (ΔERD asymmetry index)
#   m06  psi_cog (BI_pre x CoG_pre interaction term)
#   m07  tensor: pre-stim slow x fast power
#   m08  LR_pre + aperiodic controls
#   m09  BI_pre + aperiodic controls
#   m10  delta_erd + aperiodic controls
#   m11  slow-alpha phase: sin(phase) + cos(phase) main effects
#   m12  BI_pre x phase interaction (BI_pre + sin + cos + BI_pre:sin + BI_pre:cos)
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
data_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/alpha_pain_master.csv"
out_dir   <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_v2"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_fooof_r2    <- 0.80
max_fooof_error <- Inf

use_discrete    <- TRUE          # BUG FIX: was use_dicrete (typo; bam() never received it)
nthreads_to_use <- 4

# =========================================================
# HELPERS
# =========================================================
clean_names_local <- function(x) {
  x |>
    str_trim() |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("\\^2", "r2") |>
    str_replace_all("\\.", "_")
}

zscore      <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))  # BUG FIX: was an.factor (undefined)

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

# =========================================================
# LOAD DATA
# =========================================================
df <- read_csv(data_file, show_col_types = FALSE)
names(df) <- clean_names_local(names(df))

# Required columns (updated to pre-stim metric names)
required_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
  "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating",
  "pow_slow_alpha", "pow_fast_alpha",
  "bi_pre", "lr_pre", "cog_pre", "psi_cog", "delta_erd"
)

missing_required <- setdiff(required_cols, names(df))
if (length(missing_required) > 0) {
  stop("Missing required columns in alpha_pain_master.csv: ",
       paste(missing_required, collapse = ", "))
}

df <- df |>
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
# QC FILTERING
# =========================================================
df_model <- df |>
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
    !is.na(bi_pre),
    !is.na(lr_pre),
    !is.na(cog_pre),
    !is.na(psi_cog),
    !is.na(delta_erd)
  )

# Drop trials flagged as having unstable pre-stim power
if ("p5_flag" %in% names(df_model)) {
  n_before <- nrow(df_model)
  df_model <- df_model |> filter(is.na(p5_flag) | p5_flag == 0)
  message(sprintf("p5_flag: removed %d unstable-power trials.", n_before - nrow(df_model)))
}

if ("fooof_r2" %in% names(df_model)) {
  df_model <- df_model |>
    filter(is.na(fooof_r2) | fooof_r2 >= min_fooof_r2)
}

if ("fooof_error" %in% names(df_model) && is.finite(max_fooof_error)) {
  df_model <- df_model |>
    filter(is.na(fooof_error) | fooof_error <= max_fooof_error)
}

message(sprintf("Rows after QC filtering: %d", nrow(df_model)))

# =========================================================
# SCALE CONTINUOUS PREDICTORS
# =========================================================
df_model <- df_model |>
  mutate(
    age_z         = zscore(age),
    laser_power_z = zscore(laser_power),
    trial_index_z = zscore(trial_index),
    pow_slow_z    = zscore(pow_slow_alpha),
    pow_fast_z    = zscore(pow_fast_alpha),
    bi_pre_z      = zscore(bi_pre),
    lr_pre_z      = zscore(lr_pre),
    cog_pre_z     = zscore(cog_pre),
    psi_cog_z     = zscore(psi_cog),
    delta_erd_z   = zscore(delta_erd)
  )

# Optional predictors that may not be present in every run
if ("paf_cog_hz" %in% names(df_model)) {
  df_model$paf_cog_hz_z <- zscore(df_model$paf_cog_hz)  # BUG FIX: was paf_cof_hz_z (typo)
}
if ("fooof_offset" %in% names(df_model)) {
  df_model$fooof_offset_z   <- zscore(df_model$fooof_offset)
}
if ("fooof_exponent" %in% names(df_model)) {
  df_model$fooof_exponent_z <- zscore(df_model$fooof_exponent)
}

# Phase metric: derive sin and cos components if available
has_phase <- "phase_slow_rad" %in% names(df_model) &&
             !all(is.na(df_model$phase_slow_rad))

if (has_phase) {
  df_model <- df_model |>
    mutate(
      phase_sin = sin(phase_slow_rad),
      phase_cos = cos(phase_slow_rad)
    )
  message("Phase metric available: sin/cos columns created.")
} else {
  message("phase_slow_rad not found or all-NA; phase models (m11, m12) will be skipped.")
}

write_csv(df_model, file.path(out_dir, "alpha_pain_master_model_input_v2.csv"))

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

has_fooof_controls <- all(c("fooof_offset_z", "fooof_exponent_z") %in% names(df_model))

if (has_fooof_controls) {
  base_formula_aperiodic <- update(base_formula, . ~ . +
    s(fooof_offset_z, k = 10) +
    s(fooof_exponent_z, k = 10))
} else {
  message("FOOOF aperiodic controls absent; aperiodic models (m01, m08-m10) will be identical to base.")
  base_formula_aperiodic <- base_formula
}

# =========================================================
# FIT MODELS
# =========================================================
fit_bam <- function(formula, label) {
  message("Fitting ", label, " ...")
  bam(
    formula  = formula,
    data     = df_model,
    method   = "fREML",
    discrete = use_discrete,
    nthreads = nthreads_to_use
  )
}

# m00: no spectral predictor
m00_baseline <- fit_bam(base_formula, "m00_baseline")

# m01: aperiodic controls only
m01_baseline_aperiodic <- fit_bam(base_formula_aperiodic, "m01_baseline_aperiodic")

# m02: pre-stim log ratio (replaces full-epoch sf_logratio)
m02_lr_pre <- fit_bam(
  update(base_formula, . ~ . + s(lr_pre_z, k = 10)),
  "m02_lr_pre"
)

# m03: pre-stim balance index (replaces full-epoch sf_balance)
m03_bi_pre <- fit_bam(
  update(base_formula, . ~ . + s(bi_pre_z, k = 10)),
  "m03_bi_pre"
)

# m04: pre-stim centre of gravity
m04_cog_pre <- fit_bam(
  update(base_formula, . ~ . + s(cog_pre_z, k = 10)),
  "m04_cog_pre"
)

# m05: delta ERD asymmetry index
m05_delta_erd <- fit_bam(
  update(base_formula, . ~ . + s(delta_erd_z, k = 10)),
  "m05_delta_erd"
)

# m06: BI_pre x CoG_pre interaction (psi_cog)
# Main effects of bi_pre and cog_pre are included alongside their product.
m06_psi_cog <- fit_bam(
  update(base_formula, . ~ . + s(bi_pre_z, k = 10) + s(cog_pre_z, k = 10) + s(psi_cog_z, k = 10)),
  "m06_psi_cog"
)

# m07: tensor of pre-stim slow x fast power (non-linear interaction surface)
m07_tensor_prestim <- fit_bam(
  update(base_formula, . ~ . + te(pow_slow_z, pow_fast_z, k = c(10, 10))),
  "m07_tensor_prestim"
)

# m08–m10: aperiodic-controlled replication of m02, m03, m05
m08_lr_pre_aperiodic <- fit_bam(
  update(base_formula_aperiodic, . ~ . + s(lr_pre_z, k = 10)),
  "m08_lr_pre_aperiodic"
)

m09_bi_pre_aperiodic <- fit_bam(
  update(base_formula_aperiodic, . ~ . + s(bi_pre_z, k = 10)),
  "m09_bi_pre_aperiodic"
)

m10_delta_erd_aperiodic <- fit_bam(
  update(base_formula_aperiodic, . ~ . + s(delta_erd_z, k = 10)),
  "m10_delta_erd_aperiodic"
)

model_list <- list(
  m00_baseline            = m00_baseline,
  m01_baseline_aperiodic  = m01_baseline_aperiodic,
  m02_lr_pre              = m02_lr_pre,
  m03_bi_pre              = m03_bi_pre,
  m04_cog_pre             = m04_cog_pre,
  m05_delta_erd           = m05_delta_erd,
  m06_psi_cog             = m06_psi_cog,
  m07_tensor_prestim      = m07_tensor_prestim,
  m08_lr_pre_aperiodic    = m08_lr_pre_aperiodic,
  m09_bi_pre_aperiodic    = m09_bi_pre_aperiodic,
  m10_delta_erd_aperiodic = m10_delta_erd_aperiodic
)

# m11 and m12: phase models (fitted only if phase_slow_rad is available)
if (has_phase) {
  # m11: slow-alpha phase main effects (sin + cos linear decomposition)
  m11_phase_slow <- fit_bam(
    update(base_formula, . ~ . + phase_sin + phase_cos),
    "m11_phase_slow"
  )

  # m12: BI_pre x phase interaction
  # Includes main effects of BI_pre, sin(phase), cos(phase) plus their products.
  # This tests whether phase predicts pain more strongly on high-BI_pre trials.
  m12_bi_phase_interact <- fit_bam(
    update(base_formula, . ~ .
      + s(bi_pre_z, k = 10)
      + phase_sin
      + phase_cos
      + s(bi_pre_z, by = phase_sin, k = 10)
      + s(bi_pre_z, by = phase_cos, k = 10)),
    "m12_bi_phase_interact"
  )

  model_list[["m11_phase_slow"]]       <- m11_phase_slow
  model_list[["m12_bi_phase_interact"]] <- m12_bi_phase_interact
} else {
  message("Skipping m11 and m12: phase_slow_rad not available.")
}

# =========================================================
# SMOOTH EFFECT PLOTS
# =========================================================
plot_smooth <- function(model, term, file_name) {
  # Find the smooth index by matching term name
  smooth_idx <- which(
    sapply(model$smooth, function(s) {
      length(s$term) == 1 && s$term == term
    })
  )
  if (length(smooth_idx) == 0) {
    message("plot_smooth: term '", term, "' not found in model; skipping.")
    return(invisible(NULL))
  }

  png(
    filename = file.path(out_dir, file_name),
    width    = 1600,
    height   = 1200,
    res      = 200
  )
  plot(
    model,
    select    = smooth_idx[1],
    shade     = TRUE,
    shade.col = "lightblue",
    main      = paste("Effect of", term, "on Pain Rating"),
    xlab      = term,
    ylab      = "Partial Effect on Pain"
  )
  abline(h = 0, lty = 2)
  dev.off()
}

plot_smooth(m02_lr_pre,              "lr_pre_z",    "m02_lr_pre_effect.png")
plot_smooth(m03_bi_pre,              "bi_pre_z",    "m03_bi_pre_effect.png")
plot_smooth(m04_cog_pre,             "cog_pre_z",   "m04_cog_pre_effect.png")
plot_smooth(m05_delta_erd,           "delta_erd_z", "m05_delta_erd_effect.png")
plot_smooth(m08_lr_pre_aperiodic,    "lr_pre_z",    "m08_lr_pre_aperiodic_effect.png")
plot_smooth(m09_bi_pre_aperiodic,    "bi_pre_z",    "m09_bi_pre_aperiodic_effect.png")
plot_smooth(m10_delta_erd_aperiodic, "delta_erd_z", "m10_delta_erd_aperiodic_effect.png")

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
model_comparison <- imap_dfr(model_list, extract_model_info) |>
  arrange(AIC)

write_csv(model_comparison, file.path(out_dir, "model_comparison_v2.csv"))
message("Model comparison saved.")

# =========================================================
# FITTED VALUES + RESIDUALS
# =========================================================
select_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
  "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating"
)
optional_cols <- c(
  "paf_cog_hz", "pow_slow_alpha", "pow_fast_alpha",
  "bi_pre", "lr_pre", "cog_pre", "psi_cog", "delta_erd",
  "fooof_offset", "fooof_exponent",
  "phase_slow_rad", "phase_sin", "phase_cos"
)

fitted_trial_df <- df_model |>
  select(all_of(select_cols), any_of(optional_cols))

for (nm in names(model_list)) {
  fitted_trial_df[[paste0(nm, "_fitted")]] <- fitted(model_list[[nm]])
  fitted_trial_df[[paste0(nm, "_resid")]]  <- residuals(model_list[[nm]])
}

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values_v2.csv"))

# =========================================================
# SUBJECT GA SUMMARIES
# =========================================================
# Columns to average across trials per subject
ga_numeric_cols <- intersect(
  c("pain_rating", "laser_power",
    "pow_slow_alpha", "pow_fast_alpha",
    "bi_pre", "lr_pre", "cog_pre", "psi_cog", "delta_erd",
    "fooof_offset", "fooof_exponent"),
  names(fitted_trial_df)
)

subject_ga_df <- fitted_trial_df |>
  group_by(
    experiment_name, experiment_id,
    subjid, subjid_uid, global_subjid,
    age, sex, cap_size
  ) |>
  summarise(
    n_trials = n(),
    across(all_of(ga_numeric_cols), ~ mean(.x, na.rm = TRUE)),
    across(ends_with("_fitted"), ~ mean(.x, na.rm = TRUE), .names = "ga_{.col}"),
    .groups = "drop"
  )

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values_v2.csv"))

# =========================================================
# DIAGNOSTICS
# =========================================================
iwalk(model_list, function(mod, nm) {
  png(file.path(out_dir, paste0(nm, "_diagnostics.png")),
      width = 1800, height = 1400, res = 180)
  gam.check(mod)
  dev.off()
})

# =========================================================
# OBSERVED VS FITTED PLOTS
# =========================================================
for (nm in names(model_list)) {
  fitted_col <- paste0(nm, "_fitted")
  if (!fitted_col %in% names(fitted_trial_df)) next

  p <- ggplot(
    fitted_trial_df,
    aes(x = .data[[fitted_col]], y = pain_rating)
  ) +
    geom_point(alpha = 0.15) +
    geom_smooth(method = "lm", se = FALSE, colour = "steelblue") +
    labs(
      title = paste("Observed vs Fitted:", nm),
      x     = "Fitted pain rating",
      y     = "Observed pain rating"
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(out_dir, paste0(nm, "_observed_vs_fitted.png")),
    plot     = p,
    width    = 7,
    height   = 5,
    dpi      = 300
  )
}

# =========================================================
# INTERACTION SURFACE PLOTS
# =========================================================
png(file.path(out_dir, "m07_tensor_prestim_surface.png"),
    width = 1800, height = 1400, res = 200)
vis.gam(
  m07_tensor_prestim,
  view      = c("pow_slow_z", "pow_fast_z"),
  plot.type = "contour",
  too.far   = 0.05,
  main      = "Pre-stim Slow x Fast Alpha Tensor (GAMM)"
)
dev.off()

# =========================================================
# SAVE MODEL OBJECTS
# =========================================================
iwalk(model_list, function(mod, nm) {
  saveRDS(mod, file.path(out_dir, paste0(nm, ".rds")))
})

message("v2 GAMM workflow complete.  Outputs -> ", out_dir)