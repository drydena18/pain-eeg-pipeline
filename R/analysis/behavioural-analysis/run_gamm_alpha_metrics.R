# =========================================================
# run_gamm_alpha_metrics.R
# ---------------------------------------------------------
# First-pass GAMM workflow for pooled trial-level alpha-pain
# analyses across all participants / experiments.
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
out_dir   <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_fooof_r2 <- 0.80
max_fooof_error <- Inf
use_discrete <- TRUE
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

zscore <- function(x) as.numeric(scale(x))
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
# =========================================================
df_model <- df_model %>%
  mutate(
    age_z         = zscore(age),
    laser_power_z = zscore(laser_power),
    trial_index_z = zscore(trial_index)
  )

metric_candidates <- c(
  "paf_cog_hz",
  "pow_slow",
  "pow_fast",
  "pow_alpha",
  "rel_slow",
  "rel_fast",
  "sf_ratio",
  "sf_logratio",
  "sf_balance",
  "slow_frac",
  "fooof_offset",
  "fooof_exponent"
)

for (metric in metric_candidates) {
  if (metric %in% names(df_model)) {
    df_model[[paste0(metric, "_z")]] <- zscore(df_model[[metric]])
  }
}

write_csv(df_model, file.path(out_dir, "alpha_pain_master_model_input.csv"))

# =========================================================
# BASELINE FORMULAS
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
# =========================================================
fit_metric_gamm <- function(data, metric_z = NULL, model_name, include_cap = FALSE) {
  base_formula <- if (include_cap) baseline_cap_formula else baseline_formula

  if (is.null(metric_z)) {
    formula_to_fit <- base_formula
  } else {
    if (!metric_z %in% names(data)) {
      warning("Skipping ", model_name, ": metric not found -> ", metric_z)
      return(NULL)
    }

    formula_to_fit <- update(
      base_formula,
      paste(". ~ . + s(", metric_z, ", k = 10)", sep = "")
    )
  }

  message("Fitting model: ", model_name)

  bam(
    formula = formula_to_fit,
    data = data,
    method = "fREML",
    discrete = use_discrete,
    nthreads = nthreads_to_use
  )
}

# =========================================================
# MODEL SET
# =========================================================
model_specs <- tribble(
  ~model_name,              ~metric_z,         ~include_cap,
  "m00_baseline",           NA_character_,     FALSE,
  "m01_baseline_cap",       NA_character_,     TRUE,
  "m02_paf",                "paf_cog_hz_z",    FALSE,
  "m03_pow_slow",           "pow_slow_z",      FALSE,
  "m04_pow_fast",           "pow_fast_z",      FALSE,
  "m05_sf_ratio",           "sf_ratio_z",      FALSE,
  "m06_sf_logratio",        "sf_logratio_z",   FALSE,
  "m07_sf_balance",         "sf_balance_z",    FALSE,
  "m08_slow_frac",          "slow_frac_z",     FALSE,
  "m09_sf_logratio_cap",    "sf_logratio_z",   TRUE,
  "m10_sf_balance_cap",     "sf_balance_z",    TRUE,
  "m11_slow_frac_cap",      "slow_frac_z",     TRUE
)

# =========================================================
# FIT MODELS
# =========================================================
model_list <- pmap(
  model_specs,
  function(model_name, metric_z, include_cap) {
    fit_metric_gamm(
      data = df_model,
      metric_z = metric_z,
      model_name = model_name,
      include_cap = include_cap
    )
  }
)

names(model_list) <- model_specs$model_name
valid_idx <- !map_lgl(model_list, is.null)
model_list <- model_list[valid_idx]
model_specs <- model_specs[valid_idx, ]

if (length(model_list) == 0) {
  stop("No GAMMs were successfully fit.")
}

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
# MODEL COMPARISON
# =========================================================
extract_model_info <- function(mod, nm) {
  tibble(
    model_name = nm,
    AIC = AIC(mod),
    BIC = BIC(mod),
    logLik = as.numeric(logLik(mod)),
    dev_expl = summary(mod)$dev.expl,
    n = nobs(mod)
  )
}

model_comparison <- imap_dfr(model_list, extract_model_info) %>%
  arrange(AIC)

write_csv(model_comparison, file.path(out_dir, "model_comparison.csv"))

# =========================================================
# FITTED VALUES
# =========================================================
fitted_trial_df <- df_model %>%
  select(
    experiment_name,
    experiment_id,
    subjid,
    subjid_uid,
    global_subjid,
    age,
    sex,
    cap_size,
    trial,
    trial_index,
    laser_power,
    pain_rating,
    any_of(c(
      "paf_cog_hz",
      "pow_slow",
      "pow_fast",
      "sf_ratio",
      "sf_logratio",
      "sf_balance",
      "slow_frac"
    ))
  )

for (nm in names(model_list)) {
  fitted_trial_df[[paste0(nm, "_fitted")]] <- fitted(model_list[[nm]])
  fitted_trial_df[[paste0(nm, "_resid")]]  <- residuals(model_list[[nm]])
}

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values.csv"))

subject_ga_df <- fitted_trial_df %>%
  group_by(
    experiment_name,
    experiment_id,
    subjid,
    subjid_uid,
    global_subjid,
    age,
    sex,
    cap_size
  ) %>%
  summarise(
    n_trials = n(),
    mean_pain_rating = mean(pain_rating, na.rm = TRUE),
    mean_laser_power = mean(laser_power, na.rm = TRUE),
    mean_paf_cog_hz = mean(paf_cog_hz, na.rm = TRUE),
    mean_pow_slow = mean(pow_slow, na.rm = TRUE),
    mean_pow_fast = mean(pow_fast, na.rm = TRUE),
    mean_sf_ratio = mean(sf_ratio, na.rm = TRUE),
    mean_sf_logratio = mean(sf_logratio, na.rm = TRUE),
    mean_sf_balance = mean(sf_balance, na.rm = TRUE),
    mean_slow_frac = mean(slow_frac, na.rm = TRUE),
    across(
      ends_with("_fitted"),
      ~ mean(.x, na.rm = TRUE),
      .names = "ga_{.col}"
    ),
    .groups = "drop"
  )

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values.csv"))

# =========================================================
# DIAGNOSTICS
# =========================================================
iwalk(model_list, function(mod, nm) {
  png(file.path(out_dir, paste0(nm, "_diagnostics.png")), width = 1800, height = 1400, res = 180)
  gam.check(mod)
  dev.off()
})

for (nm in names(model_list)) {
  p <- ggplot(
    fitted_trial_df,
    aes(
      x = .data[[paste0(nm, "_fitted")]],
      y = pain_rating
    )
  ) +
    geom_point(alpha = 0.15) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste("Observed vs Fitted:", nm),
      x = "Fitted pain",
      y = "Observed pain"
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(out_dir, paste0(nm, "_observed_vs_fitted.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

iwalk(model_list, function(mod, nm) {
  saveRDS(mod, file.path(out_dir, paste0(nm, ".rds")))
})

message("Done.")
message("Outputs saved to: ", out_dir)