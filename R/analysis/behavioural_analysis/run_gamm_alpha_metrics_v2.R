# =========================================================
# run_gamm_alpha_metrics_v2.R
# ---------------------------------------------------------
# v2 pooled GAMM workflow for trial-level alpha-pain analyses.
#
# Focus:
#   1. Derived alpha interaction metrics
#   2. Raw slow-fast tensor interaction
#   3. Aperiodic-controlled comparator models
#   4. Save model summaries, comparisons, fitted values,
#      and interaction surface plots.
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
out_dir <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_v2"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_fooof_r2 <- 0.80
max_fooof_error <- Inf

use_dicrete = TRUE
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
safe_factor <- function(x) an.factor(as.character(x))

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

# =========================================================
# LOAD DATA
# =========================================================
df <- read_csv(data_file, show_col_types = FALSE)
names(df) <- clean_names_local(names(df))

required_cols <- c(
    "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
    "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating",
    "pow_slow", "pow_fast", "sf_logratio", "sf_balance", "slow_frac"
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
    !is.na(cap_size),
    !is.na(pow_slow),
    !is.na(pow_fast),
    !is.na(sf_logratio),
    !is.na(sf_balance),
    !is.na(slow_frac)
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
    pow_slow_z      = zscore(pow_slow),
    pow_fast_z      = zscore(pow_fast),
    sf_logratio_z   = zscore(sf_logratio),
    sf_balance_z    = zscore(sf_balance),
    slow_frac_z     = zscore(slow_frac)
  )

if ("paf_cog_hz" %in% names(df_model)) {
    df_model$paf_cof_hz_z <- zscore(df_model$paf_cog_hz)
}
if ("fooof_offset" %in% names(df_model)) {
    df_model$fooof_offset_z <- zscore(df_model$fooof_offset)
}
if ("fooof_exponent" %in% names(df_model)) {
    df_model$fooof_exponent_z <- zscore(df_model$fooof_exponent)
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
# FIT MODELS
# =========================================================
message("Fitting m00_baseline")
m00_baseline <- bam(
    formula = base_formula,
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m01_baseline_aperiodic")
m01_baseline_aperiodic <- bam(
    formula = base_formula_aperiodic,
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m02_sf_logratio")
m02_sf_logratio <- bam(
    formula = update(base_formula, . ~ . + s(sf_logratio_z, k = 10)),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m03_sf_balance")
m03_sf_balance <- bam(
    formula = update(base_formula, . ~ . + s(sf_balance_z, k = 10)),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m04_slow_frac")
m04_slow_frac <- bam(
    formula = update(base_formula, . ~ . + s(slow_frac_z, k = 10)),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m05_tensor_slow_fast")
m05_tensor_slow_fast <- bam(
    formula = update(base_formula, . ~ . + te(pow_slow_z, pow_fast_z, k = c(10, 10))),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m06_sf_logratio_aperiodic")
m06_sf_logratio_aperiodic <- bam(
    formula = update(base_formula_aperiodic, . ~ . + s(sf_logratio_z, k = 10)),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m07_sf_balance_aperiodic")
m07_sf_balance_aperiodic <- bam(
    formula = update(base_formula_aperiodic, . ~ . + s(sf_balance_z, k = 10)),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m08_slow_frac_aperiodic")
m08_slow_frac_aperiodic <- bam(
    formula = update(base_formula_aperiodic, . ~ . + s(slow_frac_z, k = 10)),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m09_tensor_slow_fast_aperiodic")
m09_tensor_slow_fast_aperiodic <- bam(
    formula = update(base_formula_aperiodic, . ~ . + te(pow_slow_z, pow_fast_z, k = c(10, 10))),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use
)

message("Fitting m10_raw_main_effects_aperiodic")
m10_raw_main_effects_aperiodic <- bam(
    formula = update(base_formula_aperiodic, . ~ . + s(pow_slow_z, k = 10) + s(pow_fast_z, k = 10)),
    data = df_model,
    method = "fREML",
    discrete = use_dicrete,
    nthreads = nthreads_to_use 
)

model_list <- list(
    m00_baseline = m00_baseline,
    m01_baseline_aperiodic = m01_baseline_aperiodic,
    m02_sf_logratio = m02_sf_logratio,
    m03_sf_balance = m03_sf_balance,
    m04_slow_frac = m04_slow_frac,
    m05_tensor_slow_fast = m05_tensor_slow_fast,
    m06_sf_logratio_aperiodic = m06_sf_logratio_aperiodic,
    m07_sf_balance_aperiodic = m07_sf_balance_aperiodic,
    m08_slow_frac_aperiodic = m08_slow_frac_aperiodic,
    m09_tensor_slow_fast_aperiodic = m09_tensor_slow_fast_aperiodic,
    m10_raw_main_effects_aperiodic = m10_raw_main_effects_aperiodic
)

# =========================================================
# SMOOTH EFFECT PLOTS
# =========================================================
plot_smooth <- function(model, term, file_name) {

    png(
        filename = file.path(out_dir, file_name),
        width = 1600,
        height = 1200,
        res = 200
    )

    plot(
        model,
        select = which(sapply(model$smooth, function(x) x$term) == term),
        shade = TRUE,
        shade.col = "lightblue",
        main = paste("Effect of", term, "on Pain Rating."),
        xlab = term,
        ylab = "Partial Effect on Pain"
    )

    abline(h = 0, lty = 2)

    dev.off()
}

# Scalar alpha metrics
plot_smooth(m02_sf_logratio, "sf_logratio_z", "m02_sf_logratio_effect.png")
plot_smooth(m03_sf_balance, "sf_balance_z", "m03_sf_balance_effect.png")
plot_smooth(m04_slow_frac, "slow_frac_z", "m04_slow_frac_effect.png")

# Aperiodic-controlled versions
plot_smooth(m06_sf_logratio_aperiodic, "sf_logratio_z", "m06_sf_logratio_aperiodic_effect.png")
plot_smooth(m07_sf_balance_aperiodic, "sf_balance_z", "m07_sf_balance_aperiodic_effect.png")
plot_smooth(m08_slow_frac_aperiodic, "slow_frac_z", "m08_slow_frac_aperiodic_effect.png")

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
model_comparison <- imap_dfr(model_list, extract_model_info) %>%
    arrange(AIC)

write_csv(model_comparison, file.path(out_dir, "model_comparison_v2.csv"))

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
      "sf_logratio",
      "sf_balance",
      "slow_frac",
      "fooof_offset",
      "fooof_exponent"
    ))
  )

for (nm in names(model_list)) {
    fitted_trial_df[[paste0(nm, "_fitted")]] <- fitted(model_list[[nm]])
    fitted_trial_df[[paste0(nm, "_resid")]] <- residuals(model_list[[nm]])
}

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values_v2.csv"))

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
    mean_pow_slow = mean(pow_slow, na.rm = TRUE),
    mean_pow_fast = mean(pow_fast, na.rm = TRUE),
    mean_sf_logratio = mean(sf_logratio, na.rm = TRUE),
    mean_sf_balance = mean(sf_balance, na.rm = TRUE),
    mean_slow_frac = mean(slow_frac, na.rm = TRUE),
    mean_fooof_offset = mean(fooof_offset, na.rm = TRUE),
    mean_fooof_exponent = mean(fooof_exponent, na.rm = TRUE),
    across(
      ends_with("_fitted"),
      ~ mean(.x, na.rm = TRUE),
      .names = "ga_{.col}"
    ),
    .groups = "drop"
  )

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values_v2.csv"))

# =========================================================
# DIAGNOSTICS
# =========================================================
iwalk(model_list, function(mod, nm) {
    png(file.path(out_dir, paste0(nm, "_diagnostics.png")), width = 1800, height = 1400, res = 180)
    gam.check(mod)
    dev.off()
})

# =========================================================
# OBSERVED VS FITTED
# =========================================================
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

# =========================================================
# INTERACTION SURFACES
# =========================================================
png(
    filename = file.path(out_dir, "m05_tensor_slow_fast_surface.png"),
    width = 1800,
    height = 1400,
    res = 200
)
vis.gam(
    m05_tensor_slow_fast,
    view = c("pow_slow_z", "pow_fast_z"),
    plot.type = "contour",
    too.far = 0.05,
    main = "Slow-Fast Alpha Tensor GAMM"
)
dev.off()

png(
  filename = file.path(out_dir, "m09_tensor_slow_fast_aperiodic_surface.png"),
  width = 1800,
  height = 1400,
  res = 200
)
vis.gam(
  m09_tensor_slow_fast_aperiodic,
  view = c("pow_slow_z", "pow_fast_z"),
  plot.type = "contour",
  too.far = 0.05,
  main = "Slow-Fast Alpha Tensor GAMM (Aperiodic Controlled)"
)
dev.off()

# =========================================================
# SAVE MODEL OBJECTS
# =========================================================
iwalk(model_list, function(mod, nm) {
  saveRDS(mod, file.path(out_dir, paste0(nm, ".rds")))
})

message("v2 GAMM workflow complete.")
message("Outputs saved to: ", out_dir)