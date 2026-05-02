# =============================================================================
# run_gamm_source.R
# -----------------------------------------------------------------------------
# Fits the source-space GAMM model family for each ROI.
#
# Model family (per ROI)
# ───────────────────────
#   m00  Baseline (covariates only)
#   m01  + s(BI_pre)          sub-band balance index
#   m02  + s(LR_pre)          log ratio
#   m03  + s(CoG_pre)         centre of gravity (peak alpha frequency proxy)
#   m04  + s(psi_cog)         BI_pre × (CoG_pre − 10) interaction
#   m05  + te(pow_slow, pow_fast)   raw tensor (slow × fast)
#   m06  + sin_phase + cos_phase    Hilbert phase at stimulus onset
#   m07  + s(delta_ERD)       sub-band ERD asymmetry
#   m08  + s(n2p2_amp)        N2-P2 LEP amplitude
#   m09  + s(BI_pre) + s(delta_ERD)           pre + post combined
#   m10  + s(BI_pre) + s(n2p2_amp)            alpha + LEP combined
#   m11  + s(BI_pre) + s(delta_ERD) + s(n2p2_amp)  full combined
#   m12  FOOOF-controlled: m01 baseline with aperiodic terms
#   m13  FOOOF-controlled: + s(BI_pre)
#   m14  FOOOF-controlled: + s(LR_pre)
#
# Output (written to <out_dir>/<roi>/)
# ─────────────────────────────────────
#   model_comparison_<roi>.csv       AIC / BIC / deviance explained per model
#   <model_name>_summary.txt         mgcv summary
#   <model_name>_diagnostics.png     gam.check plots
#   <model_name>_smooth_<term>.png   smooth effect plot per significant term
#   trial_level_fitted_<roi>.csv     per-trial fitted + residual values
#   subject_ga_fitted_<roi>.csv      subject-level GA summary
#
# One cross-ROI summary is written to the out_dir root:
#   model_comparison_all_rois.csv    all ROIs × all models in one table
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(mgcv)
library(tibble)
library(ggplot2)

# =============================================================================
# USER SETTINGS
# =============================================================================
data_file      <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/source_pain_master.csv"
out_dir        <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_source"
min_fooof_r2   <- 0.80
nthreads_use   <- 4L
use_discrete   <- TRUE
min_trials_roi <- 20L   # skip ROI × subject if fewer trials survive QC

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================
clean_names_local <- function(x) {
  x %>% str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2",  "r2") %>%
    str_replace_all("\\.",   "_")
}

zscore <- function(x) as.numeric(scale(x))

safe_factor <- function(x) as.factor(as.character(x))

extract_model_info <- function(mod, nm, roi) {
  tibble(
    roi        = roi,
    model_name = nm,
    AIC        = AIC(mod),
    BIC        = BIC(mod),
    logLik     = as.numeric(logLik(mod)),
    dev_expl   = summary(mod)$dev.expl,
    r_sq       = summary(mod)$r.sq,
    n          = nobs(mod)
  )
}

# Save a smooth-effect plot for a named term if that term is in the model.
plot_smooth_if_present <- function(model, term, file_path) {
  smooth_terms <- sapply(model$smooth, function(s) {
    if (!is.null(s$label)) s$label else paste(s$term, collapse = ":")
  })
  idx <- which(str_detect(smooth_terms, fixed(term)))[1]
  if (is.na(idx)) return(invisible(NULL))

  png(file_path, width = 1400, height = 1100, res = 180)
  plot(model, select = idx, shade = TRUE, shade.col = "lightblue",
       main = paste("Effect of", term), xlab = term,
       ylab = "Partial effect on pain rating")
  abline(h = 0, lty = 2)
  dev.off()
}

# Fit one bam() with graceful failure — returns NULL on error.
safe_bam <- function(formula, data, nm, roi) {
  tryCatch(
    bam(formula  = formula,
        data     = data,
        method   = "fREML",
        discrete = use_discrete,
        nthreads = nthreads_use),
    error = function(e) {
      message("  [WARN] Model ", nm, " (ROI ", roi, ") failed: ", conditionMessage(e))
      NULL
    }
  )
}

# =============================================================================
# LOAD DATA
# =============================================================================
if (!file.exists(data_file)) {
  stop("source_pain_master.csv not found: ", data_file,
       "\nRun merge_source_spectral.R first.")
}

df_raw <- read_csv(data_file, show_col_types = FALSE)
names(df_raw) <- clean_names_local(names(df_raw))

required_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
  "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating",
  "roi", "BI_pre", "LR_pre", "CoG_pre", "pow_slow", "pow_fast"
)
missing <- setdiff(required_cols, names(df_raw))
if (length(missing) > 0) {
  stop("source_pain_master missing required columns: ", paste(missing, collapse = ", "))
}

df_raw <- df_raw %>%
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
    pain_rating     = suppressWarnings(as.numeric(pain_rating)),
    roi             = as.character(roi)
  )

all_rois <- sort(unique(df_raw$roi))
message("ROIs found: ", paste(all_rois, collapse = ", "))

# =============================================================================
# BASE FORMULAS
# These are constructed once and shared across ROIs.
# All continuous predictors use z-scored versions; random effects
# account for subject and experiment clustering.
# =============================================================================
base_formula <- pain_rating ~
  s(laser_power_z,  k = 10) +
  s(trial_index_z,  k = 10) +
  s(age_z,          k = 10) +
  sex + cap_size +
  s(global_subjid,  bs = "re") +
  s(experiment_id,  bs = "re")

base_formula_fooof <- update(base_formula, . ~ . +
  s(fooof_offset_z,   k = 10) +
  s(fooof_exponent_z, k = 10))

# =============================================================================
# HELPER: build a named list of model definitions for one ROI's data.
# Each element is a one-sided update formula for the predictor of interest.
# Models that require columns absent from the data are silently skipped
# by safe_bam returning NULL.
# =============================================================================
model_definitions <- list(
  m00_baseline          = NULL,     # base formula only
  m01_bi_pre            = ~ . + s(BI_pre_z,    k = 10),
  m02_lr_pre            = ~ . + s(LR_pre_z,    k = 10),
  m03_cog_pre           = ~ . + s(CoG_pre_z,   k = 10),
  m04_psi_cog           = ~ . + s(psi_cog_z,   k = 10),
  m05_tensor_slow_fast  = ~ . + te(pow_slow_z, pow_fast_z, k = c(8, 8)),
  m06_phase             = ~ . + sin_phase + cos_phase,   # linear circular regressors
  m07_delta_erd         = ~ . + s(delta_ERD_z, k = 10),
  m08_n2p2_amp          = ~ . + s(n2p2_amp_z,  k = 10),
  m09_bi_delta_erd      = ~ . + s(BI_pre_z,    k = 10) + s(delta_ERD_z, k = 10),
  m10_bi_n2p2           = ~ . + s(BI_pre_z,    k = 10) + s(n2p2_amp_z,  k = 10),
  m11_full_combined     = ~ . + s(BI_pre_z,    k = 10) + s(delta_ERD_z, k = 10) +
                                s(n2p2_amp_z,  k = 10),
  m12_fooof_baseline    = NULL,     # base_formula_fooof only
  m13_fooof_bi_pre      = ~ . + s(BI_pre_z,    k = 10),
  m14_fooof_lr_pre      = ~ . + s(LR_pre_z,    k = 10)
)

# Smooth-effect terms to plot per model (NULL = no scalar smooth to highlight)
model_plot_terms <- list(
  m00_baseline          = NULL,
  m01_bi_pre            = "BI_pre_z",
  m02_lr_pre            = "LR_pre_z",
  m03_cog_pre           = "CoG_pre_z",
  m04_psi_cog           = "psi_cog_z",
  m05_tensor_slow_fast  = NULL,     # 2-D surface plotted separately
  m06_phase             = NULL,
  m07_delta_erd         = "delta_ERD_z",
  m08_n2p2_amp          = "n2p2_amp_z",
  m09_bi_delta_erd      = "BI_pre_z",
  m10_bi_n2p2           = "BI_pre_z",
  m11_full_combined     = "BI_pre_z",
  m12_fooof_baseline    = NULL,
  m13_fooof_bi_pre      = "BI_pre_z",
  m14_fooof_lr_pre      = "LR_pre_z"
)

# =============================================================================
# PER-ROI LOOP
# =============================================================================
all_roi_comparison <- list()

for (this_roi in all_rois) {
  message("\n", strrep("=", 60))
  message("  ROI: ", this_roi)
  message(strrep("=", 60))

  roi_dir <- file.path(out_dir, this_roi)
  dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Subset + QC filter ────────────────────────────────────────────────────
  df_roi <- df_raw %>%
    filter(roi == this_roi) %>%
    filter(
      !is.na(pain_rating), !is.na(laser_power),
      !is.na(trial_index), !is.na(global_subjid),
      !is.na(experiment_id), !is.na(age), !is.na(sex), !is.na(cap_size),
      !is.na(BI_pre), !is.na(LR_pre), !is.na(CoG_pre),
      !is.na(pow_slow), !is.na(pow_fast)
    )

  if ("fooof_r2" %in% names(df_roi)) {
    df_roi <- filter(df_roi, is.na(fooof_r2) | fooof_r2 >= min_fooof_r2)
  }

  if (nrow(df_roi) < min_trials_roi) {
    message("  Skipping: only ", nrow(df_roi), " trials after QC (min = ", min_trials_roi, ").")
    next
  }
  message("  Trials after QC: ", nrow(df_roi),
          "  Subjects: ", n_distinct(df_roi$subjid_uid))

  # ── Z-score continuous predictors ─────────────────────────────────────────
  # Predictors are z-scored within each ROI's available data to ensure
  # smooth k-values are well-scaled and GAMM comparison is on equal footing.
  df_roi <- df_roi %>%
    mutate(
      age_z         = zscore(age),
      laser_power_z = zscore(laser_power),
      trial_index_z = zscore(trial_index),
      pow_slow_z    = zscore(pow_slow),
      pow_fast_z    = zscore(pow_fast),
      BI_pre_z      = zscore(BI_pre),
      LR_pre_z      = zscore(LR_pre),
      CoG_pre_z     = zscore(CoG_pre)
    )

  # Optional predictors — z-score only if present and non-constant
  opt_z <- list(
    psi_cog_z      = "psi_cog",
    delta_ERD_z    = "delta_ERD",
    n2p2_amp_z     = "n2p2_amp",
    fooof_offset_z = "fooof_offset",
    fooof_exponent_z = "fooof_exponent"
  )
  for (zname in names(opt_z)) {
    raw <- opt_z[[zname]]
    if (raw %in% names(df_roi) && sum(!is.na(df_roi[[raw]])) > 1) {
      df_roi[[zname]] <- zscore(df_roi[[raw]])
    }
  }

  # sin_phase / cos_phase are already bounded [−1, 1]; no z-scoring needed
  for (pc in c("sin_phase", "cos_phase")) {
    if (!pc %in% names(df_roi)) df_roi[[pc]] <- NA_real_
  }

  write_csv(df_roi,
            file.path(roi_dir, paste0("model_input_", this_roi, ".csv")))

  # ── Fit models ────────────────────────────────────────────────────────────
  model_list <- list()

  for (nm in names(model_definitions)) {
    message("  Fitting ", nm, "...")

    # Decide which base formula to start from
    is_fooof_model <- str_starts(nm, "m12") | str_starts(nm, "m13") | str_starts(nm, "m14")
    base <- if (is_fooof_model) base_formula_fooof else base_formula

    pred_update <- model_definitions[[nm]]
    formula_use <- if (is.null(pred_update)) base else update(base, pred_update)

    mod <- safe_bam(formula_use, df_roi, nm, this_roi)
    if (!is.null(mod)) model_list[[nm]] <- mod
  }

  if (length(model_list) == 0L) {
    message("  No models fitted for ROI ", this_roi, " — skipping outputs.")
    next
  }

  # ── Save summaries + diagnostics ──────────────────────────────────────────
  iwalk(model_list, function(mod, nm) {
    capture.output(summary(mod),
                   file = file.path(roi_dir, paste0(nm, "_summary.txt")))

    png(file.path(roi_dir, paste0(nm, "_diagnostics.png")),
        width = 1800, height = 1400, res = 180)
    gam.check(mod)
    dev.off()
  })

  # ── Smooth effect plots ────────────────────────────────────────────────────
  iwalk(model_list, function(mod, nm) {
    term <- model_plot_terms[[nm]]
    if (!is.null(term)) {
      plot_smooth_if_present(
        mod, term,
        file.path(roi_dir, paste0(nm, "_smooth_", term, ".png"))
      )
    }
  })

  # ── Tensor surface plot (m05) ──────────────────────────────────────────────
  if ("m05_tensor_slow_fast" %in% names(model_list)) {
    png(file.path(roi_dir, "m05_tensor_slow_fast_surface.png"),
        width = 1800, height = 1400, res = 200)
    vis.gam(model_list[["m05_tensor_slow_fast"]],
            view = c("pow_slow_z", "pow_fast_z"),
            plot.type = "contour", too.far = 0.05,
            main = paste0("Slow × Fast tensor — ", this_roi))
    dev.off()
  }

  # ── Model comparison table ────────────────────────────────────────────────
  roi_comparison <- imap_dfr(model_list, ~extract_model_info(.x, .y, this_roi)) %>%
    mutate(delta_AIC = AIC - min(AIC)) %>%
    arrange(AIC)

  write_csv(roi_comparison,
            file.path(roi_dir, paste0("model_comparison_", this_roi, ".csv")))
  all_roi_comparison[[this_roi]] <- roi_comparison

  # ── Observed vs fitted ────────────────────────────────────────────────────
  fitted_df <- df_roi %>%
    select(
      experiment_name, experiment_id, subjid, subjid_uid, global_subjid,
      age, sex, cap_size, trial, trial_index, laser_power, pain_rating, roi,
      any_of(c("BI_pre", "LR_pre", "CoG_pre", "psi_cog",
               "pow_slow", "pow_fast", "delta_ERD", "n2p2_amp",
               "slow_phase", "sin_phase", "cos_phase",
               "fooof_offset", "fooof_exponent"))
    )

  for (nm in names(model_list)) {
    fitted_df[[paste0(nm, "_fitted")]] <- fitted(model_list[[nm]])
    fitted_df[[paste0(nm, "_resid")]]  <- residuals(model_list[[nm]])
  }

  write_csv(fitted_df,
            file.path(roi_dir, paste0("trial_level_fitted_", this_roi, ".csv")))

  # Subject GA
  subject_ga <- fitted_df %>%
    group_by(experiment_name, experiment_id, subjid, subjid_uid,
             global_subjid, age, sex, cap_size, roi) %>%
    summarise(
      n_trials         = n(),
      mean_pain_rating = mean(pain_rating, na.rm = TRUE),
      mean_laser_power = mean(laser_power, na.rm = TRUE),
      mean_BI_pre      = mean(BI_pre,      na.rm = TRUE),
      mean_LR_pre      = mean(LR_pre,      na.rm = TRUE),
      mean_CoG_pre     = mean(CoG_pre,     na.rm = TRUE),
      across(ends_with("_fitted"), ~mean(.x, na.rm = TRUE), .names = "ga_{.col}"),
      .groups = "drop"
    )
  write_csv(subject_ga,
            file.path(roi_dir, paste0("subject_ga_fitted_", this_roi, ".csv")))

  # Observed vs fitted scatter for the best-fitting model
  best_nm <- roi_comparison$model_name[1]
  if (paste0(best_nm, "_fitted") %in% names(fitted_df)) {
    p <- ggplot(fitted_df,
                aes(x = .data[[paste0(best_nm, "_fitted")]], y = pain_rating)) +
      geom_point(alpha = 0.15, size = 0.8) +
      geom_smooth(method = "lm", se = FALSE, colour = "steelblue") +
      labs(title = paste0(this_roi, " — best model: ", best_nm),
           x = "Fitted pain", y = "Observed pain") +
      theme_minimal()
    ggsave(file.path(roi_dir, paste0("obs_vs_fitted_best_", this_roi, ".png")),
           p, width = 6, height = 5, dpi = 200)
  }

  # Save model RDS objects
  iwalk(model_list, function(mod, nm) {
    saveRDS(mod, file.path(roi_dir, paste0(nm, ".rds")))
  })

  message("  Done: ", this_roi)
}

# =============================================================================
# CROSS-ROI SUMMARY
# =============================================================================
if (length(all_roi_comparison) > 0L) {
  comparison_all <- bind_rows(all_roi_comparison) %>%
    arrange(roi, AIC)

  write_csv(comparison_all,
            file.path(out_dir, "model_comparison_all_rois.csv"))

  # Best model per ROI
  best_per_roi <- comparison_all %>%
    group_by(roi) %>%
    slice_min(AIC, n = 1L) %>%
    ungroup()
  write_csv(best_per_roi, file.path(out_dir, "best_model_per_roi.csv"))

  message("\nBest model per ROI:")
  print(select(best_per_roi, roi, model_name, AIC, dev_expl))
}

message("\nSource GAMM workflow complete.")
message("Outputs saved to: ", out_dir)