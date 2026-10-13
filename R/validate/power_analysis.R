# =============================================================================
# power_analysis.R
# -----------------------------------------------------------------------------
# Validates whether the 75/25 train/holdout split provides adequate power
# across ALL significant source-space metric x ROI combinations in
# gamm_fitted — not just the best-performing one. Sizing the split to the
# hardest case (smallest effect size among significant hits) ensures you
# have sufficient power for every metric you intend to report, not just
# the most favourable one.
#
# For each significant single-metric source model in gamm_fitted:
#
#   PART 1 — Training-side power (simr::powerCurve)
#     Refits the bam() model as lmer() to extract real variance components,
#     then sweeps N from current training N to projected full-N target.
#     Run separately per cap_size group AND pooled.
#
#   PART 2 — Holdout stability (bootstrap)
#     Generates out-of-sample predictions on held-out subjects, bootstraps
#     RMSE and Pearson r to quantify metric stability at current holdout N.
#
#   PART 3 — Binding constraint summary
#     Reports which metric x ROI is the hardest case (highest N required
#     to reach 0.80/0.90 power), so the split fraction decision is sized
#     to the most conservative requirement, not the most optimistic one.
#
# Output:
#   results.duckdb: power_analysis table (one row per metric x ROI x
#                   cap_size_group x power_target)
#   <fig_dir>/: one power curve figure per metric x ROI
#   Console: full summary + binding constraint callout
#
# Run from R/gamm/ working directory:
#   setwd("R/gamm"); source("../validate/power_analysis.R")
# (sources gamm_defaults.R and helpers/*.R from there). Requires
# merge_core.R and source_core.R to have already run successfully.
#
# V2 (Oct 2026):
#   - Only gamm_fitted rows of the current feature_set are used.
#   - close_results_db() was called but never defined here (it lives in
#     gamm_combo_core.R, which was not sourced); now sourced.
#   - Holdout predictors are z-scored with the TRAINING mean/SD. They were
#     z-scored with the holdout's own mean/SD, which is not the scale the
#     model was fit on.
#   - powerCurve() had no `test =`, so simr tested the FIRST fixed effect,
#     laser_power_z: the curves were power to detect the laser effect, not
#     the metric. Now test = fixed(<metric part>, "z").
#   - The curve ran from the current N up to n_target_total without
#     simr::extend(), which powerCurve() needs to go beyond the subjects that
#     exist. Now extended first.
#   - Between / within split: the lmer refit includes both parts (linear),
#     power is for the primary part (default: subject mean), and training
#     data use the screen model's own scaling (mod$zscale).
#   - as.data.frame() on a powerCurve object fails (no method), so every
#     curve was lost to the tryCatch; values now come from summary().
#   - Per-cap-size subsets kept unused subject levels, which makes
#     extend() / powerCurve() fail ("number of levels differs"); droplevels().
#   - Holdout data prepared and predicted with the shared helpers
#     (prepare_holdout / predict_population). Note holdout_replicate.R is
#     now the replication test; the bootstrap here only sizes the split.
#   - Holdout prediction failed for every model: predict() on a
#     discrete = TRUE bam() with subject levels unseen in training stops
#     with "data is too long", and the tryCatch turned that into a silent
#     skip, so the holdout bootstrap never ran. Unseen subject (and
#     experiment) levels are now mapped to a training level and the random
#     effect is excluded, which gives the population-level prediction
#     (the placeholder level has no effect once the term is excluded;
#     tested).
#
# Methodological note on script frequency:
#   Run ONCE after first real results from source_core.R on interim data,
#   make a split-fraction decision, lock the split. Run a second time near
#   full N as a confirmatory check only. Do NOT run repeatedly and adjust
#   the split based on results — that is optional stopping and inflates
#   false positive rates.
#
# Reviewer-facing language:
#   "We conducted a simulation-based a priori power analysis (simr, R)
#   using an internal pilot design (N=39). Power curves were generated for
#   every significant metric x ROI combination identified in the screening
#   GAMM (p < 0.05, source-space, single-metric models), separately for
#   each EEG cap size group (32-, 62-, 64-channel) and pooled. The
#   train/holdout split fraction was selected to achieve >= 80% power for
#   the most conservative (smallest effect size) significant metric, with
#   holdout stability verified via subject-level bootstrap (500 resamples,
#   95% CI width for Pearson r < 0.20)."
# =============================================================================

library(DBI)
library(duckdb)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(lme4)
library(simr)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/hit_helpers.R")

# =============================================================================
# CONFIG
# =============================================================================
pa_cfg <- list(
  p_sig          = 0.05,    # significance threshold for model inclusion
  n_metrics_max  = 1L,      # single-metric models only (clean effect sizes)
  n_target_total = 678L,    # projected full N
  n_steps        = 12L,     # points on each power curve
  n_sim          = 200L,    # simr simulations per point (increase to 500+
                             # for final run; 200 is fast enough for interim)
  n_boot         = 500L,    # bootstrap resamples for holdout stability
  power_targets  = c(0.80, 0.90),
  ci_width_flag  = 0.20,    # flag holdout as unstable if CI width exceeds this
  fig_dir        = file.path(gamm_cfg$out_dir, "power_analysis")
)
dir.create(pa_cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# CONNECT
# =============================================================================
data_con    <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
results_con <- init_results_db(gamm_cfg$results_db_path)

# =============================================================================
# STEP 1: PULL ALL SIGNIFICANT SOURCE MODELS FROM gamm_fitted
# =============================================================================
message("== Step 1: Identifying significant models (feature_set ",
        gamm_cfg$feature_set, ") ==")

all_models <- dbGetQuery(results_con, "
  SELECT level, group_value, model_name, metrics, n_metrics,
         p_value, r_sq, dev_expl, aic, rds_path
  FROM gamm_fitted
  WHERE level       = 'source'
    AND n_metrics   = 1
    AND feature_set = ?
", params = list(gamm_cfg$feature_set)) %>%
  mutate(
    p_value_num = suppressWarnings(as.numeric(p_value)),
    r_sq        = suppressWarnings(as.numeric(r_sq))
  ) %>%
  filter(
    !is.na(p_value_num), p_value_num < pa_cfg$p_sig,
    !is.na(r_sq),
    !is.na(rds_path), file.exists(rds_path)
  ) %>%
  arrange(r_sq)   # ascending so binding constraint (smallest effect) is first

if (nrow(all_models) == 0) {
  stop(
    "No significant single-metric source models found in gamm_fitted ",
    "(p < ", pa_cfg$p_sig, ", level = 'source', feature_set = '",
    gamm_cfg$feature_set, "'). Run source_core.R first."
  )
}

message("  Found ", nrow(all_models), " significant model(s) to analyse:")
for (i in seq_len(nrow(all_models))) {
  message(sprintf("    [%d] %-28s @ %-14s  r_sq=%.3f  p=%.4f",
                  i, all_models$metrics[i], all_models$group_value[i],
                  all_models$r_sq[i], all_models$p_value_num[i]))
}

# =============================================================================
# STEP 2: SPLIT SUMMARY (printed once, applies to all models)
# =============================================================================
message("\n== Step 2: Split summary ==")

split_summary <- dbGetQuery(data_con, "
  SELECT cap_size, split_group, COUNT(*) AS n
  FROM subject_split
  GROUP BY cap_size, split_group
  ORDER BY cap_size, split_group
")
if (nrow(split_summary) == 0) stop("subject_split is empty — run merge_core.R first.")

message("  Current split:")
for (i in seq_len(nrow(split_summary))) {
  message(sprintf("    cap_size %-4s | %-8s | n = %d",
                  split_summary$cap_size[i],
                  split_summary$split_group[i],
                  split_summary$n[i]))
}
n_train   <- sum(split_summary$n[split_summary$split_group == "train"])
n_holdout <- sum(split_summary$n[split_summary$split_group == "holdout"])
message(sprintf("  TOTAL  train = %d  holdout = %d  (%.1f%% / %.1f%%)",
                n_train, n_holdout,
                100 * n_train   / (n_train + n_holdout),
                100 * n_holdout / (n_train + n_holdout)))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

base_filter <- function(df) {
  df %>% filter(!is.na(pain_rating), !is.na(laser_power),
                !is.na(trial_index), !is.na(age), !is.na(global_subjid))
}

# Mean / SD of every z-scored column, taken from the training data, so the
# holdout can be put on the scale the model was fit on.
scaling_params <- function(df, raw_cols) {
  setNames(lapply(raw_cols, function(cn) {
    c(mean = mean(df[[cn]], na.rm = TRUE), sd = stats::sd(df[[cn]], na.rm = TRUE))
  }), raw_cols)
}
apply_scaling <- function(x, p) (x - p[["mean"]]) / p[["sd"]]

# --- Load and z-score training data for one metric x ROI ---
# Uses the scaling stored on the screen model (mod$zscale) when present, so
# the lmer refit sees the same predictor scale as the bam(); adds the
# between / within split if the screen used it.
load_training_data <- function(data_con, roi, this_metric, mod_bam) {
  raw_metric <- str_remove(this_metric, "_z$")
  df <- load_source_metrics_wide(data_con, roi, metrics = raw_metric) %>%
    base_filter() %>%
    mutate(
      global_subjid = as.factor(as.character(global_subjid)),
      experiment_id = as.factor(as.character(experiment_id))
    )
  if (!raw_metric %in% names(df)) return(NULL)
  df <- df %>% filter(!is.na(.data[[raw_metric]]))
  need <- c("age", "laser_power", "trial_index", raw_metric)
  sc <- if (!is.null(mod_bam$zscale) && all(need %in% names(mod_bam$zscale)))
    mod_bam$zscale[need] else scaling_params(df, need)
  df <- df %>%
    mutate(
      age_z         = apply_scaling(age, sc$age),
      laser_power_z = apply_scaling(laser_power, sc$laser_power),
      trial_index_z = apply_scaling(trial_index, sc$trial_index)
    )
  df[[this_metric]] <- apply_scaling(df[[raw_metric]], sc[[raw_metric]])
  if (identical(mod_bam$gamm_meta$decomposition, "between_within"))
    df <- add_between_within(df, this_metric)
  attr(df, "scaling") <- sc
  df
}

# --- Refit as lmer for simr ---
# metric_vars: the metric column(s) as in the bam (one pooled column, or the
# between + within pair). Entered linearly: simr needs a parametric effect.
fit_lmer_for_simr <- function(df, metric_vars) {
  re_terms <- "(1 | global_subjid)"
  if (n_distinct(df$experiment_id) > 1)
    re_terms <- paste(re_terms, "+ (1 | experiment_id)")
  f <- as.formula(paste(
    "pain_rating ~ laser_power_z + trial_index_z + age_z +",
    paste(metric_vars, collapse = " + "), "+", re_terms
  ))
  tryCatch(
    lmer(f, data = df, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { message("  [WARN] lmer failed: ", conditionMessage(e)); NULL }
  )
}

# --- Run power curve for one subset (pooled or per cap_size) ---
# test_var: the effect whose power is estimated. Without `test =`, simr
# tests the FIRST fixed effect in the formula, which here is laser_power_z.
run_power_curve <- function(mod, along_var, n_current, n_target, label, n_steps, n_sim,
                            test_var) {
  breaks <- unique(round(seq(n_current, n_target, length.out = n_steps)))
  breaks <- breaks[breaks >= n_current]
  tryCatch(
    {
      # powerCurve() can only subset subjects that exist; extend() adds
      # simulated subjects so the curve can reach the projected N.
      mod_ext <- if (n_target > n_current) simr::extend(mod, along = along_var, n = n_target) else mod
      pc <- powerCurve(mod_ext, test = simr::fixed(test_var, "z"), along = along_var,
                       breaks = breaks, nsim = n_sim, progress = FALSE)
      # summary() is the documented accessor (as.data.frame() has no method
      # for powerCurve objects): one row per break.
      sm <- summary(pc)
      df <- data.frame(n = sm$nlevels, power = sm$mean, lower = sm$lower, upper = sm$upper)
      df$group <- label
      df
    },
    error = function(e) {
      message("  [WARN] powerCurve failed (", label, "): ", conditionMessage(e))
      NULL
    }
  )
}

# --- Extract N required to hit each power target ---
extract_n_required <- function(curve_df, label, power_targets) {
  if (is.null(curve_df)) return(NULL)
  map_dfr(power_targets, function(thresh) {
    hit <- curve_df %>% filter(power >= thresh) %>% slice(1)
    tibble(
      group        = label,
      power_target = thresh,
      n_required   = if (nrow(hit) > 0) as.integer(hit$n) else NA_integer_,
      achieved     = nrow(hit) > 0,
      max_power    = max(curve_df$power)
    )
  })
}

# --- Bootstrap holdout stability ---
# Held-out rows prepared and predicted exactly as in holdout_replicate.R
# (prepare_holdout / predict_population in helpers/hit_helpers.R).
bootstrap_holdout <- function(data_con, roi, this_metric, mod_bam, n_boot, ci_flag) {
  df_ho <- prepare_holdout(data_con, "source", roi, this_metric, mod_bam, gamm_cfg)
  if (is.null(df_ho) || nrow(df_ho) == 0) return(NULL)

  df_ho$predicted <- tryCatch(
    predict_population(mod_bam, df_ho),
    error = function(e) {
      message("  [WARN] holdout prediction failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(df_ho$predicted)) return(NULL)
  df_ho <- df_ho %>% filter(!is.na(predicted))

  rmse_fn <- function(a, p) sqrt(mean((a - p)^2, na.rm = TRUE))
  r_fn    <- function(a, p) cor(a, p, use = "complete.obs")

  by_subj <- split(df_ho, as.character(df_ho$global_subjid))
  boot_mat <- replicate(n_boot, {
    subs    <- sample(names(by_subj), replace = TRUE)
    boot_df <- bind_rows(by_subj[subs])
    c(rmse = rmse_fn(boot_df$pain_rating, boot_df$predicted),
      r    = r_fn(boot_df$pain_rating,    boot_df$predicted))
  })
  boot_df <- as_tibble(t(boot_mat))

  ci_width_r <- diff(quantile(boot_df$r, c(0.025, 0.975)))
  n_ho_subj  <- n_distinct(df_ho$global_subjid)

  list(
    n_holdout   = n_ho_subj,
    rmse_mean   = mean(boot_df$rmse),
    rmse_sd     = sd(boot_df$rmse),
    rmse_ci_lo  = quantile(boot_df$rmse, 0.025),
    rmse_ci_hi  = quantile(boot_df$rmse, 0.975),
    r_mean      = mean(boot_df$r),
    r_sd        = sd(boot_df$r),
    r_ci_lo     = quantile(boot_df$r, 0.025),
    r_ci_hi     = quantile(boot_df$r, 0.975),
    ci_width_r  = ci_width_r,
    stable      = ci_width_r <= ci_flag,
    cap_counts  = df_ho %>% distinct(global_subjid, cap_size) %>% count(cap_size)
  )
}

# =============================================================================
# STEP 3: MAIN LOOP — one full analysis per significant model
# =============================================================================
message("\n== Step 3: Running power analysis for each significant model ==")

all_power_rows  <- list()

for (i in seq_len(nrow(all_models))) {
  row         <- all_models[i, ]
  this_metric <- row$metrics
  this_roi    <- row$group_value
  mod_bam     <- readRDS(row$rds_path)

  message(sprintf("\n[%d/%d] %s @ %s  (r_sq=%.3f, p=%.4f)",
                  i, nrow(all_models), this_metric, this_roi,
                  row$r_sq, row$p_value_num))

  # Load training data
  df_train <- load_training_data(data_con, this_roi, this_metric, mod_bam)
  metric_vars <- model_metric_vars(mod_bam, this_metric)
  test_var    <- primary_var(mod_bam, this_metric)
  if (is.null(df_train) || nrow(df_train) == 0) {
    message("  Skipping — training data unavailable.")
    next
  }

  n_train_subj <- n_distinct(df_train$global_subjid)
  cap_counts   <- df_train %>% distinct(global_subjid, cap_size) %>% count(cap_size)
  message("  Training: ", n_train_subj, " subjects  ", nrow(df_train), " trials")
  for (j in seq_len(nrow(cap_counts))) {
    message("    cap_size ", cap_counts$cap_size[j], ": ", cap_counts$n[j])
  }

  # lmer refit
  mod_lmer <- fit_lmer_for_simr(df_train, metric_vars)
  if (is.null(mod_lmer)) next

  message("  lmer OK — subject RE SD: ",
          round(as.data.frame(VarCorr(mod_lmer))$sdcor[1], 4),
          "  residual SD: ",
          round(attr(VarCorr(mod_lmer), "sc"), 4),
          "  beta: ",
          round(fixef(mod_lmer)[test_var], 4), "  (power is for ", test_var, ")")

  # --- Pooled power curve ---
  curve_pooled <- run_power_curve(
    mod_lmer, "global_subjid",
    n_current = n_train_subj,
    n_target  = pa_cfg$n_target_total,
    label     = "pooled",
    n_steps   = pa_cfg$n_steps,
    n_sim     = pa_cfg$n_sim,
    test_var  = test_var
  )
  power_rows <- extract_n_required(curve_pooled, "pooled", pa_cfg$power_targets)

  # --- Per cap_size power curves ---
  all_curves_this_model <- list(curve_pooled)

  for (this_cap in sort(unique(as.character(df_train$cap_size)))) {
    df_cap   <- droplevels(filter(df_train, cap_size == this_cap))  # unused subject levels break extend()
    n_cap    <- n_distinct(df_cap$global_subjid)
    if (n_cap < 5) {
      message("  Skipping cap_size ", this_cap, " — only ", n_cap, " subjects.")
      next
    }

    mod_cap <- fit_lmer_for_simr(df_cap, metric_vars)
    if (is.null(mod_cap)) next

    n_target_cap <- round(pa_cfg$n_target_total *
                            (cap_counts$n[cap_counts$cap_size == this_cap] / n_train_subj))
    lbl <- paste0("cap_size_", this_cap)

    curve_cap <- run_power_curve(
      mod_cap, "global_subjid",
      n_current = n_cap,
      n_target  = max(n_target_cap, n_cap + 5L),
      label     = lbl,
      n_steps   = pa_cfg$n_steps,
      n_sim     = pa_cfg$n_sim,
      test_var  = test_var
    )
    cap_rows <- extract_n_required(curve_cap, lbl, pa_cfg$power_targets)
    if (!is.null(cap_rows)) power_rows <- bind_rows(power_rows, cap_rows)
    all_curves_this_model[[lbl]] <- curve_cap
  }

  # Tag power rows with metric/ROI identity
  if (!is.null(power_rows)) {
    power_rows <- power_rows %>%
      mutate(metric = this_metric, roi = this_roi, fitted_at = as.character(Sys.time()))
    all_power_rows[[length(all_power_rows) + 1]] <- power_rows

    message("  Power summary:")
    for (j in seq_len(nrow(power_rows))) {
      achieved_str <- if (power_rows$achieved[j]) "ACHIEVED" else "NOT ACHIEVED in range"
      message(sprintf("    %-20s  power=%.0f%%  n_required=%-5s  %s",
                      power_rows$group[j],
                      power_rows$power_target[j] * 100,
                      ifelse(is.na(power_rows$n_required[j]),
                             ">max", power_rows$n_required[j]),
                      achieved_str))
    }
  }

  # --- Holdout stability ---
  message("  Running holdout bootstrap (", pa_cfg$n_boot, " resamples)...")
  ho <- bootstrap_holdout(data_con, this_roi, this_metric, mod_bam,
                          pa_cfg$n_boot, pa_cfg$ci_width_flag)
  if (!is.null(ho)) {
    message(sprintf("  Holdout: %d subjects  r=%.3f [%.3f, %.3f]  RMSE=%.3f [%.3f, %.3f]  stable=%s",
                    ho$n_holdout,
                    ho$r_mean, ho$r_ci_lo, ho$r_ci_hi,
                    ho$rmse_mean, ho$rmse_ci_lo, ho$rmse_ci_hi,
                    ho$stable))
    if (!ho$stable) {
      message("  [FLAG] Holdout CI width (", round(ho$ci_width_r, 3),
              ") exceeds threshold (", pa_cfg$ci_width_flag,
              ") — consider increasing split_train_frac in merge_defaults.R.")
    }

    ho_row <- tibble(
      metric       = this_metric,
      roi          = this_roi,
      group        = "holdout_bootstrap",
      power_target = NA_real_,
      n_required   = ho$n_holdout,
      achieved     = ho$stable,
      max_power    = ho$r_mean,
      fitted_at    = as.character(Sys.time())
    )
    all_power_rows[[length(all_power_rows) + 1]] <- ho_row
  }

  # --- Power curve figure ---
  plot_data <- compact(all_curves_this_model) %>% bind_rows()
  if (nrow(plot_data) > 0) {
    p <- ggplot(plot_data, aes(x = n, y = power, colour = group, fill = group)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      geom_hline(yintercept = pa_cfg$power_targets, linetype = "dashed",
                 colour = "grey40", linewidth = 0.5) +
      annotate("text",
               x     = min(plot_data$n),
               y     = pa_cfg$power_targets,
               label = paste0(pa_cfg$power_targets * 100, "%"),
               hjust = -0.1, vjust = -0.4, size = 3, colour = "grey40") +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      labs(
        title    = paste0("Power curve: ", this_metric, " @ ", this_roi),
        subtitle = paste0("simr (nsim=", pa_cfg$n_sim, " per point)  |  ",
                          "lmer approximation of bam()  |  ",
                          "r²=", round(row$r_sq, 3),
                          "  p=", signif(row$p_value_num, 3)),
        x        = "N subjects (training set)",
        y        = "Estimated power",
        colour   = NULL, fill = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")

    fig_name <- paste0("power_curve_",
                       str_replace_all(this_metric, "_z$", ""), "_",
                       str_replace_all(this_roi, "-", "_"), ".png")
    fig_path <- file.path(pa_cfg$fig_dir, fig_name)
    ggsave(fig_path, p, width = 8, height = 5, dpi = 150)
    message("  Figure saved: ", fig_path)
  }
}

# =============================================================================
# STEP 4: BINDING CONSTRAINT SUMMARY
# =============================================================================
message("\n", strrep("=", 60))
message("== BINDING CONSTRAINT SUMMARY ==")
message(strrep("=", 60))

all_power_df <- bind_rows(all_power_rows) %>%
  filter(!is.na(power_target))   # exclude holdout_bootstrap rows for this summary

if (nrow(all_power_df) > 0) {

  # Per power target: which metric x ROI x group requires the MOST subjects?
  for (thresh in pa_cfg$power_targets) {
    message(sprintf("\nFor %.0f%% power:", thresh * 100))
    sub <- all_power_df %>%
      filter(power_target == thresh) %>%
      arrange(desc(n_required))
    if (nrow(sub) == 0) { message("  No data."); next }

    binding <- sub %>% filter(!is.na(n_required)) %>% slice(1)
    if (nrow(binding) > 0) {
      message(sprintf("  Binding constraint: %s @ %s (%s)  n_required = %d",
                      binding$metric, binding$roi, binding$group, binding$n_required))
    }

    not_achieved <- sub %>% filter(!achieved)
    if (nrow(not_achieved) > 0) {
      message("  NOT achieved within projected N range:")
      for (j in seq_len(nrow(not_achieved))) {
        message(sprintf("    %s @ %s (%s)  max_power = %.2f",
                        not_achieved$metric[j], not_achieved$roi[j],
                        not_achieved$group[j], not_achieved$max_power[j]))
      }
    } else {
      message("  All significant models achieve this threshold within projected N.")
    }
  }

  # Verdict
  message("\nVerdict:")
  binding_80 <- all_power_df %>%
    filter(power_target == 0.80, !is.na(n_required)) %>%
    arrange(desc(n_required)) %>% slice(1)

  if (nrow(binding_80) > 0 && !is.na(binding_80$n_required)) {
    current_train <- sum(split_summary$n[split_summary$split_group == "train"])
    if (binding_80$n_required <= current_train) {
      message(sprintf("  Current training N (%d) is SUFFICIENT for 80%% power",
                      current_train))
      message(sprintf("  across all significant models (binding constraint: n=%d).",
                      binding_80$n_required))
    } else {
      message(sprintf("  [ACTION REQUIRED] Current training N (%d) is INSUFFICIENT.",
                      current_train))
      message(sprintf("  Need n=%d for 80%% power on the binding metric (%s @ %s).",
                      binding_80$n_required, binding_80$metric, binding_80$roi))
      new_frac <- round(binding_80$n_required / pa_cfg$n_target_total + 0.01, 2)
      message(sprintf("  Consider increasing split_train_frac in merge_defaults.R to at least %.2f.",
                      new_frac))
    }
  }
}

# =============================================================================
# STEP 5: WRITE ALL RESULTS TO results.duckdb
# =============================================================================
message("\n== Step 5: Writing to results.duckdb ==")

invisible(dbExecute(results_con, "
  CREATE TABLE IF NOT EXISTS power_analysis (
    fitted_at    VARCHAR,
    metric       VARCHAR,
    roi          VARCHAR,
    group_label  VARCHAR,
    power_target DOUBLE,
    n_required   INTEGER,
    achieved     BOOLEAN,
    max_power    DOUBLE
  )
"))

final_rows <- bind_rows(all_power_rows) %>%
  rename(group_label = group)

# Upsert all rows for this run (delete this metric x ROI, then insert).
# No primary key: holdout_bootstrap rows have power_target = NULL, which a
# PRIMARY KEY column cannot hold.
if (nrow(final_rows) > 0) {
  keys <- distinct(final_rows, metric, roi)
  for (k in seq_len(nrow(keys))) {
    invisible(dbExecute(results_con,
      "DELETE FROM power_analysis WHERE metric = ? AND roi = ?",
      params = list(keys$metric[k], keys$roi[k])
    ))
  }
  dbAppendTable(results_con, "power_analysis",
                as.data.frame(final_rows %>% select(fitted_at, metric, roi, group_label,
                                                    power_target, n_required, achieved, max_power)))
  message("  Written ", nrow(final_rows), " rows to power_analysis table.")
}

disconnect_db(data_con)
close_results_db(results_con)

message("\nPower analysis complete.")
message("Results: ", gamm_cfg$results_db_path, "  (table: power_analysis)")
message("Figures: ", pa_cfg$fig_dir)