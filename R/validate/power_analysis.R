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
# For each significant single-metric source model in gamm_fitted (current
# feature_set and decomposition):
#
#   PART 1 — Training-side power (simr::powerCurve)
#     Refits the bam() model as lmer() on the screen model's own rows
#     (mod$model: same subjects, same QC, same scaling) to extract real
#     variance components, then sweeps training N from the current N to the
#     projected full N (the whole sample in training = split fraction 1).
#     Run separately per cap_size group AND pooled, at the observed effect
#     and at shrunken effects (effect_scales; see below).
#
#   PART 2 — Holdout stability (bootstrap)
#     Generates out-of-sample predictions on held-out subjects, bootstraps
#     RMSE and Pearson r to quantify metric stability at current holdout N.
#
#   PART 3 — Binding constraint summary
#     Reports which metric x ROI is the hardest case (highest total training
#     N required to reach 0.80/0.90 power), at the most conservative effect
#     scale, so the split fraction decision is sized to the most
#     conservative requirement, not the most optimistic one.
#
# Effect scales (winner's curse): every model analysed here was SELECTED
# because its effect looked large in this same training sample, so its
# estimated effect is biased upwards and power at that effect is optimistic.
# Each curve is therefore also run with the metric's coefficient multiplied
# by the smaller values in pa_cfg$effect_scales (default 0.5), and the
# verdict uses the smallest scale.
#
# Output:
#   results.duckdb: power_analysis table (one row per metric x ROI x
#                   group x effect_scale x power_target, plus one
#                   holdout_bootstrap row per metric x ROI)
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
# V3 (Oct 2026):
#   - Training data = the screen model's own frame (mod$model): exactly the
#     rows, QC (p5_flag / FOOOF filters, missing sex / cap size dropped) and
#     scaling the hit was found on. The reloaded data used before skipped
#     the screen's QC.
#   - The lmer refit includes sex and cap size (when they vary), as the
#     screen model does.
#   - effect_scales: power at the observed effect AND at shrunken effects
#     (winner's curse); the verdict uses the smallest scale.
#   - Per-cap-size requirements are converted to TOTAL training N
#     (n_required_total = n_required / that cap size's share of training
#     subjects) before they are compared with the current training N; the
#     cap-size N was compared with the total training N before.
#   - Models are processed weakest evidence first (largest p). They were
#     ordered by whole-model R², which is dominated by laser power and the
#     subject intercepts, not by the metric.
#   - power_analysis gains feature_set, decomposition, effect_scale,
#     n_current, n_required_total and dedicated holdout_* columns (the
#     holdout r used to be stored in max_power). A table in the old layout is
#     renamed to power_analysis_pre_v3, not dropped.
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
#   each EEG cap size group (32-, 62-, 64-channel) and pooled, at the
#   observed effect and at half the observed effect to allow for
#   selection-inflated estimates. The train/holdout split fraction was
#   selected to achieve >= 80% power for the most conservative (smallest
#   effect size) significant metric at the reduced effect, with holdout
#   stability verified via subject-level bootstrap (500 resamples, 95% CI
#   width for Pearson r < 0.20)."
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
  p_sig          = gamm_cfg$hits_p_sig,  # significance threshold for model inclusion
  n_target_total = 678L,    # projected full N (= training N if the split fraction were 1)
  n_steps        = 12L,     # points on each power curve
  n_sim          = 200L,    # simr simulations per point (increase to 500+
                            # for final run; 200 is fast enough for interim)
  n_boot         = 500L,    # bootstrap resamples for holdout stability
  power_targets  = c(0.80, 0.90),
  # Multipliers on the metric's estimated coefficient (winner's curse). 1 =
  # observed effect; the smallest value drives the verdict. Each extra scale
  # adds a full set of power curves (runtime).
  effect_scales  = c(1.0, 0.5),
  ci_width_flag  = 0.20,    # flag holdout as unstable if CI width exceeds this
  fig_dir        = file.path(gamm_cfg$out_dir, "power_analysis")
)
dir.create(pa_cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
decomp <- gamm_cfg$metric_decomposition

# =============================================================================
# CONNECT
# =============================================================================
data_con    <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
results_con <- init_results_db(gamm_cfg$results_db_path)

# =============================================================================
# STEP 1: PULL ALL SIGNIFICANT SOURCE MODELS FROM gamm_fitted
# =============================================================================
message("== Step 1: Identifying significant models (feature_set ",
        gamm_cfg$feature_set, ", decomposition ", decomp, ") ==")

all_models <- load_hits(results_con, gamm_cfg$feature_set, p_sig = pa_cfg$p_sig,
                        levels = "source", decomposition = decomp) %>%
  filter(rds_ok) %>%
  arrange(desc(p_num))   # weakest evidence first

if (nrow(all_models) == 0) {
  stop(
    "No significant single-metric source models found in gamm_fitted ",
    "(p < ", pa_cfg$p_sig, ", level = 'source', feature_set = '",
    gamm_cfg$feature_set, "', decomposition = '", decomp, "'). Run source_core.R first."
  )
}

message("  Found ", nrow(all_models), " significant model(s) to analyse:")
for (i in seq_len(nrow(all_models))) {
  message(sprintf("    [%d] %-28s @ %-14s  p=%.4g",
                  i, all_models$metrics[i], all_models$group_value[i], all_models$p_num[i]))
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

# --- Training data for one hit: the screen model's own frame ---
# Same rows, QC and z-scaling as the bam() the hit came from. Factor columns
# are re-levelled (droplevels) so subsets behave in extend() / powerCurve().
training_frame <- function(mod_bam) {
  df <- as.data.frame(mod_bam$model)
  for (v in intersect(c("global_subjid", "experiment_id", "sex", "cap_size"), names(df)))
    df[[v]] <- factor(as.character(df[[v]]))
  droplevels(df)
}

# --- Refit as lmer for simr ---
# metric_vars: the metric column(s) as in the bam (one pooled column, or the
# between + within pair). Entered linearly: simr needs a parametric effect.
# Covariates mirror the screen's base model; a factor with one level in this
# subset is left out (contrasts would fail).
fit_lmer_for_simr <- function(df, metric_vars) {
  fixed <- c("laser_power_z", "trial_index_z", "age_z")
  fixed <- fixed[fixed %in% names(df)]
  for (v in c("sex", "cap_size"))
    if (v %in% names(df) && nlevels(droplevels(df[[v]])) > 1) fixed <- c(fixed, v)
  re_terms <- "(1 | global_subjid)"
  if ("experiment_id" %in% names(df) && n_distinct(df$experiment_id) > 1)
    re_terms <- paste(re_terms, "+ (1 | experiment_id)")
  f <- as.formula(paste(
    "pain_rating ~", paste(c(fixed, metric_vars), collapse = " + "), "+", re_terms
  ))
  tryCatch(
    lmer(f, data = df, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { message("  [WARN] lmer failed: ", conditionMessage(e)); NULL }
  )
}

# --- The fitted lmer with the metric's coefficient multiplied by `scale` ---
scale_effect <- function(mod, test_var, scale) {
  if (isTRUE(all.equal(scale, 1))) return(mod)
  fixef(mod)[test_var] <- fixef(mod)[test_var] * scale   # simr's fixef<-
  mod
}

# --- Run power curve for one subset (pooled or per cap_size) ---
# test_var: the effect whose power is estimated. Without `test =`, simr
# tests the FIRST fixed effect in the formula, which here is laser_power_z.
run_power_curve <- function(mod, along_var, n_current, n_target, label, n_steps, n_sim,
                            test_var, effect_scale) {
  breaks <- unique(round(seq(n_current, n_target, length.out = n_steps)))
  breaks <- breaks[breaks >= n_current]
  tryCatch(
    {
      mod_s <- scale_effect(mod, test_var, effect_scale)
      # powerCurve() can only subset subjects that exist; extend() adds
      # simulated subjects so the curve can reach the projected N.
      mod_ext <- if (n_target > n_current) simr::extend(mod_s, along = along_var, n = n_target) else mod_s
      pc <- powerCurve(mod_ext, test = simr::fixed(test_var, "z"), along = along_var,
                       breaks = breaks, nsim = n_sim, progress = FALSE)
      # summary() is the documented accessor (as.data.frame() has no method
      # for powerCurve objects): one row per break.
      sm <- summary(pc)
      df <- data.frame(n = sm$nlevels, power = sm$mean, lower = sm$lower, upper = sm$upper)
      df$group <- label
      df$effect_scale <- effect_scale
      df
    },
    error = function(e) {
      message("  [WARN] powerCurve failed (", label, ", effect x", effect_scale, "): ",
              conditionMessage(e))
      NULL
    }
  )
}

# --- Extract N required to hit each power target ---
# share: this group's share of the training subjects (1 for pooled), so a
# per-cap-size requirement can be expressed as total training N.
extract_n_required <- function(curve_df, label, power_targets, share, n_current) {
  if (is.null(curve_df)) return(NULL)
  map_dfr(power_targets, function(thresh) {
    hit <- curve_df %>% filter(power >= thresh) %>% slice(1)
    n_req <- if (nrow(hit) > 0) as.integer(hit$n) else NA_integer_
    tibble(
      group            = label,
      effect_scale     = curve_df$effect_scale[1],
      power_target     = thresh,
      n_current        = as.integer(n_current),
      n_required       = n_req,
      n_required_total = as.integer(ceiling(n_req / share)),
      achieved         = nrow(hit) > 0,
      max_power        = max(curve_df$power)
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
  if (n_distinct(df_ho$global_subjid) < 3) return(NULL)

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

  ci_width_r <- unname(diff(quantile(boot_df$r, c(0.025, 0.975), na.rm = TRUE)))
  n_ho_subj  <- n_distinct(df_ho$global_subjid)

  list(
    n_holdout   = n_ho_subj,
    rmse_mean   = mean(boot_df$rmse),
    rmse_ci_lo  = unname(quantile(boot_df$rmse, 0.025)),
    rmse_ci_hi  = unname(quantile(boot_df$rmse, 0.975)),
    r_mean      = mean(boot_df$r, na.rm = TRUE),
    r_ci_lo     = unname(quantile(boot_df$r, 0.025, na.rm = TRUE)),
    r_ci_hi     = unname(quantile(boot_df$r, 0.975, na.rm = TRUE)),
    ci_width_r  = ci_width_r,
    stable      = ci_width_r <= ci_flag
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

  message(sprintf("\n[%d/%d] %s @ %s  (p=%.4g)",
                  i, nrow(all_models), this_metric, this_roi, row$p_num))

  df_train    <- training_frame(mod_bam)
  metric_vars <- model_metric_vars(mod_bam, this_metric)
  test_var    <- primary_var(mod_bam, this_metric)
  if (nrow(df_train) == 0 || !all(metric_vars %in% names(df_train))) {
    message("  Skipping — model frame lacks the metric columns.")
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

  power_rows <- NULL
  all_curves_this_model <- list()

  for (es in pa_cfg$effect_scales) {
    # --- Pooled power curve ---
    curve_pooled <- run_power_curve(
      mod_lmer, "global_subjid",
      n_current = n_train_subj,
      n_target  = pa_cfg$n_target_total,
      label     = "pooled",
      n_steps   = pa_cfg$n_steps,
      n_sim     = pa_cfg$n_sim,
      test_var  = test_var,
      effect_scale = es
    )
    power_rows <- bind_rows(power_rows, extract_n_required(curve_pooled, "pooled",
                                                           pa_cfg$power_targets, 1, n_train_subj))
    all_curves_this_model[[length(all_curves_this_model) + 1]] <- curve_pooled

    # --- Per cap_size power curves ---
    for (this_cap in sort(unique(as.character(df_train$cap_size)))) {
      df_cap   <- droplevels(df_train[as.character(df_train$cap_size) == this_cap, ])  # unused levels break extend()
      n_cap    <- n_distinct(df_cap$global_subjid)
      if (n_cap < 5) {
        message("  Skipping cap_size ", this_cap, " — only ", n_cap, " subjects.")
        next
      }

      mod_cap <- fit_lmer_for_simr(df_cap, metric_vars)
      if (is.null(mod_cap)) next

      share        <- n_cap / n_train_subj
      n_target_cap <- round(pa_cfg$n_target_total * share)
      lbl <- paste0("cap_size_", this_cap)

      curve_cap <- run_power_curve(
        mod_cap, "global_subjid",
        n_current = n_cap,
        n_target  = max(n_target_cap, n_cap + 5L),
        label     = lbl,
        n_steps   = pa_cfg$n_steps,
        n_sim     = pa_cfg$n_sim,
        test_var  = test_var,
        effect_scale = es
      )
      power_rows <- bind_rows(power_rows, extract_n_required(curve_cap, lbl, pa_cfg$power_targets,
                                                             share, n_cap))
      all_curves_this_model[[length(all_curves_this_model) + 1]] <- curve_cap
    }
  }

  # Tag power rows with metric/ROI identity
  if (!is.null(power_rows) && nrow(power_rows) > 0) {
    power_rows <- power_rows %>%
      mutate(metric = this_metric, roi = this_roi, fitted_at = as.character(Sys.time()))
    all_power_rows[[length(all_power_rows) + 1]] <- power_rows

    message("  Power summary:")
    for (j in seq_len(nrow(power_rows))) {
      achieved_str <- if (power_rows$achieved[j]) "ACHIEVED" else "NOT ACHIEVED in range"
      message(sprintf("    %-14s effect x%.2f  power=%.0f%%  n_required=%-5s (total training N %-5s)  %s",
                      power_rows$group[j], power_rows$effect_scale[j],
                      power_rows$power_target[j] * 100,
                      ifelse(is.na(power_rows$n_required[j]), ">max", power_rows$n_required[j]),
                      ifelse(is.na(power_rows$n_required_total[j]), ">max", power_rows$n_required_total[j]),
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
              ") — consider a larger holdout share (split_train_frac in merge_defaults.R).")
    }

    ho_row <- tibble(
      metric        = this_metric,
      roi           = this_roi,
      group         = "holdout_bootstrap",
      n_current     = as.integer(ho$n_holdout),
      holdout_r     = ho$r_mean,
      holdout_r_lo  = ho$r_ci_lo,
      holdout_r_hi  = ho$r_ci_hi,
      holdout_rmse  = ho$rmse_mean,
      holdout_stable = ho$stable,
      fitted_at     = as.character(Sys.time())
    )
    all_power_rows[[length(all_power_rows) + 1]] <- ho_row
  }

  # --- Power curve figure ---
  plot_data <- compact(all_curves_this_model) %>% bind_rows()
  if (nrow(plot_data) > 0) {
    plot_data$effect <- factor(sprintf("effect x%.2g", plot_data$effect_scale),
                               levels = sprintf("effect x%.2g", sort(unique(plot_data$effect_scale), TRUE)))
    p <- ggplot(plot_data, aes(x = n, y = power, colour = group, fill = group)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      geom_hline(yintercept = pa_cfg$power_targets, linetype = "dashed",
                 colour = "grey40", linewidth = 0.5) +
      facet_wrap(~ effect, nrow = 1) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      labs(
        title    = paste0("Power curve: ", metric_label(this_metric), ", ", unit_label(this_roi)),
        subtitle = paste0("simr (nsim=", pa_cfg$n_sim, " per point) | lmer approximation of bam() | ",
                          "screen p ", fmt_p(row$p_num), " | power for ", test_var),
        x        = "N subjects in the group (training set)",
        y        = "Estimated power",
        colour   = NULL, fill = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", plot.background = element_rect(fill = "white", colour = NA))

    fig_name <- paste0("power_curve_",
                       str_replace_all(this_metric, "_z$", ""), "_",
                       str_replace_all(this_roi, "-", "_"),
                       if (decomp == "pooled") "_pooled" else "", ".png")
    fig_path <- file.path(pa_cfg$fig_dir, fig_name)
    ggsave(fig_path, p, width = 4 + 4 * length(pa_cfg$effect_scales), height = 5, dpi = 150)
    message("  Figure saved: ", fig_path)
  }
}

# =============================================================================
# STEP 4: BINDING CONSTRAINT SUMMARY
# =============================================================================
message("\n", strrep("=", 60))
message("== BINDING CONSTRAINT SUMMARY ==")
message(strrep("=", 60))

all_rows_df  <- bind_rows(all_power_rows)
all_power_df <- if (nrow(all_rows_df) > 0 && "power_target" %in% names(all_rows_df))
  filter(all_rows_df, !is.na(power_target)) else tibble()   # exclude holdout_bootstrap rows

if (nrow(all_power_df) > 0) {
  verdict_scale <- min(pa_cfg$effect_scales)

  for (es in sort(unique(all_power_df$effect_scale), decreasing = TRUE)) {
    for (thresh in pa_cfg$power_targets) {
      message(sprintf("\nFor %.0f%% power, effect x%.2f:", thresh * 100, es))
      sub <- all_power_df %>%
        filter(power_target == thresh, effect_scale == es) %>%
        arrange(desc(n_required_total))
      if (nrow(sub) == 0) { message("  No data."); next }

      binding <- sub %>% filter(!is.na(n_required_total)) %>% slice(1)
      if (nrow(binding) > 0) {
        message(sprintf("  Binding constraint: %s @ %s (%s)  total training N required = %d",
                        binding$metric, binding$roi, binding$group, binding$n_required_total))
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
  }

  # Verdict (most conservative effect scale, 80% power)
  message(sprintf("\nVerdict (80%% power, effect x%.2f):", verdict_scale))
  v80 <- all_power_df %>% filter(power_target == 0.80, effect_scale == verdict_scale)
  if (any(!v80$achieved)) {
    message("  [ACTION REQUIRED] ", sum(!v80$achieved), " metric x ROI x group curve(s) never reach ",
            "80% power within the projected N, even with every subject in training.")
  }
  binding_80 <- v80 %>% filter(!is.na(n_required_total)) %>%
    arrange(desc(n_required_total)) %>% slice(1)

  if (nrow(binding_80) > 0) {
    current_train <- sum(split_summary$n[split_summary$split_group == "train"])
    if (binding_80$n_required_total <= current_train) {
      message(sprintf("  Current training N (%d) is SUFFICIENT for 80%% power", current_train))
      message(sprintf("  across all models that reach it (binding constraint: total N = %d).",
                      binding_80$n_required_total))
    } else {
      message(sprintf("  [ACTION REQUIRED] Current training N (%d) is INSUFFICIENT.", current_train))
      message(sprintf("  Need total training N = %d for 80%% power on the binding metric (%s @ %s, %s).",
                      binding_80$n_required_total, binding_80$metric, binding_80$roi, binding_80$group))
      new_frac  <- min(1, round(binding_80$n_required_total / pa_cfg$n_target_total + 0.01, 2))
      cur_frac  <- n_train / (n_train + n_holdout)
      if (new_frac <= cur_frac) {
        message(sprintf("  At full N (%d) that needs a training fraction of %.2f; the current split (%.2f) ",
                        pa_cfg$n_target_total, new_frac, cur_frac),
                "already meets it once recruitment is complete - no change to the split needed.")
      } else {
        message(sprintf("  At full N (%d) that needs a training fraction of at least %.2f, above the current %.2f: ",
                        pa_cfg$n_target_total, new_frac, cur_frac),
                "raise split_train_frac in merge_defaults.R (applies to subjects assigned from now on).")
      }
    }
  }
}

# =============================================================================
# STEP 5: WRITE ALL RESULTS TO results.duckdb
# =============================================================================
message("\n== Step 5: Writing to results.duckdb ==")

# An older layout (no feature_set column) is kept under another name, not
# dropped: power curves are expensive to recompute.
if (dbExistsTable(results_con, "power_analysis") &&
    !"feature_set" %in% dbListFields(results_con, "power_analysis")) {
  old_name <- "power_analysis_pre_v3"
  if (dbExistsTable(results_con, old_name)) old_name <- paste0(old_name, "_", format(Sys.time(), "%Y%m%d%H%M%S"))
  invisible(dbExecute(results_con, sprintf("ALTER TABLE power_analysis RENAME TO %s", old_name)))
  message("  Older power_analysis table renamed to ", old_name, ".")
}

invisible(dbExecute(results_con, "
  CREATE TABLE IF NOT EXISTS power_analysis (
    fitted_at        VARCHAR,
    feature_set      VARCHAR,
    decomposition    VARCHAR,
    metric           VARCHAR,
    roi              VARCHAR,
    group_label      VARCHAR,
    effect_scale     DOUBLE,
    power_target     DOUBLE,
    n_current        INTEGER,
    n_required       INTEGER,
    n_required_total INTEGER,
    achieved         BOOLEAN,
    max_power        DOUBLE,
    holdout_r        DOUBLE,
    holdout_r_lo     DOUBLE,
    holdout_r_hi     DOUBLE,
    holdout_rmse     DOUBLE,
    holdout_stable   BOOLEAN
  )
"))

final_rows <- all_rows_df
if (nrow(final_rows) > 0) {
  final_rows <- final_rows %>%
    rename(group_label = group) %>%
    mutate(feature_set = gamm_cfg$feature_set, decomposition = decomp)
  cols <- dbListFields(results_con, "power_analysis")
  for (cn in setdiff(cols, names(final_rows))) final_rows[[cn]] <- NA
  final_rows <- as.data.frame(final_rows[, cols])

  # Upsert all rows for this run (delete this metric x ROI, then insert).
  # No primary key: holdout_bootstrap rows have power_target = NULL, which a
  # PRIMARY KEY column cannot hold.
  keys <- distinct(final_rows, metric, roi)
  for (k in seq_len(nrow(keys))) {
    invisible(dbExecute(results_con,
      "DELETE FROM power_analysis WHERE feature_set = ? AND decomposition = ? AND metric = ? AND roi = ?",
      params = list(gamm_cfg$feature_set, decomp, keys$metric[k], keys$roi[k])
    ))
  }
  dbAppendTable(results_con, "power_analysis", final_rows)
  message("  Written ", nrow(final_rows), " rows to power_analysis table.")
}

disconnect_db(data_con)
close_results_db(results_con)

message("\nPower analysis complete.")
message("Results: ", gamm_cfg$results_db_path, "  (table: power_analysis)")
message("Figures: ", pa_cfg$fig_dir)