# =============================================================================
# robustness_core.R
# =============================================================================
# Robustness checks for GAMM hits. Run from R/gamm/ after source_core.R /
# channel_core.R and holdout_replicate.R:
#   setwd("R/gamm"); source("../robustness/robustness_core.R")
#
# Which hits: single-metric models of the current feature_set and
# decomposition with p < hits_p_sig, at gamm_cfg$robustness$levels; with
# only_replicated = TRUE (default) only those that replicated in
# holdout_replication. Which tests: gamm_cfg$robustness$tests.
#
# Every prediction test refits the screen model AND its nested baseline
# (metric terms removed) on the same subset, predicts people the refit never
# saw at the population level (subject random effect excluded; experiment
# effect kept only for experiments in the refit), and scores the change in
# prediction error, metric - baseline (negative = the metric helps), exactly
# as holdout_replicate.R does. Data are the screen model's own rows
# (mod$model: same QC, same z-scaling) and held-out rows prepared with
# prepare_holdout() (training scaling).
#
#   loeo          Leave one EXPERIMENT out. Refit on the training subjects of
#                 every other experiment; predict all subjects (training and
#                 held-out) of the left-out experiment. One row per
#                 experiment (variant "exp_<id>") + variant "all" (all
#                 left-out predictions pooled). Tests transfer to new
#                 hardware / paradigm / lab. If the left-out experiment's
#                 cap size occurs in no other experiment, that fold drops the
#                 cap_size term from both models (noted). Caveat: training
#                 subjects of the left-out experiment helped SELECT the hit
#                 in the screen.
#   kfold         Grouped k-fold CV inside the training set (whole subjects
#                 per fold, folds stratified by cap size). Rows "fold_01".. +
#                 "all". Generalisation to new people from the same
#                 experiments, without touching the held-out set.
#   permutation   Screen model and baseline fixed. On the held-out set, the
#                 metric's primary part is shuffled at SUBJECT level (each
#                 person's mean given to another person of the same
#                 experiment; within-person deviations kept), or within
#                 subject for primary_part = "within". stat = observed
#                 primary change in error (%), stat_lo/hi = 2.5/97.5% of the
#                 null, p_value = (1 + #null <= observed) / (n_perm + 1).
#                 Null: the metric carries no subject-specific information.
#   split_half    Random halves of the training subjects (stratified by cap
#                 size), metric model refit on each; stat = median
#                 correlation of the two fitted metric curves (primary part,
#                 1st-99th percentile grid), stat_lo/hi = 2.5/97.5% over
#                 repetitions. Stability of the functional form. Curves
#                 shrunk to flat (edf ~ 0) give no correlation; counted in
#                 note.
#   fooof_sensitivity  Drop training and held-out subjects whose
#                 subject-level fooof_r2 at this unit is below each
#                 threshold (variant "fooof_r2_<thr>"), refit, predict the
#                 remaining held-out subjects.
#
# For prediction tests stat = primary change in error (%), stat_lo/hi its
# 95% CI and p_value its one-sided bootstrap p (primary level from
# gamm_cfg$holdout_primary); trial_* and subj_* hold both levels.
#
# Writes: results.duckdb  robustness (one row per hit x test x variant)
#         <fig_dir>/<level>/<unit>/[pooled/]<metric>_robustness.png
#
# V1 (Oct 2026): replaces run_robustness.R, run_cv_internal.R,
# run_cv_external.R, cv_helpers.R and cv_split.R, which never ran (wrong
# paths and packages, undefined objects, trial-level permutation,
# holdout-scaled z-scores, no between/within columns). External validation
# lives in R/gamm/holdout_replicate.R; leave-one-subject-out is replaced by
# grouped k-fold (LOSO is ~N refits per hit at full N).
# =============================================================================

library(dplyr)
library(tibble)
library(mgcv)
library(ggplot2)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/hit_helpers.R")
source("../robustness/helpers/robustness_helpers.R")

cfg     <- gamm_cfg
rc      <- cfg$robustness
decomp  <- cfg$metric_decomposition
primary <- match.arg(cfg$holdout_primary, c("subject", "trial"))
tests   <- match.arg(rc$tests, c("loeo", "kfold", "permutation", "split_half", "fooof_sensitivity"),
                     several.ok = TRUE)
set.seed(rc$seed)

data_con    <- connect_data_db(cfg$data_db_path, read_only = TRUE)
results_con <- init_results_db(cfg$results_db_path)
init_robustness_table(results_con)

# =============================================================================
# Hits
# =============================================================================
hits <- load_hits(results_con, cfg$feature_set, p_sig = cfg$hits_p_sig,
                  levels = rc$levels, decomposition = decomp) %>% filter(rds_ok)
if (isTRUE(rc$only_replicated)) {
  if (!dbExistsTable(results_con, "holdout_replication") ||
      !"decomposition" %in% dbListFields(results_con, "holdout_replication"))
    stop("robustness$only_replicated = TRUE but there is no holdout_replication table for this ",
         "version: run holdout_replicate.R first (or set only_replicated = FALSE).")
  rep_rows <- dbGetQuery(results_con, "
    SELECT level, group_value, metric AS metrics FROM holdout_replication
    WHERE feature_set = ? AND decomposition = ? AND replicated",
    params = list(cfg$feature_set, decomp))
  hits <- semi_join(hits, rep_rows, by = c("level", "group_value", "metrics"))
}
if (nrow(hits) == 0) stop("No hits to test (feature_set ", cfg$feature_set, ", decomposition ",
                          decomp, if (isTRUE(rc$only_replicated)) ", replicated only" else "", ").")
message("Robustness: ", nrow(hits), " hit(s); tests: ", paste(tests, collapse = ", "),
        "; primary criterion = ", primary, "-level prediction")

pcol <- function(lvl) if (lvl == "subject") c("subj_pct", "subj_pct_lo", "subj_pct_hi", "subj_p") else
  c("trial_pct", "trial_pct_lo", "trial_pct_hi", "trial_p")

row_id <- function(h, test, variant) {
  tibble(fitted_at = as.character(Sys.time()), feature_set = cfg$feature_set,
         decomposition = decomp, level = h$level, group_value = h$group_value,
         metric = h$metrics, test = test, variant = variant)
}
with_scores <- function(id, sc, n_train, note = NA_character_) {
  if (is.null(sc)) return(mutate(id, status = "no_predictions", n_train_subj = n_train, note = note))
  pc <- pcol(primary)
  bind_cols(id, sc) %>%
    mutate(status = "ok", n_train_subj = as.integer(n_train),
           stat = .data[[pc[1]]], stat_lo = .data[[pc[2]]], stat_hi = .data[[pc[3]]],
           p_value = .data[[pc[4]]], note = note)
}
status_row <- function(id, status, note = NA_character_) mutate(id, status = status, note = note)
n_subj <- function(df) n_distinct(as.character(df$global_subjid))

# =============================================================================
# Tests
# =============================================================================
test_loeo <- function(h, mod, train, ho) {
  exps <- sort(unique(c(as.character(train$experiment_id), as.character(ho$experiment_id))))
  rows <- list(); pooled <- list()
  for (E in exps) {
    id <- row_id(h, "loeo", paste0("exp_", E))
    tr <- train[as.character(train$experiment_id) != E, ]
    te <- bind_rows(to_pred_frame(train[as.character(train$experiment_id) == E, ]),
                    to_pred_frame(ho[as.character(ho$experiment_id) == E, names(train)]))
    if (n_subj(te) < rc$min_test_subjects) { rows[[E]] <- status_row(id, "too_few_test_subjects"); next }
    if (n_subj(tr) < rc$min_test_subjects) { rows[[E]] <- status_row(id, "too_few_train_subjects"); next }
    # A cap size that occurs only in the left-out experiment has no estimate
    # in the refit: drop the cap_size term from both models for this fold.
    cap_new <- !all(unique(te$cap_size) %in% as.character(tr$cap_size))
    fits <- tryCatch(refit_pair(mod, tr, h$metrics, cfg, drop_vars = if (cap_new) "cap_size" else character(0)),
                     error = function(e) e)
    if (inherits(fits, "error")) { rows[[E]] <- status_row(id, "fit_failed", conditionMessage(fits)); next }
    pm <- predict_population(fits$metric, te); pb <- predict_population(fits$base, te)
    rows[[E]] <- with_scores(id, score_predictions(te, pm, pb, rc$n_boot), n_subj(tr),
                             if (cap_new) paste0("cap size ", paste(unique(te$cap_size), collapse = "/"),
                                                 " absent from the refit: cap_size term dropped") else NA_character_)
    pooled[[E]] <- list(te = te, pm = pm, pb = pb)
  }
  if (length(pooled) > 0) {
    te <- bind_rows(lapply(pooled, `[[`, "te"))
    rows[["all"]] <- with_scores(row_id(h, "loeo", "all"),
                                 score_predictions(te, unlist(lapply(pooled, `[[`, "pm")),
                                                   unlist(lapply(pooled, `[[`, "pb")), rc$n_boot),
                                 NA_integer_, sprintf("%d experiment(s) predicted", length(pooled)))
  }
  bind_rows(rows)
}

test_kfold <- function(h, mod, train) {
  sc_df <- distinct(to_pred_frame(train)[, c("global_subjid", "cap_size")])
  folds <- stratified_folds(sc_df, rc$kfold_k, rc$seed)
  f_row <- folds[as.character(train$global_subjid)]
  rows <- list(); pooled <- list()
  for (k in sort(unique(folds))) {
    id <- row_id(h, "kfold", sprintf("fold_%02d", k))
    tr <- train[f_row != k, ]; te <- to_pred_frame(train[f_row == k, ])
    fits <- tryCatch(refit_pair(mod, tr, h$metrics, cfg), error = function(e) e)
    if (inherits(fits, "error")) { rows[[k]] <- status_row(id, "fit_failed", conditionMessage(fits)); next }
    pm <- predict_population(fits$metric, te); pb <- predict_population(fits$base, te)
    rows[[k]] <- with_scores(id, score_predictions(te, pm, pb, rc$n_boot), n_subj(tr))
    pooled[[k]] <- list(te = te, pm = pm, pb = pb)
  }
  pooled <- Filter(Negate(is.null), pooled)
  if (length(pooled) > 0) {
    te <- bind_rows(lapply(pooled, `[[`, "te"))
    rows[[length(rows) + 1]] <- with_scores(
      row_id(h, "kfold", "all"),
      score_predictions(te, unlist(lapply(pooled, `[[`, "pm")), unlist(lapply(pooled, `[[`, "pb")), rc$n_boot),
      n_subj(train), sprintf("%d-fold, whole subjects, stratified by cap size", rc$kfold_k))
  }
  bind_rows(rows)
}

test_permutation <- function(h, mod, ho) {
  id <- row_id(h, "permutation", "subject_level")
  if (n_subj(ho) < rc$min_test_subjects) return(status_row(id, "too_few_test_subjects"))
  base <- tryCatch(refit_baseline(mod, h$metrics, cfg), error = function(e) NULL)
  if (is.null(base)) return(status_row(id, "baseline_fit_failed"))
  pb <- predict_population(base, ho)
  obs <- primary_pct(ho, predict_population(mod, ho), pb, primary)
  part <- if (identical(decomp, "between_within")) mod$gamm_meta$primary_part else "pooled"
  null <- vapply(seq_len(rc$n_perm), function(i) {
    hp <- permute_primary(ho, h$metrics, decomp, part, isTRUE(rc$perm_within_experiment))
    primary_pct(hp, predict_population(mod, hp), pb, primary)
  }, numeric(1))
  mutate(id, status = "ok", n_test_subj = n_subj(ho), n_test_trials = nrow(ho),
         stat = obs, stat_lo = unname(quantile(null, 0.025, na.rm = TRUE)),
         stat_hi = unname(quantile(null, 0.975, na.rm = TRUE)),
         p_value = (1 + sum(null <= obs, na.rm = TRUE)) / (sum(!is.na(null)) + 1),
         note = sprintf("%d permutations of the %s part (%s), %s-level error", rc$n_perm, part,
                        if (part == "within") "within subject" else
                          if (isTRUE(rc$perm_within_experiment)) "subjects within experiment" else "subjects",
                        primary))
}

test_split_half <- function(h, mod, train) {
  id <- row_id(h, "split_half", "all")
  pv <- primary_var(mod, h$metrics)
  grid <- setNames(data.frame(metric_grid(mod, pv, n = 100)), pv)
  sc_df <- distinct(to_pred_frame(train)[, c("global_subjid", "cap_size")])
  rs <- vapply(seq_len(rc$n_split_half), function(r) {
    half <- stratified_folds(sc_df, 2L, rc$seed + r)[as.character(train$global_subjid)]
    curves <- lapply(1:2, function(k) {
      fit <- tryCatch(refit_pair(mod, train[half == k, ], h$metrics, cfg, baseline = FALSE)$metric,
                      error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      predict_from_lp(fit, marginal_lpmatrix(fit, grid))$fit
    })
    if (any(vapply(curves, is.null, logical(1)))) return(NA_real_)
    if (stats::sd(curves[[1]]) < 1e-8 || stats::sd(curves[[2]]) < 1e-8) return(NA_real_)
    stats::cor(curves[[1]], curves[[2]])
  }, numeric(1))
  if (all(is.na(rs))) return(status_row(id, "no_curves", "every half-fit shrank the term to flat or failed"))
  mutate(id, status = "ok", n_train_subj = as.integer(n_subj(train) %/% 2),
         stat = median(rs, na.rm = TRUE),
         stat_lo = unname(quantile(rs, 0.025, na.rm = TRUE)),
         stat_hi = unname(quantile(rs, 0.975, na.rm = TRUE)),
         note = sprintf("%d repetitions; %d with a flat or failed half", length(rs), sum(is.na(rs))))
}

test_fooof <- function(h, mod, train, ho) {
  unit_col <- level_spec(h$level)$unit_col
  fq <- bind_rows(load_subject_metrics_wide(data_con, h$level, unit = h$group_value,
                                            metrics = "fooof_r2", split = "train"),
                  load_subject_metrics_wide(data_con, h$level, unit = h$group_value,
                                            metrics = "fooof_r2", split = "holdout"))
  if (!"fooof_r2" %in% names(fq))
    return(status_row(row_id(h, "fooof_sensitivity", "all"), "no_fooof_r2"))
  r2 <- setNames(fq$fooof_r2, fq$global_subjid)
  bind_rows(lapply(rc$fooof_thresholds, function(thr) {
    id <- row_id(h, "fooof_sensitivity", sprintf("fooof_r2_%.2f", thr))
    keep <- function(df) { v <- r2[as.character(df$global_subjid)]; df[!is.na(v) & v >= thr, ] }
    tr <- keep(train); te <- keep(ho)
    note <- sprintf("kept %d of %d training and %d of %d held-out subjects", n_subj(tr), n_subj(train),
                    n_subj(te), n_subj(ho))
    if (n_subj(tr) < rc$min_test_subjects || n_subj(te) < rc$min_test_subjects)
      return(status_row(id, "too_few_subjects", note))
    fits <- tryCatch(refit_pair(mod, tr, h$metrics, cfg), error = function(e) e)
    if (inherits(fits, "error")) return(status_row(id, "fit_failed", conditionMessage(fits)))
    with_scores(id, score_predictions(te, predict_population(fits$metric, te),
                                      predict_population(fits$base, te), rc$n_boot),
                n_subj(tr), note)
  }))
}

# =============================================================================
# Figure: primary change in error per test / variant for one hit
# =============================================================================
plot_robustness <- function(h, res) {
  fd <- res %>% filter(status == "ok", test %in% c("loeo", "kfold", "fooof_sensitivity"),
                       !(test == "kfold" & variant != "all")) %>%
    mutate(lab = paste(test, variant, sep = ": "),
           lab = factor(lab, levels = rev(unique(lab))))
  if (nrow(fd) == 0) return(NULL)
  sub <- c(
    if (any(res$test == "permutation" & res$status == "ok"))
      with(res[res$test == "permutation", ], sprintf("permutation p %s", fmt_p(p_value))),
    if (any(res$test == "split_half" & res$status == "ok"))
      with(res[res$test == "split_half", ], sprintf("split-half curve r = %.2f [%.2f, %.2f]", stat, stat_lo, stat_hi)))
  ggplot(fd, aes(x = stat, y = lab)) +
    geom_vline(xintercept = 0, colour = INK_3, linewidth = 0.4, linetype = "dashed") +
    geom_linerange(aes(xmin = stat_lo, xmax = stat_hi), colour = cfg$col_main, linewidth = 0.6) +
    geom_point(colour = cfg$col_main, size = 2.4) +
    labs(title = sprintf("Robustness: %s, %s", metric_label(h$metrics), unit_label(h$group_value)),
         subtitle = paste(sub, collapse = " | "),
         x = sprintf("Change in %s-level prediction error, metric vs baseline (%%)", primary), y = NULL,
         caption = "Negative = the metric improves prediction for people the refit never saw. 95% CI from resampling test subjects.") +
    theme_gamm()
}

# =============================================================================
# Run
# =============================================================================
for (i in seq_len(nrow(hits))) {
  h <- hits[i, ]
  message(sprintf("\n[%d/%d] %-8s %-12s %s", i, nrow(hits), h$level, h$group_value, h$metrics))
  mod   <- readRDS(h$rds_path)
  train <- as.data.frame(mod$model)
  ho    <- prepare_holdout(data_con, h$level, h$group_value, h$metrics, mod, cfg)
  if (is.null(ho)) { message("  [WARN] model has no $zscale - re-run the screen."); next }
  ho    <- to_pred_frame(ho)

  res <- list()
  run <- function(name, expr) {
    if (!name %in% tests) return(NULL)
    message("  ", name, " ...")
    t0 <- Sys.time()
    out <- tryCatch(expr, error = function(e) status_row(row_id(h, name, "all"), "error", conditionMessage(e)))
    message(sprintf("    done (%.0f s)", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    out
  }
  res$loeo  <- run("loeo",  test_loeo(h, mod, train, ho))
  res$kfold <- run("kfold", test_kfold(h, mod, train))
  res$perm  <- run("permutation", test_permutation(h, mod, ho))
  res$half  <- run("split_half", test_split_half(h, mod, train))
  res$fooof <- run("fooof_sensitivity", test_fooof(h, mod, train, ho))
  res <- bind_rows(res)
  write_robustness_rows(results_con, res)

  for (j in seq_len(nrow(res))) {
    r <- res[j, ]
    message(sprintf("    %-18s %-16s %-22s stat %+8.2f  [%+.2f, %+.2f]  p %s",
                    r$test, r$variant, r$status, r$stat %||% NA, r$stat_lo %||% NA, r$stat_hi %||% NA,
                    fmt_p(r$p_value %||% NA)))
  }

  p <- tryCatch(plot_robustness(h, res), error = function(e) NULL)
  if (!is.null(p)) {
    out_dir <- fig_unit_dir(cfg$fig_dir, h$level, h$group_value, decomp)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(out_dir, paste0(sub("_z$", "", h$metrics), "_robustness.png")), p,
           width = cfg$fig_width * 1.3, height = max(3, 1.8 + 0.3 * nrow(p$data)),
           dpi = cfg$fig_dpi, bg = "white")
  }
}

disconnect_db(data_con)
close_results_db(results_con)
message("\nTable: ", cfg$results_db_path, " (robustness)")