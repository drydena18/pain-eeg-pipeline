# =============================================================================
# holdout_replicate.R
# =============================================================================
# Does each GAMM hit PREDICT pain in people the model has never seen?
# Run from R/gamm/ after source_core.R / channel_core.R.
#
# This is the replication step the uncorrected screen relies on. The
# criterion is prediction (research question: can a candidate biomarker
# predict pain sensitivity?), not significance in a refit.
#
# gamm_cfg$holdout_part chooses WHICH part of a between/within-split metric
# is tested:
#
# holdout_part = "between" (original behaviour)
#   For every single-metric model whose BETWEEN-subject part has p <
#   hits_p_sig (current feature_set and decomposition):
#   1. Metric model   = the screen model (fit on TRAIN subjects).
#      Baseline model = the same model without the metric term, refit on
#      exactly the same training rows (nested; differs only by the metric).
#   2. Held-out subjects' trials for that unit, with the screen's QC, are
#      z-scored with the TRAINING mean/SD.
#   3. Both models predict every held-out trial at the population level:
#      subject random effect excluded (new people); experiment effect kept
#      for experiments seen in training.
#   4. Prediction error, metric vs baseline, at two levels:
#        trial   : mean squared error over trials
#        subject : squared error of each subject's MEAN rating
#      Change in MSE = metric - baseline (negative = metric predicts better),
#      with a 95% CI from resampling held-out SUBJECTS (2000 draws) and a
#      one-sided bootstrap p for "no improvement".
#   5. Primary criterion (holdout_primary, default "subject"): replicates if
#      the primary error is at least holdout_min_pct % lower than baseline
#      AND its Benjamini-Hochberg-adjusted p across all hits is < holdout_alpha.
#
# holdout_part = "within" (V4)
#   Hits = models whose WITHIN-subject part has p < hits_p_sig, or exactly the
#   pre-registered holdout_targets (tested whatever their screen p).
#   1. Baseline = the screen model WITHOUT the within term only (between term
#      kept), refit on the same training rows: the comparison isolates the
#      within-person part.
#   2-3. As above (training z-scale; each held-out person's between / within
#      split from their own trials; population-level predictions).
#   4. Primary error = SUBJECT-CENTRED trial error: each held-out person's
#      own mean is removed from the observed ratings and from both models'
#      predictions, so only the ordering of that person's trials is scored.
#      Change in error, subject-resampling CI and one-sided p as above.
#   5. The same model refit on the held-out subjects ALONE: p, edf and
#      direction of the within term in an independent sample.
#   6. Criterion (holdout_within_criterion):
#        "refit" (default): replicates if the within term has BH-adjusted
#          p < holdout_alpha in the held-out refit AND the same direction as
#          in training.
#        "prediction": replicates if the subject-centred error is lower than
#          baseline by at least holdout_min_pct_within % AND BH p < alpha.
#      Why refit by default: a within-person effect of ~0.2 rating points
#      lowers trial error by only ~0.1 %, so the prediction test has little
#      power (simulation, 53 trials/subject, residual SD 1.78: ~25 % power
#      with 40 held-out subjects, ~45 % with 128; refit ~38 % and ~90 %;
#      both ~5-7 % false positives under the null).
#   BH runs across the tested hits (a single pre-registered target = its own p).
#
# Writes: results.duckdb  holdout_replication (one row per hit and tested part)
#         <fig_dir>/holdout_replication_<level>[_pooled][_within].png
# Reads:  data.duckdb (holdout subjects, read-only), gamm_fitted, the .rds
#
# V4 (Nov 2026): holdout_part / holdout_targets / holdout_min_pct_within;
#   tested_part and the within_* / refit_* columns (an older table is
#   rebuilt - it is derived from gamm_fitted and the .rds files).
# V3 (Oct 2026): keyed on decomposition.
# =============================================================================

library(dplyr)
library(tibble)
library(mgcv)
library(ggplot2)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/hit_helpers.R")

cfg <- gamm_cfg
decomp <- cfg$metric_decomposition
part <- match.arg(if (is.null(cfg$holdout_part)) "between" else cfg$holdout_part, c("between", "within"))
if (part == "within" && decomp != "between_within")
  stop("holdout_part = 'within' needs metric_decomposition = 'between_within'.")
within_crit <- match.arg(if (is.null(cfg$holdout_within_criterion)) "refit" else cfg$holdout_within_criterion,
                         c("refit", "prediction"))
primary <- if (part == "within") paste0("within_", within_crit) else
  match.arg(cfg$holdout_primary, c("subject", "trial"))
targets <- cfg$holdout_targets
set.seed(cfg$holdout_seed)

data_con    <- connect_data_db(cfg$data_db_path, read_only = TRUE)
results_con <- init_results_db(cfg$results_db_path)

ensure_derived_table(results_con, "holdout_replication",
                     c("decomposition", "tested_part", "within_p", "refit_p"), "
  CREATE TABLE IF NOT EXISTS holdout_replication (
    fitted_at          VARCHAR,
    feature_set        VARCHAR,
    decomposition      VARCHAR,
    tested_part        VARCHAR,
    level              VARCHAR,
    group_value        VARCHAR,
    metric             VARCHAR,
    status             VARCHAR,
    pre_registered     BOOLEAN,
    screen_p           DOUBLE,
    n_holdout_subj     INTEGER,
    n_holdout_trials   INTEGER,
    trial_mse_base     DOUBLE,
    trial_mse_metric   DOUBLE,
    trial_dmse         DOUBLE,
    trial_dmse_lo      DOUBLE,
    trial_dmse_hi      DOUBLE,
    trial_pct          DOUBLE,
    trial_p            DOUBLE,
    subj_mse_base      DOUBLE,
    subj_mse_metric    DOUBLE,
    subj_dmse          DOUBLE,
    subj_dmse_lo       DOUBLE,
    subj_dmse_hi       DOUBLE,
    subj_pct           DOUBLE,
    subj_p             DOUBLE,
    subj_r_base        DOUBLE,
    subj_r_metric      DOUBLE,
    within_mse_base    DOUBLE,
    within_mse_metric  DOUBLE,
    within_dmse        DOUBLE,
    within_dmse_lo     DOUBLE,
    within_dmse_hi     DOUBLE,
    within_pct         DOUBLE,
    within_p           DOUBLE,
    train_direction    DOUBLE,
    refit_p            DOUBLE,
    refit_edf          DOUBLE,
    refit_direction    DOUBLE,
    refit_n_subj       INTEGER,
    refit_dropped      VARCHAR,
    primary_level      VARCHAR,
    primary_p_fdr      DOUBLE,
    replicated         BOOLEAN,
    PRIMARY KEY (feature_set, decomposition, tested_part, level, group_value, metric)
  )")

hits <- load_hits(results_con, cfg$feature_set, p_sig = cfg$hits_p_sig, decomposition = decomp,
                  part = part, targets = targets)
if (nrow(hits) == 0)
  stop(if (is.null(targets)) paste0("No single-metric models with a significant ", part, " part in gamm_fitted.")
       else "None of holdout_targets was found in gamm_fitted (check level / unit / metric names).")
message("Held-out replication: ", nrow(hits), " ", if (!is.null(targets)) "pre-registered " else "",
        "hit(s); tested part = ", part, "; primary criterion = ", primary, "-level prediction")

# -----------------------------------------------------------------------------
# One hit
# -----------------------------------------------------------------------------
evaluate_hit <- function(h) {
  row <- tibble(fitted_at = as.character(Sys.time()), feature_set = cfg$feature_set,
                decomposition = decomp, tested_part = part,
                level = h$level, group_value = h$group_value, metric = h$metrics,
                status = "ok", pre_registered = !is.null(targets), screen_p = h$p_test)
  if (!h$rds_ok) return(mutate(row, status = "rds_missing"))
  mod <- readRDS(h$rds_path)

  ho <- prepare_holdout(data_con, h$level, h$group_value, h$metrics, mod, cfg)
  if (is.null(ho)) return(mutate(row, status = "model_has_no_zscale (re-run the screen)"))
  if (nrow(ho) == 0) return(mutate(row, status = "no_holdout_data"))

  within_var <- if (part == "within") metric_parts(h$metrics, "between_within")[["within"]] else NULL
  base <- tryCatch(refit_baseline(mod, h$metrics, cfg, drop = within_var), error = function(e) NULL)
  if (is.null(base)) return(mutate(row, status = "baseline_fit_failed"))

  ho$pred_metric <- predict_population(mod, ho)
  ho$pred_base   <- predict_population(base, ho)
  ho <- ho %>% filter(!is.na(pred_metric), !is.na(pred_base))
  n_subj <- n_distinct(ho$global_subjid)
  if (n_subj < cfg$holdout_min_subjects)
    return(mutate(row, status = "too_few_holdout_subjects", n_holdout_subj = n_subj))

  ho <- ho %>% mutate(e_base = (pain_rating - pred_base)^2,
                      e_metric = (pain_rating - pred_metric)^2)

  # Trial level: per-subject sums so the bootstrap resamples subjects
  ts <- ho %>% group_by(global_subjid) %>%
    summarise(num = sum(e_metric - e_base), den = n(), .groups = "drop")
  bt <- boot_subjects(ts$num, ts$den, cfg$holdout_n_boot)

  # Subject level: each subject's mean observed vs mean predicted rating
  ss <- ho %>% group_by(global_subjid) %>%
    summarise(obs = mean(pain_rating), pb = mean(pred_base), pm = mean(pred_metric),
              .groups = "drop") %>%
    mutate(d = (obs - pm)^2 - (obs - pb)^2)
  bs <- boot_subjects(ss$d, rep(1, nrow(ss)), cfg$holdout_n_boot)

  mse_tb <- mean(ho$e_base); mse_tm <- mean(ho$e_metric)
  mse_sb <- mean((ss$obs - ss$pb)^2); mse_sm <- mean((ss$obs - ss$pm)^2)
  safe_cor <- function(a, b) if (stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_

  row <- mutate(row,
         n_holdout_subj = n_subj, n_holdout_trials = nrow(ho),
         trial_mse_base = mse_tb, trial_mse_metric = mse_tm,
         trial_dmse = bt[["est"]], trial_dmse_lo = bt[["lo"]], trial_dmse_hi = bt[["hi"]],
         trial_pct = 100 * bt[["est"]] / mse_tb, trial_p = bt[["p"]],
         subj_mse_base = mse_sb, subj_mse_metric = mse_sm,
         subj_dmse = bs[["est"]], subj_dmse_lo = bs[["lo"]], subj_dmse_hi = bs[["hi"]],
         subj_pct = 100 * bs[["est"]] / mse_sb, subj_p = bs[["p"]],
         subj_r_base = safe_cor(ss$obs, ss$pb), subj_r_metric = safe_cor(ss$obs, ss$pm))

  if (part == "within") {
    we <- within_person_errors(ho)
    bw <- boot_subjects(we$num, we$den, cfg$holdout_n_boot)
    mse_wb <- sum(we$sse_b) / sum(we$den)
    rf <- tryCatch(refit_holdout(mod, ho, within_var, cfg), error = function(e) {
      message("    [WARN] held-out refit: ", conditionMessage(e)); NULL })
    row <- mutate(row,
      within_mse_base = mse_wb, within_mse_metric = sum(we$sse_m) / sum(we$den),
      within_dmse = bw[["est"]], within_dmse_lo = bw[["lo"]], within_dmse_hi = bw[["hi"]],
      within_pct = 100 * bw[["est"]] / mse_wb, within_p = bw[["p"]],
      train_direction = term_direction(mod, within_var),
      refit_p = if (is.null(rf)) NA_real_ else rf$p,
      refit_edf = if (is.null(rf)) NA_real_ else rf$edf,
      refit_direction = if (is.null(rf)) NA_real_ else rf$direction,
      refit_n_subj = if (is.null(rf)) NA_integer_ else as.integer(rf$n_subj),
      refit_dropped = if (is.null(rf)) NA_character_ else rf$dropped)
  }
  row
}

res <- vector("list", nrow(hits))
for (i in seq_len(nrow(hits))) {
  h <- hits[i, ]
  message(sprintf("  [%d/%d] %-8s %-12s %s", i, nrow(hits), h$level, h$group_value, h$metrics))
  res[[i]] <- tryCatch(evaluate_hit(h), error = function(e) {
    message("    [WARN] ", conditionMessage(e))
    tibble(fitted_at = as.character(Sys.time()), feature_set = cfg$feature_set,
           decomposition = decomp, tested_part = part,
           level = h$level, group_value = h$group_value, metric = h$metrics,
           status = paste("error:", conditionMessage(e)), pre_registered = !is.null(targets),
           screen_p = h$p_test)
  })
}
res <- bind_rows(res)

# Primary criterion + BH across all evaluated hits
crit <- switch(primary,
  subject = list(p = "subj_p",   d = "subj_dmse",   pct = "subj_pct",   lo = "subj_dmse_lo",
                 hi = "subj_dmse_hi",   base = "subj_mse_base",   floor = cfg$holdout_min_pct),
  trial   = list(p = "trial_p",  d = "trial_dmse",  pct = "trial_pct",  lo = "trial_dmse_lo",
                 hi = "trial_dmse_hi",  base = "trial_mse_base",  floor = cfg$holdout_min_pct),
  within_prediction = list(p = "within_p", d = "within_dmse", pct = "within_pct", lo = "within_dmse_lo",
                 hi = "within_dmse_hi", base = "within_mse_base",
                 floor = if (is.null(cfg$holdout_min_pct_within)) 0 else cfg$holdout_min_pct_within),
  # refit: the p of the within term in the held-out-only refit; the forest
  # plot still shows the subject-centred prediction change
  within_refit = list(p = "refit_p", d = "within_dmse", pct = "within_pct", lo = "within_dmse_lo",
                 hi = "within_dmse_hi", base = "within_mse_base", floor = NA))
for (cn in unlist(crit[c("p", "d", "pct", "lo", "hi", "base")]))
  if (!cn %in% names(res)) res[[cn]] <- NA_real_
res$primary_level <- primary
res$primary_p_fdr <- NA_real_
ok <- !is.na(res[[crit$p]])
res$primary_p_fdr[ok] <- p.adjust(res[[crit$p]][ok], method = "BH")
# Replicates:
#   prediction criteria: an improvement of at least `floor` % that also
#     survives BH;
#   within_refit: the within term is significant (BH) in the held-out-only
#     refit AND points the same way as in training.
res$replicated <- if (primary == "within_refit") {
  ok & res$primary_p_fdr < cfg$holdout_alpha &
    !is.na(res$refit_direction) & !is.na(res$train_direction) &
    sign(res$refit_direction) == sign(res$train_direction) & res$train_direction != 0
} else {
  ok & res[[crit$d]] < 0 & res[[crit$pct]] <= -crit$floor & res$primary_p_fdr < cfg$holdout_alpha
}

# Write (all schema columns, in order)
cols <- dbListFields(results_con, "holdout_replication")
for (cn in setdiff(cols, names(res))) res[[cn]] <- NA
res_out <- as.data.frame(res[, cols])
for (i in seq_len(nrow(res_out))) {
  invisible(dbExecute(results_con,
    "DELETE FROM holdout_replication WHERE feature_set = ? AND decomposition = ? AND tested_part = ?
       AND level = ? AND group_value = ? AND metric = ?",
    params = list(res_out$feature_set[i], res_out$decomposition[i], res_out$tested_part[i],
                  res_out$level[i], res_out$group_value[i], res_out$metric[i])))
}
invisible(dbAppendTable(results_con, "holdout_replication", res_out))

# -----------------------------------------------------------------------------
# Forest plot per level: primary % change in error with 95% CI
# -----------------------------------------------------------------------------
crit_name <- c(subject = "Subject-level (mean rating)", trial = "Trial-level",
               within_prediction = "Within-person (subject-centred trial)",
               within_refit = "Within-person (subject-centred trial)")[[primary]]
rep_rule <- if (primary == "within_refit")
  paste0("Replicated = within term p < ", cfg$holdout_alpha, " (BH) in a refit on held-out subjects only, ",
         "same direction as training.") else
  paste0("Replicated = at least ", crit$floor, "% lower error with BH-adjusted p < ", cfg$holdout_alpha, ".")
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
for (lv in unique(res$level[ok])) {
  fd <- res[ok & res$level == lv, ] %>%
    mutate(pct = .data[[crit$pct]],
           lo = 100 * .data[[crit$lo]] / .data[[crit$base]],
           hi = 100 * .data[[crit$hi]] / .data[[crit$base]],
           lab = sprintf("%s, %s", metric_label(metric), unit_label(group_value)),
           lab = factor(lab, levels = lab[order(pct, decreasing = TRUE)]),
           grp = ifelse(replicated, "Replicated", "Not replicated"))
  n_rep <- sum(fd$replicated)
  p <- ggplot(fd, aes(x = pct, y = lab, colour = grp)) +
    geom_vline(xintercept = 0, colour = INK_3, linewidth = 0.4, linetype = "dashed") +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.6) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = c("Replicated" = cfg$col_main, "Not replicated" = INK_3),
                        name = NULL, drop = FALSE) +
    labs(title = sprintf("Held-out prediction (%s part), %s level: %d of %d %s replicate",
                         part, lv, n_rep, nrow(fd), if (!is.null(targets)) "pre-registered tests" else "hits"),
         subtitle = sprintf("%s prediction error, metric model vs baseline (negative = metric predicts better)",
                            crit_name),
         x = "Change in mean squared error (%)", y = NULL,
         caption = paste0("Models fit on training subjects; predictions for held-out subjects without ",
                          "their random effect.\n95% CI from resampling held-out subjects. ", rep_rule)) +
    theme_gamm() +
    theme(legend.position = "top", legend.justification = "left",
          panel.grid.major.y = element_blank())
  ggsave(file.path(cfg$fig_dir, sprintf("holdout_replication_%s%s%s.png", lv,
                                        if (decomp == "pooled") "_pooled" else "",
                                        if (part == "within") "_within" else "")), p,
         width = cfg$fig_width * 1.3, height = max(3, 1.6 + 0.32 * nrow(fd)),
         dpi = cfg$fig_dpi, bg = "white")
}

disconnect_db(data_con)
close_results_db(results_con)

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------
message("\nStatus: ", paste(sprintf("%s = %d", names(table(res$status)), table(res$status)),
                            collapse = ", "))
message(sprintf("Replicated (%s part, %s-level, BH < %.2f): %d of %d evaluated %s",
                part, primary, cfg$holdout_alpha, sum(res$replicated), sum(ok),
                if (!is.null(targets)) "pre-registered tests" else "hits"))
show <- res[ok, ] %>% arrange(desc(replicated), .data[[crit$p]])
for (i in seq_len(min(30, nrow(show)))) {
  if (part == "within") {
    message(sprintf(paste0("  %-3s %-8s %-12s %-28s within %+7.3f%% [%+.3f, %+.3f]  p %s (FDR %s)\n",
                           "      held-out refit: p %s, edf %.2f, direction %+.2f (training %+.2f), %s held-out subjects"),
                    ifelse(show$replicated[i], "YES", "no"), show$level[i], show$group_value[i],
                    show$metric[i], show$within_pct[i],
                    100 * show$within_dmse_lo[i] / show$within_mse_base[i],
                    100 * show$within_dmse_hi[i] / show$within_mse_base[i],
                    sub("^= ", "", fmt_p(show$within_p[i])), sub("^= ", "", fmt_p(show$primary_p_fdr[i])),
                    sub("^= ", "", fmt_p(show$refit_p[i])), show$refit_edf[i], show$refit_direction[i],
                    show$train_direction[i], show$refit_n_subj[i]))
  } else {
    message(sprintf("  %-3s %-8s %-12s %-28s trial %+6.1f%%  subject %+6.1f%%  (FDR %s)",
                    ifelse(show$replicated[i], "YES", "no"), show$level[i], show$group_value[i],
                    show$metric[i], show$trial_pct[i], show$subj_pct[i],
                    sub("^= ", "", fmt_p(show$primary_p_fdr[i]))))
  }
}
message("\nTable: ", cfg$results_db_path, " (holdout_replication)")