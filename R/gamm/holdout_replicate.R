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
# For every significant single-metric model (p < hits_p_sig, current
# feature_set and decomposition, source and channel):
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
#                  (does the metric track trial-to-trial pain?)
#        subject : squared error of each subject's MEAN rating
#                  (does the metric predict who is more pain-sensitive,
#                   beyond laser power, age, sex, cap size, experiment?)
#      Change in MSE = metric - baseline (negative = metric predicts better),
#      with a 95% CI from resampling held-out SUBJECTS (2000 draws) and a
#      one-sided bootstrap p for "no improvement".
#   5. Primary criterion (gamm_cfg$holdout_primary, default "subject"):
#      a hit REPLICATES if its primary error is at least holdout_min_pct %
#      lower than baseline AND its Benjamini-Hochberg-adjusted p across all
#      hits is < holdout_alpha.
#
# Writes: results.duckdb  holdout_replication (one row per hit)
#         <fig_dir>/holdout_replication_<level>[_pooled].png  (forest plot)
# Reads:  data.duckdb (holdout subjects, read-only), gamm_fitted, the .rds
#
# V3 (Oct 2026): keyed on decomposition (column added; an older table is
# rebuilt - it is derived from gamm_fitted and the .rds files).
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
primary <- match.arg(cfg$holdout_primary, c("subject", "trial"))
set.seed(cfg$holdout_seed)

data_con    <- connect_data_db(cfg$data_db_path, read_only = TRUE)
results_con <- init_results_db(cfg$results_db_path)

ensure_derived_table(results_con, "holdout_replication", c("decomposition"), "
  CREATE TABLE IF NOT EXISTS holdout_replication (
    fitted_at          VARCHAR,
    feature_set        VARCHAR,
    decomposition      VARCHAR,
    level              VARCHAR,
    group_value        VARCHAR,
    metric             VARCHAR,
    status             VARCHAR,
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
    primary_level      VARCHAR,
    primary_p_fdr      DOUBLE,
    replicated         BOOLEAN,
    PRIMARY KEY (feature_set, decomposition, level, group_value, metric)
  )")

hits <- load_hits(results_con, cfg$feature_set, p_sig = cfg$hits_p_sig, decomposition = decomp)
if (nrow(hits) == 0) stop("No significant single-metric models in gamm_fitted.")
message("Held-out replication: ", nrow(hits), " hit(s); primary criterion = ", primary,
        "-level prediction")

# -----------------------------------------------------------------------------
# One hit
# -----------------------------------------------------------------------------
evaluate_hit <- function(h) {
  row <- tibble(fitted_at = as.character(Sys.time()), feature_set = cfg$feature_set,
                decomposition = decomp,
                level = h$level, group_value = h$group_value, metric = h$metrics,
                status = "ok", screen_p = h$p_num)
  if (!h$rds_ok) return(mutate(row, status = "rds_missing"))
  mod <- readRDS(h$rds_path)

  ho <- prepare_holdout(data_con, h$level, h$group_value, h$metrics, mod, cfg)
  if (is.null(ho)) return(mutate(row, status = "model_has_no_zscale (re-run the screen)"))
  if (nrow(ho) == 0) return(mutate(row, status = "no_holdout_data"))

  base <- tryCatch(refit_baseline(mod, h$metrics, cfg), error = function(e) NULL)
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

  mutate(row,
         n_holdout_subj = n_subj, n_holdout_trials = nrow(ho),
         trial_mse_base = mse_tb, trial_mse_metric = mse_tm,
         trial_dmse = bt[["est"]], trial_dmse_lo = bt[["lo"]], trial_dmse_hi = bt[["hi"]],
         trial_pct = 100 * bt[["est"]] / mse_tb, trial_p = bt[["p"]],
         subj_mse_base = mse_sb, subj_mse_metric = mse_sm,
         subj_dmse = bs[["est"]], subj_dmse_lo = bs[["lo"]], subj_dmse_hi = bs[["hi"]],
         subj_pct = 100 * bs[["est"]] / mse_sb, subj_p = bs[["p"]],
         subj_r_base = safe_cor(ss$obs, ss$pb), subj_r_metric = safe_cor(ss$obs, ss$pm))
}

res <- vector("list", nrow(hits))
for (i in seq_len(nrow(hits))) {
  h <- hits[i, ]
  message(sprintf("  [%d/%d] %-8s %-12s %s", i, nrow(hits), h$level, h$group_value, h$metrics))
  res[[i]] <- tryCatch(evaluate_hit(h), error = function(e) {
    message("    [WARN] ", conditionMessage(e))
    tibble(fitted_at = as.character(Sys.time()), feature_set = cfg$feature_set,
           decomposition = decomp,
           level = h$level, group_value = h$group_value, metric = h$metrics,
           status = paste("error:", conditionMessage(e)), screen_p = h$p_num)
  })
}
res <- bind_rows(res)

# Primary criterion + BH across all evaluated hits
pcol <- if (primary == "subject") "subj_p" else "trial_p"
dcol <- if (primary == "subject") "subj_dmse" else "trial_dmse"
for (cn in c("trial_p", "subj_p", "trial_dmse", "subj_dmse")) if (!cn %in% names(res)) res[[cn]] <- NA_real_
res$primary_level <- primary
res$primary_p_fdr <- NA_real_
ok <- !is.na(res[[pcol]])
res$primary_p_fdr[ok] <- p.adjust(res[[pcol]][ok], method = "BH")
pcol_pct <- if (primary == "subject") "subj_pct" else "trial_pct"
if (!pcol_pct %in% names(res)) res[[pcol_pct]] <- NA_real_
# Replicates = a meaningful improvement (at least holdout_min_pct % lower
# error) that also survives BH. Without the size floor, a negligible but
# consistent improvement (e.g. -0.01 %) can be "significant".
res$replicated <- ok & res[[dcol]] < 0 &
  res[[pcol_pct]] <= -cfg$holdout_min_pct &
  res$primary_p_fdr < cfg$holdout_alpha

# Write (all schema columns, in order)
cols <- dbListFields(results_con, "holdout_replication")
for (cn in setdiff(cols, names(res))) res[[cn]] <- NA
res_out <- as.data.frame(res[, cols])
for (i in seq_len(nrow(res_out))) {
  invisible(dbExecute(results_con,
    "DELETE FROM holdout_replication WHERE feature_set = ? AND decomposition = ? AND level = ?
       AND group_value = ? AND metric = ?",
    params = list(res_out$feature_set[i], res_out$decomposition[i], res_out$level[i],
                  res_out$group_value[i], res_out$metric[i])))
}
invisible(dbAppendTable(results_con, "holdout_replication", res_out))

# -----------------------------------------------------------------------------
# Forest plot per level: primary % change in MSE with 95% CI
# -----------------------------------------------------------------------------
pct  <- if (primary == "subject") "subj_pct" else "trial_pct"
lo_c <- if (primary == "subject") "subj_dmse_lo" else "trial_dmse_lo"
hi_c <- if (primary == "subject") "subj_dmse_hi" else "trial_dmse_hi"
base_c <- if (primary == "subject") "subj_mse_base" else "trial_mse_base"

dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
for (lv in unique(res$level[ok])) {
  fd <- res[ok & res$level == lv, ] %>%
    mutate(pct = .data[[pct]],
           lo = 100 * .data[[lo_c]] / .data[[base_c]],
           hi = 100 * .data[[hi_c]] / .data[[base_c]],
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
    labs(title = sprintf("Held-out prediction, %s level: %d of %d hits replicate",
                         lv, n_rep, nrow(fd)),
         subtitle = sprintf("%s-level prediction error, metric model vs baseline (negative = metric predicts better)",
                            tools::toTitleCase(primary)),
         x = "Change in mean squared error (%)", y = NULL,
         caption = paste0("Models fit on training subjects; predictions for held-out subjects without ",
                          "their random effect.\n95% CI from resampling held-out subjects. ",
                          "Replicated = at least ", cfg$holdout_min_pct, "% lower error with BH-adjusted p < ",
                          cfg$holdout_alpha, ".")) +
    theme_gamm() +
    theme(legend.position = "top", legend.justification = "left",
          panel.grid.major.y = element_blank())
  ggsave(file.path(cfg$fig_dir, sprintf("holdout_replication_%s%s.png", lv,
                                        if (decomp == "pooled") "_pooled" else "")), p,
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
message(sprintf("Replicated (%s-level, BH < %.2f): %d of %d evaluated hits",
                primary, cfg$holdout_alpha, sum(res$replicated), sum(ok)))
show <- res[ok, ] %>% arrange(desc(replicated), .data[[pcol]])
for (i in seq_len(min(30, nrow(show)))) {
  message(sprintf("  %-3s %-8s %-12s %-28s trial %+6.1f%%  subject %+6.1f%%  (FDR %s)",
                  ifelse(show$replicated[i], "YES", "no"), show$level[i], show$group_value[i],
                  show$metric[i], show$trial_pct[i], show$subj_pct[i],
                  sub("^= ", "", fmt_p(show$primary_p_fdr[i]))))
}
message("\nTable: ", cfg$results_db_path, " (holdout_replication)")