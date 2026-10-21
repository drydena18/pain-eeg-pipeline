# =============================================================================
# compare_levels.R
# =============================================================================
# Channel space vs source space: does source localization help? Run from
# R/gamm/ after source_core.R, channel_core.R and holdout_replicate.R.
# Replaces R/analysis/behavioural_analysis/compare_channel_source_gamm.R
# (CSV outputs and metric names from before v3.1).
#
# For every screened metric (single-metric models, current feature_set and
# decomposition):
#
#   1. Screen summary per level: units tested, units with a hit (p <
#      hits_p_sig), the share of units with a hit, the smallest p and its
#      unit; held-out summary: hits evaluated / replicated in
#      holdout_replication, best held-out change in error.
#      The share, not the count, is compared: 64 channels vs 18 ROIs.
#
#   2. Paired held-out comparison, for metrics with a hit at BOTH levels.
#      The best ROI and the best channel are chosen by SCREEN p (training
#      data only, so the held-out set is not used for selection). Both
#      models predict the same held-out subjects (population level); per
#      subject, squared error of the mean rating, source - channel
#      (negative = source predicts better). 95% CI and two-sided p from
#      resampling subjects. Caveat: the best of 64 channels has more chances
#      to fit well by luck than the best of 18 ROIs (the figure says so).
#
#   3. Curve concordance, same pairs: correlation of the two fitted
#      primary-part curves (z-scale, over the overlap of their 1st-99th
#      percentile ranges) - do the two spaces describe the same shape?
#
# Writes: <out_dir>/compare/level_summary.csv, level_paired_holdout.csv
#         <fig_dir>/level_comparison[_pooled].png
# Reads:  results.duckdb (gamm_fitted, holdout_replication), the .rds files,
#         data.duckdb (held-out subjects, read-only)
#
# V1 (Oct 2026).
# =============================================================================

library(dplyr)
library(tibble)
library(tidyr)
library(mgcv)
library(ggplot2)
library(patchwork)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/hit_helpers.R")

cfg    <- gamm_cfg
cc     <- cfg$compare
decomp <- cfg$metric_decomposition
set.seed(cfg$holdout_seed)
out_dir <- file.path(cfg$out_dir, "compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)

data_con    <- connect_data_db(cfg$data_db_path, read_only = TRUE)
results_con <- init_results_db(cfg$results_db_path)

# =============================================================================
# 1. Screen + held-out summary per metric x level
# =============================================================================
fits <- dbGetQuery(results_con, "
  SELECT level, group_value, metrics, p_value, rds_path
  FROM gamm_fitted
  WHERE feature_set = ? AND decomposition = ? AND n_metrics = 1",
  params = list(cfg$feature_set, decomp)) %>%
  as_tibble() %>%
  mutate(p_num = suppressWarnings(as.numeric(p_value)))
if (nrow(fits) == 0) stop("No single-metric models in gamm_fitted for feature_set ",
                          cfg$feature_set, ", decomposition ", decomp, ".")
if (n_distinct(fits$level) < 2) message("[WARN] Only one level in gamm_fitted (",
                                        paste(unique(fits$level), collapse = ", "),
                                        "): the comparison needs both screens.")

ho <- if (dbExistsTable(results_con, "holdout_replication") &&
          "decomposition" %in% dbListFields(results_con, "holdout_replication")) {
  as_tibble(dbGetQuery(results_con, "
    SELECT level, group_value, metric AS metrics, status, subj_pct, trial_pct, replicated
    FROM holdout_replication WHERE feature_set = ? AND decomposition = ?",
    params = list(cfg$feature_set, decomp)))
} else tibble(level = character(), group_value = character(), metrics = character(),
              status = character(), subj_pct = numeric(), trial_pct = numeric(), replicated = logical())
pct_col <- if (cfg$holdout_primary == "subject") "subj_pct" else "trial_pct"

screen_sum <- fits %>%
  group_by(metrics, level) %>%
  summarise(n_units = n(),
            n_hits = sum(p_num < cfg$hits_p_sig, na.rm = TRUE),
            share_hits = n_hits / n_units,
            min_p = suppressWarnings(min(p_num, na.rm = TRUE)),
            best_unit = group_value[which.min(replace(p_num, is.na(p_num), Inf))],
            .groups = "drop")
ho_sum <- ho %>% filter(status == "ok") %>%
  group_by(metrics, level) %>%
  summarise(n_evaluated = n(), n_replicated = sum(as.logical(replicated), na.rm = TRUE),
            best_holdout_pct = min(.data[[pct_col]], na.rm = TRUE),
            best_holdout_unit = group_value[which.min(.data[[pct_col]])],
            .groups = "drop")
summary_tab <- screen_sum %>% left_join(ho_sum, by = c("metrics", "level")) %>%
  mutate(n_evaluated = coalesce(n_evaluated, 0L), n_replicated = coalesce(n_replicated, 0L)) %>%
  arrange(metrics, level)
readr::write_csv(summary_tab, file.path(out_dir, "level_summary.csv"))
message("Screen / held-out summary: ", file.path(out_dir, "level_summary.csv"))

# =============================================================================
# 2 + 3. Paired held-out comparison and curve concordance
# =============================================================================
best <- fits %>% filter(p_num < cfg$hits_p_sig, file.exists(rds_path)) %>%
  group_by(metrics, level) %>% slice_min(p_num, n = 1, with_ties = FALSE) %>% ungroup()
pair_metrics <- best %>% count(metrics) %>% filter(n == 2) %>% pull(metrics)
message("Metrics with a hit at both levels: ", length(pair_metrics))

subject_errors <- function(mod, level, unit, metric) {
  d <- prepare_holdout(data_con, level, unit, metric, mod, cfg)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$pred <- predict_population(mod, d)
  d %>% filter(!is.na(pred)) %>% group_by(global_subjid) %>%
    summarise(obs = mean(pain_rating), pred = mean(pred), n = n(), .groups = "drop") %>%
    mutate(se = (obs - pred)^2)
}

curve_cor <- function(m_s, m_c, metric) {
  vs <- primary_var(m_s, metric); vc <- primary_var(m_c, metric)
  qs <- stats::quantile(m_s$model[[vs]], c(0.01, 0.99), na.rm = TRUE)
  qc <- stats::quantile(m_c$model[[vc]], c(0.01, 0.99), na.rm = TRUE)
  lo <- max(qs[1], qc[1]); hi <- min(qs[2], qc[2])
  if (!(hi > lo)) return(NA_real_)
  g <- seq(lo, hi, length.out = cc$grid_n)
  fs <- predict_from_lp(m_s, marginal_lpmatrix(m_s, setNames(data.frame(g), vs)))$fit
  fc <- predict_from_lp(m_c, marginal_lpmatrix(m_c, setNames(data.frame(g), vc)))$fit
  if (stats::sd(fs) < 1e-8 || stats::sd(fc) < 1e-8) return(NA_real_)
  stats::cor(fs, fc)
}

paired <- lapply(pair_metrics, function(m) {
  bs <- best %>% filter(metrics == m, level == "source")
  bc <- best %>% filter(metrics == m, level == "channel")
  row <- tibble(metric = m, source_unit = bs$group_value, channel_unit = bc$group_value,
                source_screen_p = bs$p_num, channel_screen_p = bc$p_num)
  m_s <- readRDS(bs$rds_path); m_c <- readRDS(bc$rds_path)
  row$curve_r <- tryCatch(curve_cor(m_s, m_c, m), error = function(e) NA_real_)
  es <- subject_errors(m_s, "source", bs$group_value, m)
  ec <- subject_errors(m_c, "channel", bc$group_value, m)
  if (is.null(es) || is.null(ec)) return(mutate(row, status = "no_holdout_data"))
  j <- inner_join(es, ec, by = "global_subjid", suffix = c("_s", "_c"))
  if (nrow(j) < cfg$holdout_min_subjects) return(mutate(row, status = "too_few_common_subjects",
                                                        n_subj = nrow(j)))
  d <- j$se_s - j$se_c
  B <- cc$n_boot
  bt <- replicate(B, mean(sample(d, replace = TRUE)))
  mutate(row, status = "ok", n_subj = nrow(j),
         mse_source = mean(j$se_s), mse_channel = mean(j$se_c),
         diff = mean(d), diff_lo = unname(quantile(bt, 0.025)), diff_hi = unname(quantile(bt, 0.975)),
         diff_pct = 100 * mean(d) / mean(j$se_c),
         p_two_sided = min(1, 2 * min((1 + sum(bt >= 0)) / (B + 1), (1 + sum(bt <= 0)) / (B + 1))))
})
paired <- bind_rows(paired)
if (nrow(paired) > 0) {
  readr::write_csv(paired, file.path(out_dir, "level_paired_holdout.csv"))
  message("Paired held-out comparison: ", file.path(out_dir, "level_paired_holdout.csv"))
}

# =============================================================================
# Figure
# =============================================================================
lv_cols <- c(source = cfg$col_main, channel = "#eb6834")
fa <- summary_tab %>%
  mutate(lab = metric_label(metrics),
         lab = factor(lab, levels = unique(lab[order(ave(share_hits, metrics, FUN = max))])),
         y = as.numeric(lab) + ifelse(level == "source", 0.15, -0.15))   # the two levels side by side
pA <- ggplot(fa, aes(x = share_hits, y = y, colour = level)) +
  geom_point(size = 2.3, alpha = 0.9) +
  scale_colour_manual(values = lv_cols, name = NULL) +
  scale_y_continuous(breaks = seq_along(levels(fa$lab)), labels = levels(fa$lab)) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(subtitle = sprintf("A  Share of units with a screen hit (p < %.2g)", cfg$hits_p_sig),
       x = "Units with a hit", y = NULL) +
  theme_gamm() + theme(legend.position = "top", legend.justification = "left")

plots <- pA
if (nrow(paired) > 0 && any(paired$status == "ok")) {
  pB <- paired %>% filter(status == "ok") %>%
    mutate(lab = sprintf("%s\n%s vs %s", metric_label(metric), unit_label(source_unit), channel_unit),
           lo = 100 * diff_lo / mse_channel, hi = 100 * diff_hi / mse_channel) %>%
    ggplot(aes(x = diff_pct, y = reorder(lab, -diff_pct))) +
    geom_vline(xintercept = 0, colour = INK_3, linetype = "dashed", linewidth = 0.4) +
    geom_linerange(aes(xmin = lo, xmax = hi), colour = INK_2, linewidth = 0.6) +
    geom_point(colour = INK_2, size = 2.3) +
    labs(subtitle = "B  Held-out error, best ROI vs best channel (by screen p)",
         x = "Source - channel, % of channel error (negative = source better)", y = NULL) +
    theme_gamm()
  plots <- pA | pB
}
n_u <- summary_tab %>% group_by(level) %>% summarise(n = max(n_units), .groups = "drop")
n_of <- function(lv) if (lv %in% n_u$level) n_u$n[n_u$level == lv] else 0L
plots <- plots + plot_annotation(
  title = "Channel space vs source space",
  caption = paste0("Feature set ", cfg$feature_set, ", ", decomp, " decomposition. B: 95% CI from resampling held-out subjects; ",
                   "the best of ", n_of("channel"), " channels vs the best of ", n_of("source"),
                   " ROIs favours the level with more units."),
  theme = theme_gamm())
ggsave(file.path(cfg$fig_dir, sprintf("level_comparison%s.png", if (decomp == "pooled") "_pooled" else "")),
       plots, width = cfg$fig_width * 2, height = max(4, 1.5 + 0.28 * n_distinct(summary_tab$metrics)),
       dpi = cfg$fig_dpi, bg = "white")

# =============================================================================
# Console
# =============================================================================
message("\nShare of units with a hit, mean over metrics: ",
        paste(sprintf("%s %.1f%%", c("source", "channel"),
                      100 * c(mean(summary_tab$share_hits[summary_tab$level == "source"]),
                              mean(summary_tab$share_hits[summary_tab$level == "channel"]))),
              collapse = " | "))
if (nrow(paired) > 0) for (i in seq_len(nrow(paired))) {
  r <- paired[i, ]
  if (!identical(r$status, "ok")) { message(sprintf("  %-28s %s", r$metric, r$status)); next }
  message(sprintf("  %-28s %-10s vs %-5s  source - channel %+6.1f%% [%+.3f, %+.3f]  p %s  curve r = %.2f",
                  r$metric, r$source_unit, r$channel_unit, r$diff_pct, r$diff_lo, r$diff_hi,
                  fmt_p(r$p_two_sided), r$curve_r))
}
disconnect_db(data_con)
close_results_db(results_con)