# =============================================================================
# classical_tests.R
# =============================================================================
# Classical inferential tests that complement the GAMM pipeline, on the v3.1
# tables in data.duckdb, at channel level (electrode clusters / electrodes)
# or source level (ROIs). Run from R/gamm/.
# Replaces R/analysis/behavioural_analysis/run_classical_tests.R (CSV inputs
# and metric names from before v3.1: BI_pre, CoG_pre, delta_ERD, ...).
#
#   Test 1  Paired t-test   Does alpha power drop after the laser (ERD)?
#                           Subject-level pre vs post power per unit x band
#                           (slow 8-10, fast 10-12, total 8-12 Hz), from the
#                           subject-level table (trial-mean powers).
#   Test 2  rmANOVA         Do subject-level metrics differ across units?
#                           (classical$rmanova_metrics; default
#                           pre_sf_balance, pre_paf_cog_hz, erd_asym;
#                           Greenhouse-Geisser via afex if installed) +
#                           BH-corrected pairwise paired t-tests. With the
#                           two S1 clusters this is the contralateral vs
#                           ipsilateral comparison.
#   Test 3  Rayleigh test   Is pre-stimulus (and post) slow-alpha phase
#                           non-uniform across a subject's trials?
#                           Per subject x unit, from the trial-level table.
#
# Level and units (V2):
#   classical$level  "channel" or "source"; NULL = "channel" if it is in
#                    gamm_cfg$analysis_levels, else "source". Must be in
#                    analysis_levels.
#   channel level    units = those selected by channel_units / cluster_names
#                    / channel_include, narrowed to classical$units if set.
#   source level     units = classical$rois.
#
# Dropped from the old script: channel-vs-source per-subject R² (now
# compare_levels.R, on held-out prediction) and the ANOVA of GAMM deviance
# across ROIs (the old per-ROI fitted-value CSVs no longer exist; per-ROI
# screen and held-out results are in compare_levels.R's summary).
#
# Subjects: classical$split ("all" by default: these tests describe the
# signal, they are not the predictive claim; use "train" to keep the
# held-out subjects untouched).
#
# Writes (<classical$out_dir>/<level>/): test1_erd_ttest.csv,
#         test2_rmanova_<metric>.csv, test2_posthoc_<metric>.csv,
#         test3_rayleigh_<phase>.csv, test3_rayleigh_summary.csv,
#         classical_tests_report.txt, figures. Unit column: `unit`.
#
# V2 (Nov 2026): channel level (clusters) as well as source; outputs per
# level subfolder; the unit column is `unit` (was `roi`).
# V1 (Oct 2026).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/hit_helpers.R")

cc <- gamm_cfg$classical

# -- Level and units ----
# Pure function (no database) so the selection can be tested.
select_classical_units <- function(available, level, cfg) {
  cc <- cfg$classical
  if (level == "source") return(intersect(cc$rois, available))
  base <- select_units(apply_channel_units(available, cfg$channel_units %||% "all",
                                           cfg$cluster_names %||% character(0)),
                       cfg$channel_include)
  if (is.null(cc$units)) return(base)
  missing <- setdiff(cc$units, base)
  if (length(missing) > 0)
    message("  [WARN] classical$units not available under channel_units / channel_include: ",
            paste(missing, collapse = ", "))
  intersect(cc$units, base)
}

classical_level <- function(cfg) {
  lv <- get_analysis_levels(cfg)
  level <- cfg$classical$level %||% (if ("channel" %in% lv) "channel" else "source")
  if (!level %in% lv)
    stop("classical$level = \"", level, "\" is not in gamm_cfg$analysis_levels (",
         paste(lv, collapse = ", "), ").", call. = FALSE)
  level
}

level    <- classical_level(gamm_cfg)
unit_col <- level_spec(level)$unit_col
unit_lbl <- if (level == "source") "ROI" else "unit"
out_dir  <- file.path(cc$out_dir, level)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
afex_ok <- requireNamespace("afex", quietly = TRUE)

data_con <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
report <- c("CLASSICAL TESTS", strrep("=", 60),
            paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
            paste("Feature set:", gamm_cfg$feature_set, "| level:", level, "| subjects:", cc$split), "")

sig_label <- function(p) dplyr::case_when(is.na(p) ~ "n/a", p < 0.001 ~ "***", p < 0.01 ~ "**",
                                          p < 0.05 ~ "*", TRUE ~ "ns")

paired_t <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(tibble(n = length(x), t = NA_real_, df = NA_real_, p_value = NA_real_,
                                   mean_diff = NA_real_, cohens_dz = NA_real_))
  tt <- t.test(x, y, paired = TRUE)
  tibble(n = length(x), t = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value,
         mean_diff = unname(tt$estimate), cohens_dz = mean(x - y) / sd(x - y))
}

units <- select_classical_units(list_available_units(data_con, level, split = cc$split), level, gamm_cfg)
if (length(units) == 0) stop("No ", unit_lbl, "s to test at ", level, " level (check classical$",
                             if (level == "source") "rois" else "units / channel_units", ").")
message("Level: ", level, " | ", unit_lbl, "s: ", paste(units, collapse = ", "))
ylab_units <- unit_label(units)

subj <- load_subject_metrics_wide(data_con, level, unit = units, split = cc$split)
names(subj)[names(subj) == unit_col] <- "unit"

# =============================================================================
# TEST 1: ERD (pre vs post power)
# =============================================================================
message("\nTEST 1: pre vs post alpha power")
bands <- c(slow = "pow_slow_alpha", fast = "pow_fast_alpha", total = "pow_alpha_total")
t1 <- bind_rows(lapply(units, function(u) bind_rows(lapply(names(bands), function(b) {
  pre <- paste0("pre_", bands[[b]]); post <- paste0("post_", bands[[b]])
  d <- subj %>% filter(unit == u)
  if (!all(c(pre, post) %in% names(d))) return(NULL)
  paired_t(d[[post]], d[[pre]]) %>% mutate(unit = u, band = b, .before = 1)
}))))
if (nrow(t1) > 0) {
  t1 <- t1 %>% mutate(p_adj = p.adjust(p_value, "BH"), sig = sig_label(p_adj))
  readr::write_csv(t1, file.path(out_dir, "test1_erd_ttest.csv"))
  p1 <- ggplot(t1, aes(x = factor(band, c("slow", "fast", "total")),
                       y = factor(unit, rev(units), rev(ylab_units)), fill = cohens_dz)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sig), size = 3.2) +
    scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#eb6834", midpoint = 0, name = "Cohen's dz") +
    labs(title = "Post - pre alpha power (subject level)",
         subtitle = "Negative = power drops after the laser (ERD). Stars: BH-adjusted paired t-test.",
         x = "Band", y = NULL) +
    theme_gamm() + theme(legend.position = "right")
  ggsave(file.path(out_dir, "test1_erd_plot.png"), p1, width = 6.5, height = 1.5 + 0.35 * length(units),
         dpi = gamm_cfg$fig_dpi, bg = "white")
  report <- c(report, "TEST 1  Post vs pre power (paired t, BH)", strrep("-", 40),
              sprintf("  %-12s %-6s n=%-4d t=%7.2f  dz=%6.2f  p(BH)=%.4g %s",
                      t1$unit, t1$band, t1$n, t1$t, t1$cohens_dz, t1$p_adj, t1$sig), "")
} else message("  [SKIP] pre_/post_ power columns not in the subject-level table.")

# =============================================================================
# TEST 2: rmANOVA across units
# =============================================================================
message("\nTEST 2: rmANOVA across ", unit_lbl, "s")
report <- c(report, paste0("TEST 2  Metric differences across ", unit_lbl, "s (rmANOVA)"), strrep("-", 40))
t2_plot <- list()
for (m in cc$rmanova_metrics) {
  if (!m %in% names(subj)) { message("  [SKIP] ", m, " not in the subject-level table."); next }
  d <- subj %>% select(global_subjid, unit, value = all_of(m)) %>% filter(!is.na(value)) %>%
    group_by(global_subjid) %>% filter(n_distinct(unit) == length(units)) %>% ungroup() %>%
    mutate(global_subjid = factor(global_subjid), unit = factor(unit, units))
  if (n_distinct(d$global_subjid) < 3 || length(units) < 2) {
    message("  [SKIP] ", m, ": too few complete subjects / ", unit_lbl, "s."); next
  }
  tab <- tryCatch({
    if (afex_ok) {
      fit <- afex::aov_ez(id = "global_subjid", dv = "value", data = as.data.frame(d), within = "unit",
                          anova_table = list(correction = "GG"))
      a <- as.data.frame(fit$anova_table)
      tibble(metric = m, n_subjects = n_distinct(d$global_subjid), num_df = a$`num Df`, den_df = a$`den Df`,
             F_stat = a$F, p_value = a$`Pr(>F)`, correction = "Greenhouse-Geisser")
    } else {
      s <- summary(aov(value ~ unit + Error(global_subjid / unit), data = d))$`Error: global_subjid:unit`[[1]]
      tibble(metric = m, n_subjects = n_distinct(d$global_subjid), num_df = s["unit", "Df"],
             den_df = s["Residuals", "Df"], F_stat = s["unit", "F value"], p_value = s["unit", "Pr(>F)"],
             correction = "none (install afex for Greenhouse-Geisser)")
    }
  }, error = function(e) { message("  [WARN] ", m, ": ", conditionMessage(e)); NULL })
  if (is.null(tab)) next
  w <- d %>% pivot_wider(names_from = unit, values_from = value)
  ph <- bind_rows(lapply(combn(units, 2, simplify = FALSE), function(pr)
    paired_t(w[[pr[1]]], w[[pr[2]]]) %>% mutate(unit_1 = pr[1], unit_2 = pr[2], .before = 1))) %>%
    mutate(p_adj = p.adjust(p_value, "BH"), sig = sig_label(p_adj))
  readr::write_csv(tab, file.path(out_dir, paste0("test2_rmanova_", m, ".csv")))
  readr::write_csv(ph,  file.path(out_dir, paste0("test2_posthoc_", m, ".csv")))
  report <- c(report, sprintf("  %-16s n=%d  F(%.1f, %.1f) = %.2f  p = %.4g  [%s]  %d of %d pairs differ (BH)",
                              m, tab$n_subjects, tab$num_df, tab$den_df, tab$F_stat, tab$p_value,
                              tab$correction, sum(ph$p_adj < cc$alpha, na.rm = TRUE), nrow(ph)))
  t2_plot[[m]] <- d %>% mutate(metric = metric_label(m))
}
report <- c(report, "")
if (length(t2_plot) > 0) {
  p2 <- ggplot(bind_rows(t2_plot), aes(x = value, y = factor(unit, rev(units), rev(ylab_units)))) +
    geom_boxplot(fill = gamm_cfg$col_main, alpha = 0.25, colour = INK_2, outlier.size = 0.6, width = 0.6) +
    facet_wrap(~ metric, scales = "free_x") +
    labs(title = paste0("Subject-level metrics across ", unit_lbl, "s"), x = NULL, y = NULL) +
    theme_gamm()
  ggsave(file.path(out_dir, "test2_rmanova_plot.png"), p2, width = 3.5 * length(t2_plot) + 1.5,
         height = 1.5 + 0.35 * length(units), dpi = gamm_cfg$fig_dpi, bg = "white")
}

# =============================================================================
# TEST 3: Rayleigh test of phase uniformity
# -----------------------------------------------------------------------------
# z = n R^2, R = mean resultant length; p from the Zar (1999) approximation
#   p = exp(sqrt(1 + 4n + 4(n^2 - (nR)^2)) - (1 + 2n)), accurate for n >= 10.
# =============================================================================
message("\nTEST 3: Rayleigh test")
rayleigh <- function(th) {
  th <- th[!is.na(th)]; n <- length(th)
  if (n < 10) return(c(n = n, R = NA, z = NA, p = NA, mean_phase = NA))
  C <- mean(cos(th)); S <- mean(sin(th)); R <- sqrt(C^2 + S^2)
  p <- exp(sqrt(1 + 4 * n + 4 * (n^2 - (n * R)^2)) - (1 + 2 * n))
  c(n = n, R = R, z = n * R^2, p = min(1, max(0, p)), mean_phase = atan2(S, C))
}
report <- c(report, paste0("TEST 3  Phase non-uniformity (Rayleigh, per subject x ", unit_lbl,
                           ", BH within phase column)"), strrep("-", 40))
t3_sum <- list()
for (ph_col in cc$phase_metrics) {
  r3 <- bind_rows(lapply(units, function(u) {
    d <- load_metrics_wide(data_con, level, u, metrics = ph_col, split = cc$split)
    if (!ph_col %in% names(d)) return(NULL)
    d %>% group_by(global_subjid) %>%
      summarise(st = list(rayleigh(.data[[ph_col]])), .groups = "drop") %>%
      mutate(unit = u, n_trials = sapply(st, `[[`, "n"), R = sapply(st, `[[`, "R"),
             rayleigh_z = sapply(st, `[[`, "z"), p_value = sapply(st, `[[`, "p"),
             mean_phase = sapply(st, `[[`, "mean_phase")) %>%
      select(-st)
  }))
  if (nrow(r3) == 0) { message("  [SKIP] ", ph_col, " not in the trial-level table."); next }
  r3 <- r3 %>% mutate(p_adj = p.adjust(p_value, "BH"), sig = sig_label(p_adj))
  readr::write_csv(r3, file.path(out_dir, paste0("test3_rayleigh_", ph_col, ".csv")))
  s3 <- r3 %>% group_by(unit) %>%
    summarise(n_subjects = sum(!is.na(p_value)), n_sig = sum(p_adj < cc$alpha, na.rm = TRUE),
              pct_sig = 100 * n_sig / n_subjects, median_R = median(R, na.rm = TRUE), .groups = "drop") %>%
    mutate(phase = ph_col, .before = 1)
  t3_sum[[ph_col]] <- s3
  report <- c(report, sprintf("  %-16s %-12s %3d of %3d subjects non-uniform (%.0f%%)  median R = %.3f",
                              s3$phase, s3$unit, s3$n_sig, s3$n_subjects, s3$pct_sig, s3$median_R))
}
if (length(t3_sum) > 0) {
  s3 <- bind_rows(t3_sum)
  readr::write_csv(s3, file.path(out_dir, "test3_rayleigh_summary.csv"))
  p3 <- ggplot(s3, aes(x = pct_sig, y = factor(unit, rev(units), rev(ylab_units)), colour = phase)) +
    geom_vline(xintercept = 100 * cc$alpha, colour = INK_3, linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = setNames(c(gamm_cfg$col_main, "#eb6834")[seq_along(cc$phase_metrics)],
                                          cc$phase_metrics), name = NULL) +
    labs(title = "Slow-alpha phase: subjects with non-uniform phase across trials",
         subtitle = sprintf("Rayleigh test per subject, BH-adjusted; dashed line = %.0f%% expected by chance",
                            100 * cc$alpha),
         x = "Subjects with p(BH) < alpha (%)", y = NULL) +
    theme_gamm() + theme(legend.position = "top", legend.justification = "left")
  ggsave(file.path(out_dir, "test3_rayleigh_plot.png"), p3, width = 6.5, height = 2 + 0.35 * length(units),
         dpi = gamm_cfg$fig_dpi, bg = "white")
}

writeLines(report, file.path(out_dir, "classical_tests_report.txt"))
disconnect_db(data_con)
message("\n", paste(report, collapse = "\n"))
message("\nOutputs: ", out_dir)