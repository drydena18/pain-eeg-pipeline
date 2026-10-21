# =============================================================================
# plot_gamm_hits.R
# =============================================================================
# Figures for every significant single-metric GAMM (p < gamm_cfg$hits_p_sig,
# current feature_set and decomposition, source and channel). Run from
# R/gamm/ after source_core.R / channel_core.R, and after sex_diff_core.R if
# the sex figures are wanted.
#
# Per hit:
#   <fig_dir>/<level>/<unit>/[pooled/]<metric>_effect.png
#     Predicted pain rating across the metric, in the metric's original
#     units, with a 95% confidence band and a rug of the observed data.
#     With the between/within split: A = subject-mean curve (the tested
#     part, coloured), B = within-person curve (grey), same y-scale.
#   <fig_dir>/<level>/<unit>/[pooled/]<metric>_sex.png  (if gamm_sex_diff has it)
#     A: predicted curve for each sex.  B: the difference curve
#     (other - reference) with its 95% band; x-ranges where the band
#     excludes 0 are shaded.
#
# What "predicted" means here: covariates (laser power, trial index, age)
# at their mean; averaged over sex and cap size in proportion to the sample
# (sex held fixed in the sex figure); subject and experiment random effects
# excluded. So a curve is the expected rating for a typical subject, and the
# y-axis is in rating units, unlike mgcv's plot.gam(), which shows the
# centred partial effect.
#
# Reads only: results.duckdb (gamm_fitted, gamm_sex_diff) and the .rds files.
#
# V3 (Oct 2026): hits and sex rows filtered on decomposition; pooled figures
# go to a pooled/ sub-folder.
# =============================================================================

library(dplyr)
library(tibble)
library(mgcv)
library(ggplot2)
library(patchwork)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/hit_helpers.R")

cfg <- gamm_cfg
decomp <- cfg$metric_decomposition
# Colours, sex labels and theme_gamm() live in helpers/hit_helpers.R.

# Up to n observed values (trials) for the rug, inside the plotted range.
rug_sample <- function(mod, m, rng, n = 1500) {
  z <- mod$model[[m]]
  z <- z[z >= rng[1] & z <= rng[2]]
  if (length(z) > n) z <- sample(z, n)
  tibble(x = to_original_units(mod, m, z)$x)
}

x_axis_title <- function(mod, m) {
  if (to_original_units(mod, m, 0)$raw) metric_label(m)
  else paste0(metric_label(m), " (z-score)")
}

# -----------------------------------------------------------------------------
# Main effect figure
# -----------------------------------------------------------------------------
PART_TITLES <- c(between = "Between people (subject mean)",
                 within  = "Within person (trial deviation from own mean)",
                 pooled  = "")

# One curve. Rug: one tick per subject for a subject-mean term, otherwise a
# sample of trials.
effect_curve <- function(mod, var) {
  g <- metric_grid(mod, var)
  nd <- setNames(data.frame(g), var)
  pr <- predict_from_lp(mod, marginal_lpmatrix(mod, nd)) %>%
    mutate(x = to_original_units(mod, var, g)$x)
  rug <- if (grepl("_between_z$", var)) {
    z <- unique(mod$model[, c("global_subjid", var)])[[var]]
    tibble(x = to_original_units(mod, var, z[z >= min(g) & z <= max(g)])$x)
  } else rug_sample(mod, var, range(g))
  list(pr = pr, rug = rug)
}

effect_panel <- function(cv, mod, var, col, subtitle, ylim) {
  ggplot(cv$pr, aes(x, fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = col, alpha = 0.18) +
    geom_line(colour = col, linewidth = 1) +
    geom_rug(data = cv$rug, aes(x = x), inherit.aes = FALSE, sides = "b",
             colour = INK_2, alpha = if (grepl("_between_z$", var)) 0.5 else 0.12,
             length = grid::unit(0.025, "npc")) +
    coord_cartesian(ylim = ylim) +
    labs(subtitle = subtitle, x = x_axis_title(mod, var), y = "Predicted pain rating") +
    theme_gamm()
}

stat_text <- function(edf, f, p) sprintf("edf = %.1f | F = %.1f | p %s", edf, f, fmt_p(p))

plot_effect <- function(hit, mod) {
  m <- hit$metrics
  vp <- primary_var(mod, m)
  vo <- other_var(mod, m)
  split <- !is.na(vo)

  n_line <- sprintf("%s subjects, %s trials%s",
                    format(hit$n_subjects, big.mark = ","), format(hit$n_obs, big.mark = ","),
                    if (!is.na(hit$delta_aic)) sprintf(" | AIC change vs baseline = %.1f", hit$delta_aic) else "")
  caption <- paste0("Covariates at their mean; averaged over sex and cap size; subject and ",
                    "experiment effects excluded.\nBand: 95% CI. Ticks: ",
                    if (split) "subjects (A), sample of trials (B). " else "observed trials (sample). ",
                    tools::toTitleCase(hit$level), " level, training subjects only.",
                    if (split) paste0("\nIn-sample AIC barely changes for a between-person effect (subject ",
                                      "intercepts already fit each person's mean): judge it by held-out prediction.") else "")
  title <- sprintf("%s, %s", metric_label(m), unit_label(hit$group_value))

  cp <- effect_curve(mod, vp)
  if (!split) {
    return(effect_panel(cp, mod, vp, cfg$col_main, NULL, range(c(cp$pr$lo, cp$pr$hi))) +
             labs(title = title,
                  subtitle = paste(stat_text(hit$edf_num, hit$f_num, hit$p_num), "|", n_line),
                  caption = caption))
  }

  co <- effect_curve(mod, vo)
  ylim <- range(c(cp$pr$lo, cp$pr$hi, co$pr$lo, co$pr$hi))   # same y-scale in both panels
  part_p <- mod$gamm_meta$primary_part
  part_o <- setdiff(c("between", "within"), part_p)
  pA <- effect_panel(cp, mod, vp, cfg$col_main,
                     sprintf("A  %s\n%s  (tested)", PART_TITLES[[part_p]],
                             stat_text(hit$edf_num, hit$f_num, hit$p_num)), ylim)
  pB <- effect_panel(co, mod, vo, INK_2,
                     sprintf("B  %s\n%s", PART_TITLES[[part_o]],
                             stat_text(hit$edf_other_num, hit$f_other_num, hit$p_other_num)), ylim)
  (pA | pB) +
    plot_annotation(title = title, subtitle = n_line, caption = caption, theme = theme_gamm()) &
    theme(plot.subtitle = element_text(colour = INK_2))
}

# -----------------------------------------------------------------------------
# Sex figure: curves by sex + difference curve
# -----------------------------------------------------------------------------
excl_zero_runs <- function(x, lo, hi) {
  flag <- lo > 0 | hi < 0
  if (!any(flag)) return(tibble(xmin = numeric(), xmax = numeric()))
  r <- rle(flag)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  tibble(xmin = x[starts[r$values]], xmax = x[ends[r$values]])
}

plot_sex <- function(hit, sx, fit) {
  m <- primary_var(fit, hit$metrics)
  ref <- sx$sex_reference
  oth <- sx$sex_other
  g <- metric_grid(fit, m)
  nd <- setNames(data.frame(g), m)
  xr <- to_original_units(fit, m, g)$x

  X_ref <- marginal_lpmatrix(fit, nd, fix = list(sex_o = ref))
  X_oth <- marginal_lpmatrix(fit, nd, fix = list(sex_o = oth))
  n_by <- c(setNames(sx$n_subj_ref, ref), setNames(sx$n_subj_other, oth))
  sex_lab <- function(s) sprintf("%s (n = %d)", sex_name(s), as.integer(n_by[as.character(s)]))

  curves <- bind_rows(
    predict_from_lp(fit, X_ref) %>% mutate(x = xr, sex = ref),
    predict_from_lp(fit, X_oth) %>% mutate(x = xr, sex = oth)
  ) %>% mutate(sex = factor(sex, levels = c(ref, oth)))
  cols <- cfg$col_sex[c(ref, oth)]
  if (any(is.na(cols))) cols <- setNames(c("#eb6834", "#2a78d6"), c(ref, oth))
  ends <- curves %>% group_by(sex) %>% slice_max(x, n = 1) %>% ungroup() %>%
    mutate(lab = vapply(as.character(sex), sex_lab, character(1)))

  dx <- diff(range(xr))
  pA <- ggplot(curves, aes(x, fit, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_text(data = ends, aes(label = lab), hjust = 0, nudge_x = 0.02 * dx,
              size = 3.3, colour = INK_2) +
    scale_colour_manual(values = cols, labels = sex_lab, name = NULL) +
    scale_fill_manual(values = cols, labels = sex_lab, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.28))) +
    labs(subtitle = "A  Predicted pain by sex", x = x_axis_title(fit, m),
         y = "Predicted pain rating") +
    theme_gamm() +
    theme(legend.position = "top", legend.justification = "left",
          plot.subtitle = element_text(colour = INK, face = "bold"))

  dpr <- predict_from_lp(fit, X_oth - X_ref) %>% mutate(x = xr)
  runs <- excl_zero_runs(dpr$x, dpr$lo, dpr$hi)
  diff_name <- sprintf("%s - %s", sex_name(oth), sex_name(ref))

  pB <- ggplot(dpr, aes(x, fit)) +
    { if (nrow(runs) > 0) geom_rect(data = runs, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                    inherit.aes = FALSE, fill = INK_3, alpha = 0.12) } +
    geom_hline(yintercept = 0, colour = INK_3, linewidth = 0.4, linetype = "dashed") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = cfg$col_diff, alpha = 0.18) +
    geom_line(colour = cfg$col_diff, linewidth = 1) +
    labs(subtitle = sprintf("B  Difference (%s)", diff_name), x = x_axis_title(fit, m),
         y = "Difference in predicted pain") +
    theme_gamm() +
    theme(plot.subtitle = element_text(colour = INK, face = "bold"))

  stat_line <- if (sx$term_type == "smooth") {
    sprintf("Difference curve: edf = %.1f | F = %.1f | p %s (FDR %s) | AIC change vs single-curve model = %.1f",
            sx$diff_edf, sx$diff_f, fmt_p(sx$diff_p), fmt_p(sx$diff_p_fdr), sx$delta_aic)
  } else {
    sprintf("Slope difference (linear term): F = %.1f | p %s (FDR %s) | AIC change = %.1f",
            sx$diff_f, fmt_p(sx$diff_p), fmt_p(sx$diff_p_fdr), sx$delta_aic)
  }

  (pA | pB) +
    plot_annotation(
      title = sprintf("%s, %s", metric_label(m), unit_label(hit$group_value)),
      subtitle = stat_line,
      caption = paste0("Covariates at their mean; averaged over cap size; subject and experiment ",
                       "effects excluded. Bands: 95% CI.\nShaded x-ranges in B: CI excludes 0. ",
                       "The difference-curve p tests whether the curve shapes differ."),
      theme = theme_gamm()
    )
}

# =============================================================================
# RUN
# =============================================================================
results_con <- dbConnect(duckdb::duckdb(), dbdir = cfg$results_db_path, read_only = TRUE)
hits <- load_hits(results_con, cfg$feature_set, p_sig = cfg$hits_p_sig, decomposition = decomp)
sexdiff <- if (dbExistsTable(results_con, "gamm_sex_diff") &&
               "decomposition" %in% dbListFields(results_con, "gamm_sex_diff")) {
  as_tibble(dbGetQuery(results_con,
                       "SELECT * FROM gamm_sex_diff WHERE feature_set = ? AND decomposition = ?",
                       params = list(cfg$feature_set, decomp)))
} else tibble()
dbDisconnect(results_con, shutdown = TRUE)

if (nrow(hits) == 0) stop("No significant single-metric models to plot.")
message("Plotting ", nrow(hits), " significant model(s)",
        if (nrow(sexdiff) > 0) " (+ sex figures where available)" else
          " (no gamm_sex_diff table: run sex_diff_core.R for sex figures)")

set.seed(1)
n_main <- n_sex <- 0L
for (i in seq_len(nrow(hits))) {
  h <- hits[i, ]
  if (!h$rds_ok) { message("  [WARN] missing model file: ", h$rds_path); next }
  out_dir <- fig_unit_dir(cfg$fig_dir, h$level, h$group_value, decomp)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- sub("_z$", "", h$metrics)
  mod <- readRDS(h$rds_path)

  p <- tryCatch(plot_effect(h, mod), error = function(e) {
    message("  [WARN] ", h$group_value, " ", stem, ": ", conditionMessage(e)); NULL })
  if (!is.null(p)) {
    two <- inherits(p, "patchwork")
    ggsave(file.path(out_dir, paste0(stem, "_effect.png")), p,
           width = cfg$fig_width * if (two) 1.8 else 1,
           height = cfg$fig_height * if (two) 1.15 else 1, dpi = cfg$fig_dpi, bg = "white")
    n_main <- n_main + 1L
  }

  sx <- if (nrow(sexdiff) > 0) sexdiff %>% filter(level == h$level, group_value == h$group_value,
                                                   metric == h$metrics, status == "ok") else tibble()
  if (nrow(sx) == 1 && file.exists(sx$rds_path)) {
    fit <- readRDS(sx$rds_path)
    ps <- tryCatch(plot_sex(h, sx, fit), error = function(e) {
      message("  [WARN] sex figure ", h$group_value, " ", stem, ": ", conditionMessage(e)); NULL })
    if (!is.null(ps)) {
      ggsave(file.path(out_dir, paste0(stem, "_sex.png")), ps,
             width = cfg$fig_width * 1.8, height = cfg$fig_height * 1.15,
             dpi = cfg$fig_dpi, bg = "white")
      n_sex <- n_sex + 1L
    }
  }
}
message("Done: ", n_main, " effect figure(s), ", n_sex, " sex figure(s) under ", cfg$fig_dir)