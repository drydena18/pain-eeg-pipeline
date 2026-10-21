# =============================================================================
# hit_helpers.R
# -----------------------------------------------------------------------------
# Shared by sex_diff_core.R, plot_gamm_hits.R, holdout_replicate.R,
# compare_levels.R, ../validate/power_analysis.R and
# ../robustness/robustness_core.R.
#
#   load_hits()            significant single-metric models of the current
#                          feature set and decomposition, from gamm_fitted
#   ensure_derived_table() (re)create a follow-up table whose layout changed
#   fig_unit_dir()         figure folder for one unit (per decomposition)
#   metric_label()         "pre_paf_cog_hz" -> "Pre-stimulus PAF (CoG, Hz)"
#   unit_label()           "Insula_lh" -> "Insula (left)"
#   smooth_labels_re()     labels of the random-effect smooths in a model
#   marginal_lpmatrix()    population-level prediction matrix: covariates at
#                          their mean, random effects excluded, parametric
#                          factors (sex, cap size) averaged over the sample
#   predict_from_lp()      fit + SE + 95% CI from a prediction matrix
#   metric_grid()          prediction grid over the 1st-99th percentile
#
# All predictions use the model's own frame (mod$model), so they describe
# exactly the rows the model was fitted on.
#
# V3 (Oct 2026): load_hits() filters on decomposition (gamm_fitted now holds
# both) and warns about rows fitted with discrete = TRUE; bam_discrete
# defaults to FALSE in refit_baseline(); ensure_derived_table() and
# fig_unit_dir() added.
# =============================================================================

library(DBI)
library(dplyr)
library(tibble)

# =============================================================================
# Significant models
# =============================================================================
# decomposition: "between_within" / "pooled"; NULL = both (not recommended:
# the follow-up tables are keyed per decomposition).
load_hits <- function(results_con, feature_set, p_sig = 0.05,
                      levels = c("source", "channel"), decomposition = NULL) {
  dec_sql <- if (is.null(decomposition)) "" else
    sprintf("AND decomposition = %s", as.character(dbQuoteString(results_con, decomposition)))
  rows <- dbGetQuery(results_con, sprintf("
    SELECT level, group_value, model_name, metrics, decomposition, primary_part,
           edf, f_stat, p_value, edf_other, f_other, p_other,
           aic, dev_expl, r_sq, n_obs, n_subjects, bam_discrete, rds_path
    FROM gamm_fitted
    WHERE feature_set = ? AND n_metrics = 1
      AND level IN (%s) %s",
    paste(sprintf("'%s'", levels), collapse = ", "), dec_sql),
    params = list(feature_set))

  base <- dbGetQuery(results_con, sprintf("
    SELECT level, group_value, decomposition, aic AS aic_baseline, n_obs AS n_obs_baseline
    FROM gamm_fitted
    WHERE feature_set = ? AND model_name = 'm00_baseline' %s", dec_sql),
    params = list(feature_set))

  out <- rows %>%
    as_tibble() %>%
    mutate(bam_discrete = as.logical(bam_discrete),
           p_num = suppressWarnings(as.numeric(p_value)),
           edf_num = suppressWarnings(as.numeric(edf)),
           f_num = suppressWarnings(as.numeric(f_stat)),
           p_other_num = suppressWarnings(as.numeric(p_other)),
           edf_other_num = suppressWarnings(as.numeric(edf_other)),
           f_other_num = suppressWarnings(as.numeric(f_other))) %>%
    filter(!is.na(p_num), p_num < p_sig) %>%
    left_join(base, by = c("level", "group_value", "decomposition")) %>%
    # dAIC is only meaningful when both models saw the same rows
    mutate(delta_aic = if_else(!is.na(n_obs_baseline) & n_obs == n_obs_baseline,
                               aic - aic_baseline, NA_real_),
           rds_ok = !is.na(rds_path) & file.exists(rds_path)) %>%
    arrange(level, group_value, p_num)

  n_bad <- sum(out$bam_discrete %in% TRUE & out$decomposition == "between_within")
  if (n_bad > 0)
    warning(n_bad, " hit(s) come from discrete = TRUE fits with the between/within split; ",
            "their between-subject p-values are invalid. Re-run the screen with bam_discrete = FALSE.",
            call. = FALSE, immediate. = TRUE)
  out
}

# A follow-up table (gamm_sex_diff, holdout_replication, ...) is derived
# entirely from gamm_fitted and the .rds files, so a table with an older
# layout is simply dropped and recreated.
ensure_derived_table <- function(con, table, required_cols, create_sql) {
  if (dbExistsTable(con, table) && !all(required_cols %in% dbListFields(con, table))) {
    message("  Rebuilding ", table, " (older layout).")
    invisible(dbExecute(con, sprintf("DROP TABLE %s", table)))
  }
  invisible(dbExecute(con, create_sql))
}

# <fig_dir>/<level>/<unit>/ for between_within; .../<unit>/pooled/ for pooled
fig_unit_dir <- function(fig_dir, level, unit, decomposition = "between_within") {
  d <- file.path(fig_dir, level, unit)
  if (identical(decomposition, "pooled")) d <- file.path(d, "pooled")
  d
}

# =============================================================================
# Labels
# =============================================================================
METRIC_FAMILY_LABELS <- c(
  pow_slow_alpha  = "slow alpha power (8-10 Hz)",
  pow_fast_alpha  = "fast alpha power (10-12 Hz)",
  pow_alpha_total = "alpha power (8-12 Hz)",
  paf_cog_hz      = "PAF (centre of gravity, Hz)",
  rel_slow_alpha  = "relative slow alpha power",
  rel_fast_alpha  = "relative fast alpha power",
  sf_ratio        = "slow / fast ratio",
  sf_logratio     = "slow / fast log-ratio",
  sf_balance      = "slow-fast balance index",
  psi_cog         = "psi (balance x CoG)"
)
METRIC_WINDOW_LABELS <- c(whole = "Whole-epoch", pre = "Pre-stimulus",
                          post = "Post-stimulus", delta = "Change (post - pre) in")
METRIC_ERD_LABELS <- c(
  erd_slow            = "ERD, slow alpha",
  erd_fast            = "ERD, fast alpha",
  erd_pow_alpha_total = "ERD, alpha",
  erd_asym            = "ERD asymmetry (slow - fast)"
)

metric_label <- function(x) {
  x <- sub("_z$", "", x)
  part <- ifelse(grepl("_between$", x), " (subject mean)",
          ifelse(grepl("_within$", x), " (deviation from subject mean)", ""))
  x <- sub("_(between|within)$", "", x)
  paste0(vapply(x, function(m) {
    if (m %in% names(METRIC_ERD_LABELS)) return(unname(METRIC_ERD_LABELS[m]))
    w <- sub("_.*$", "", m)
    fam <- sub("^[a-z]+_", "", m)
    if (w %in% names(METRIC_WINDOW_LABELS) && fam %in% names(METRIC_FAMILY_LABELS)) {
      return(paste(METRIC_WINDOW_LABELS[[w]], METRIC_FAMILY_LABELS[[fam]]))
    }
    m
  }, character(1), USE.NAMES = FALSE), part)
}

# The model column that defines the hit: the primary part of a split metric
# (e.g. "x_between_z"), or the metric itself for a pooled model.
primary_var <- function(mod, metric) {
  meta <- mod$gamm_meta
  dec <- if (is.null(meta$decomposition)) "pooled" else meta$decomposition
  if (dec != "between_within") return(metric)
  metric_parts(metric, dec)[[meta$primary_part]]
}
other_var <- function(mod, metric) {
  meta <- mod$gamm_meta
  if (is.null(meta$decomposition) || meta$decomposition != "between_within") return(NA_character_)
  metric_parts(metric, meta$decomposition)[[setdiff(c("between", "within"), meta$primary_part)]]
}
model_metric_vars <- function(mod, metric) {
  meta <- mod$gamm_meta
  dec <- if (is.null(meta$decomposition)) "pooled" else meta$decomposition
  unname(metric_parts(metric, dec))
}

unit_label <- function(x) {
  x <- sub("_lh$", " (left)", x)
  sub("_rh$", " (right)", x)
}

# =============================================================================
# Prediction machinery
# =============================================================================
smooth_labels_re <- function(mod) {
  labs <- vapply(mod$smooth, function(s) s$label, character(1))
  is_re <- vapply(mod$smooth, function(s) inherits(s, "random.effect"), logical(1))
  labs[is_re]
}

# Grid of n points across the 1st-99th percentile of a model column, so a
# handful of extreme trials don't stretch the x-axis.
metric_grid <- function(mod, var, n = 200, probs = c(0.01, 0.99)) {
  q <- stats::quantile(mod$model[[var]], probs, na.rm = TRUE, names = FALSE)
  seq(q[1], q[2], length.out = n)
}

# grid           : data.frame with the varying column(s), e.g. metric_z
# fix            : named list of factor values to hold fixed, e.g.
#                  list(sex_o = "M")
# marginal_over  : parametric factors averaged over, weighted by their
#                  frequency in the model frame (unless fixed)
#
# Other numeric covariates are set to 0 (their mean: all are z-scored).
# Random-effect smooths are excluded, so the prediction is for a typical
# subject in a typical experiment, not for any particular one.
marginal_lpmatrix <- function(mod, grid, fix = list(),
                              marginal_over = c("sex", "sex_o", "cap_size")) {
  mf <- mod$model
  re <- smooth_labels_re(mod)
  re_vars <- unique(unlist(lapply(mod$smooth[vapply(mod$smooth, function(s)
    inherits(s, "random.effect"), logical(1))], function(s) s$term)))

  avg <- setdiff(intersect(marginal_over, names(mf)), names(fix))
  avg <- avg[vapply(avg, function(v) is.factor(mf[[v]]), logical(1))]
  combos <- if (length(avg) > 0) {
    mf %>% count(across(all_of(avg))) %>% mutate(w = n / sum(n))
  } else {
    tibble(w = 1)
  }

  build_nd <- function(combo) {
    nd <- grid
    for (v in setdiff(names(mf), c(names(grid), names(mf)[1]))) {
      col <- mf[[v]]
      if (is.factor(col)) {
        val <- if (v %in% names(fix)) fix[[v]]
               else if (v %in% names(combo)) as.character(combo[[v]])
               else levels(col)[1]                     # RE factors: excluded anyway
        nd[[v]] <- factor(rep(val, nrow(nd)), levels = levels(col), ordered = is.ordered(col))
        if (!is.null(attr(col, "contrasts"))) contrasts(nd[[v]]) <- attr(col, "contrasts")
      } else {
        nd[[v]] <- 0
      }
    }
    nd
  }

  X <- 0
  for (i in seq_len(nrow(combos))) {
    nd <- build_nd(combos[i, , drop = FALSE])
    Xi <- suppressWarnings(predict(mod, newdata = nd, type = "lpmatrix", exclude = re))
    X <- X + combos$w[i] * Xi
  }
  X
}

predict_from_lp <- function(mod, X) {
  fit <- as.numeric(X %*% stats::coef(mod))
  se  <- sqrt(pmax(0, rowSums((X %*% mod$Vp) * X)))
  tibble(fit = fit, se = se, lo = fit - 1.96 * se, hi = fit + 1.96 * se)
}

# z -> original units, using the scaling stored on the model by the screen.
# Returns list(x = values, raw = TRUE/FALSE).
# A within-subject deviation has no mean to add back: it converts to a
# deviation in the metric's original units.
to_original_units <- function(mod, metric_z, z) {
  is_within <- grepl("_within_z$", metric_z)
  raw_name <- sub("_(between|within)$", "", sub("_z$", "", metric_z))
  sc <- mod$zscale[[raw_name]]
  if (is.null(sc)) return(list(x = z, raw = FALSE))
  list(x = z * sc[["sd"]] + if (is_within) 0 else sc[["mean"]], raw = TRUE)
}

# Returns the relation too, for "p %s": "= 0.003", "= 4.1e-05", "< 2e-16"
# (mgcv reports p below machine precision as 0).
fmt_p <- function(p) {
  ifelse(is.na(p), "= NA",
  ifelse(p < 2e-16, "< 2e-16",
  ifelse(p < 0.001, paste("=", formatC(p, format = "e", digits = 1)),
         paste("=", formatC(p, format = "f", digits = 3)))))
}

# =============================================================================
# Figure style (shared by plot_gamm_hits.R and holdout_replicate.R)
# =============================================================================
INK      <- "#0b0b0b"
INK_2    <- "#52514e"
INK_3    <- "#8a8984"
GRIDLINE <- "#e6e5e1"
SEX_NAMES <- c(F = "Female", M = "Male")
sex_name <- function(s) ifelse(s %in% names(SEX_NAMES), SEX_NAMES[s], s)

theme_gamm <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background     = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major    = ggplot2::element_line(colour = GRIDLINE, linewidth = 0.3),
      axis.line           = ggplot2::element_line(colour = INK_3, linewidth = 0.3),
      axis.ticks          = ggplot2::element_line(colour = INK_3, linewidth = 0.3),
      axis.text           = ggplot2::element_text(colour = INK_2),
      axis.title          = ggplot2::element_text(colour = INK),
      plot.title          = ggplot2::element_text(colour = INK, face = "bold"),
      plot.subtitle       = ggplot2::element_text(colour = INK_2, margin = ggplot2::margin(b = 8)),
      plot.caption        = ggplot2::element_text(colour = INK_3, size = ggplot2::rel(0.75), hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position     = "none"
    )
}


# =============================================================================
# Held-out prediction (holdout_replicate.R)
# =============================================================================

# Population-level prediction for subjects the model has never seen.
#   - The subject random effect is always excluded (new people).
#   - The experiment random effect is kept for experiments seen in training
#     (it is known before a new person is tested) and excluded otherwise.
#   - Rows whose sex / cap size level was not in training get NA.
# predict() on a discrete bam() fails on unseen factor levels even when
# the term is excluded, so unseen levels are mapped to a training level
# first; with the term excluded that level has no effect on the prediction.
predict_population <- function(mod, newdata) {
  mf <- mod$model
  is_re <- vapply(mod$smooth, function(s) inherits(s, "random.effect"), logical(1))
  re_term <- function(v) {
    hit <- vapply(mod$smooth[is_re], function(s) identical(s$term, v), logical(1))
    if (any(hit)) mod$smooth[is_re][[which(hit)[1]]]$label else character(0)
  }
  subj_lab <- re_term("global_subjid")
  exp_lab  <- re_term("experiment_id")

  out <- rep(NA_real_, nrow(newdata))
  ok <- rep(TRUE, nrow(newdata))
  for (v in c("sex", "sex_o", "cap_size")) {
    if (v %in% names(mf) && is.factor(mf[[v]]))
      ok <- ok & as.character(newdata[[v]]) %in% levels(mf[[v]])
  }
  exp_seen <- if ("experiment_id" %in% names(mf))
    as.character(newdata$experiment_id) %in% levels(mf$experiment_id) else rep(TRUE, nrow(newdata))

  for (seen in unique(exp_seen[ok])) {
    idx <- which(ok & exp_seen == seen)
    nd <- newdata[idx, , drop = FALSE]
    for (v in names(mf)[-1]) {
      if (!is.factor(mf[[v]])) next
      val <- as.character(nd[[v]])
      val[!val %in% levels(mf[[v]])] <- levels(mf[[v]])[1]
      nd[[v]] <- factor(val, levels = levels(mf[[v]]), ordered = is.ordered(mf[[v]]))
    }
    excl <- c(subj_lab, if (!seen) exp_lab)
    out[idx] <- as.numeric(suppressWarnings(predict(mod, newdata = nd, exclude = excl,
                                                    type = "response")))
  }
  out
}

# Held-out rows for one unit x metric, prepared exactly like the screen's
# training rows (same QC), then z-scored with the TRAINING mean/SD stored on
# the model (mod$zscale). Returns NULL if the model predates $zscale.
prepare_holdout <- function(data_con, level, unit, metric_z, mod, cfg) {
  raw <- sub("_z$", "", metric_z)
  sc <- mod$zscale
  need <- c("age", "laser_power", "trial_index", raw)
  if (is.null(sc) || !all(need %in% names(sc))) return(NULL)

  load_cols <- raw
  if (isTRUE(cfg$exclude_p5_flag)) load_cols <- c(load_cols, "p5_flag")
  df <- load_metrics_wide(data_con, level, unit, metrics = load_cols, split = "holdout")
  if (nrow(df) == 0 || !raw %in% names(df)) return(df[0, ])

  df <- df %>%
    filter(!is.na(pain_rating), !is.na(laser_power), !is.na(trial_index),
           !is.na(global_subjid), !is.na(experiment_id), !is.na(age),
           !is.na(sex), !is.na(cap_size), !is.na(.data[[raw]]))
  if (isTRUE(cfg$exclude_p5_flag) && "p5_flag" %in% names(df)) {
    df <- df %>% filter(is.na(p5_flag) | p5_flag < 0.5)
  }
  if (is.numeric(cfg$min_fooof_r2)) {
    fq <- load_subject_metrics_wide(data_con, level, unit = unit, metrics = "fooof_r2",
                                    split = "holdout")
    if ("fooof_r2" %in% names(fq)) {
      bad <- fq$global_subjid[is.na(fq$fooof_r2) | fq$fooof_r2 < cfg$min_fooof_r2]
      df <- df %>% filter(!global_subjid %in% bad)
    }
  }

  z <- function(x, v) (x - sc[[v]][["mean"]]) / sc[[v]][["sd"]]
  df <- df %>%
    mutate(age_z = z(age, "age"),
           laser_power_z = z(laser_power, "laser_power"),
           trial_index_z = z(trial_index, "trial_index"),
           !!metric_z := z(.data[[raw]], raw),
           across(c(global_subjid, experiment_id, sex, cap_size), as.character))
  # Split: a held-out person's mean comes from their own trials, exactly as
  # for training subjects.
  if (identical(mod$gamm_meta$decomposition, "between_within"))
    df <- add_between_within(df, metric_z)
  df
}

# Refit the screen model without the metric on the SAME training rows, so
# metric and baseline models differ only by the metric term.
refit_baseline <- function(mod, metric_z, cfg) {
  labs <- attr(stats::terms(mod$formula), "term.labels")
  mv <- model_metric_vars(mod, metric_z)
  keep <- labs[!vapply(labs, function(l) any(mv %in% all.vars(str2lang(l))), logical(1))]
  f0 <- stats::reformulate(keep, response = "pain_rating")
  bam(f0, data = mod$model, method = "fREML", select = TRUE,
      discrete = isTRUE(cfg$bam_discrete), nthreads = cfg$nthreads)
}

# Cluster (subject) bootstrap of a mean difference.
#   num, den : per-subject numerator / denominator (trial level: sum of
#              squared-error differences / n trials; subject level: one
#              difference per subject and den = 1)
# Returns estimate, 95% CI and the one-sided p for "difference >= 0".
boot_subjects <- function(num, den, B = 2000) {
  n <- length(num)
  est <- sum(num) / sum(den)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B)
  bs <- rowSums(matrix(num[idx], nrow = B)) / rowSums(matrix(den[idx], nrow = B))
  c(est = est, lo = unname(stats::quantile(bs, 0.025)), hi = unname(stats::quantile(bs, 0.975)),
    p = (1 + sum(bs >= 0)) / (B + 1))
}