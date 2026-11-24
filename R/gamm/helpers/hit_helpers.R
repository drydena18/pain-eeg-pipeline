# =============================================================================
# hit_helpers.R
# -----------------------------------------------------------------------------
# Shared by sex_diff_core.R, plot_gamm_hits.R, holdout_replicate.R,
# compare_levels.R, ../validate/power_analysis.R and
# ../robustness/robustness_core.R.
#
#   load_hits()            significant single-metric models of the current
#                          feature set and decomposition, from gamm_fitted
#   filter_channel_units() keep the channel-level rows selected by
#                          channel_units / cluster_names / channel_include
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
# V5 (Nov 2026): within-person held-out testing.
#   - load_hits(part = "within") selects hits on the WITHIN part (p_other);
#     `targets` (as_targets()) restricts to pre-registered level x unit x
#     metric rows, whatever their screen p. p_test = the p of the tested part.
#   - refit_baseline(drop = ...) removes only the named model variables
#     (within-person test: drop x_within_z, keep x_between_z).
#   - within_person_errors(): subject-centred trial errors (does the model
#     rank a NEW person's trials correctly, ignoring their overall level?).
#   - refit_holdout(): the same model refit on the held-out subjects alone
#     (independent-sample p and direction of one term); term_direction().
#
# V4 (Nov 2026): prepare_holdout() fills missing age / sex with the model's
# own $covariate_fill (training median age, sex "U") when the screen used
# missing_covariates = "impute".
# load_hits() returns only levels in gamm_cfg$analysis_levels
# (an explicit `levels` argument is intersected with it) and applies
# gamm_cfg$channel_units /
# cluster_names / channel_include to channel-level hits (apply_channel_units()
# in db_helpers.R), so every follow-up script works on the same units as the
# screen - an older all-electrode screen stored in gamm_fitted no longer
# leaks electrode hits into a clusters-only follow-up. Pass cfg = NULL to
# get every stored hit. unit_label() also spells out cluster names.
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
# cfg: the GAMM config (defaults to the global gamm_cfg that every follow-up
# script sources). Its analysis_levels decides which levels are returned:
# levels = NULL means "all of analysis_levels"; levels given explicitly
# (robustness$levels, ...) are intersected with analysis_levels. Its
# channel_units / cluster_names / channel_include select the channel-level
# units. cfg = NULL: both levels, every unit.
# part   : "between" (default; p_value = the primary part) or "within"
#          (p_other, between/within models only).
# targets: NULL, or pre-registered tests (see as_targets()); only those rows
#          are returned, whatever their screen p.
load_hits <- function(results_con, feature_set, p_sig = 0.05,
                      levels = NULL, decomposition = NULL,
                      cfg = get0("gamm_cfg", envir = globalenv(), ifnotfound = NULL),
                      part = c("between", "within"), targets = NULL) {
  part <- match.arg(part)
  scope <- if (is.null(cfg)) ANALYSIS_LEVELS else get_analysis_levels(cfg)
  if (is.null(levels)) {
    levels <- scope
  } else {
    off <- setdiff(levels, scope)
    if (length(off) > 0)
      message("  load_hits: level(s) ", paste(off, collapse = ", "),
              " skipped (not in gamm_cfg$analysis_levels).")
    levels <- intersect(levels, scope)
  }
  if (length(levels) == 0) levels <- "__none__"   # valid SQL, returns no rows
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
           f_other_num = suppressWarnings(as.numeric(f_other)),
           p_test = if (part == "within") p_other_num else p_num,
           tested_part = part) %>%
    { if (is.null(targets)) filter(., !is.na(p_test), p_test < p_sig)
      else semi_join(., as_targets(targets), by = c("level", "group_value", "metrics")) } %>%
    left_join(base, by = c("level", "group_value", "decomposition")) %>%
    # dAIC is only meaningful when both models saw the same rows
    mutate(delta_aic = if_else(!is.na(n_obs_baseline) & n_obs == n_obs_baseline,
                               aic - aic_baseline, NA_real_),
           rds_ok = !is.na(rds_path) & file.exists(rds_path)) %>%
    arrange(level, group_value, p_test)
  if (part == "within") out <- out %>% filter(decomposition == "between_within")

  # -- Channel-level units: same selection as the screen ----
  out <- filter_channel_units(out, cfg, what = "hit")
  if (!is.null(targets)) {
    miss <- anti_join(as_targets(targets), out, by = c("level", "group_value", "metrics"))
    if (nrow(miss) > 0)
      warning("Pre-registered target(s) not found in gamm_fitted (or outside channel_units): ",
              paste(sprintf("%s/%s/%s", miss$level, miss$group_value, miss$metrics), collapse = ", "),
              call. = FALSE, immediate. = TRUE)
  }

  n_bad <- sum(out$bam_discrete %in% TRUE & out$decomposition == "between_within")
  if (n_bad > 0)
    warning(n_bad, " hit(s) come from discrete = TRUE fits with the between/within split; ",
            "their between-subject p-values are invalid. Re-run the screen with bam_discrete = FALSE.",
            call. = FALSE, immediate. = TRUE)
  out
}

# Pre-registered tests -> tibble(level, group_value, metrics).
# Accepts a data.frame with level / unit (or group_value) / metric (or
# metrics) columns, or a list of named vectors:
#   list(c(level = "channel", unit = "S1_contra", metric = "pre_sf_logratio_z"))
as_targets <- function(targets) {
  if (is.null(targets) || length(targets) == 0)
    return(tibble(level = character(), group_value = character(), metrics = character()))
  if (!is.data.frame(targets))
    targets <- dplyr::bind_rows(lapply(targets, function(x) as.list(x)))
  nm <- names(targets)
  pick <- function(a, b) if (a %in% nm) targets[[a]] else if (b %in% nm) targets[[b]] else
    stop("holdout_targets needs a '", a, "' column/field")
  tibble(level = as.character(pick("level", "level")),
         group_value = as.character(pick("unit", "group_value")),
         metrics = as.character(pick("metric", "metrics")))
}

# Keep only the channel-level rows the screen would select under cfg
# (channel_units, cluster_names, channel_include); source rows untouched.
# df needs `level` and `group_value` columns. cfg = NULL keeps everything.
filter_channel_units <- function(df, cfg = get0("gamm_cfg", envir = globalenv(), ifnotfound = NULL),
                                 what = "row") {
  if (is.null(cfg) || nrow(df) == 0 || !any(df$level == "channel")) return(df)
  ch <- unique(df$group_value[df$level == "channel"])
  keep_ch <- select_units(
    apply_channel_units(ch, cfg$channel_units %||% "all",
                        cfg$cluster_names %||% character(0), warn = FALSE),
    cfg$channel_include)
  n_drop <- sum(df$level == "channel" & !df$group_value %in% keep_ch)
  if (n_drop > 0)
    message("  ", n_drop, " channel-level ", what, "(s) outside channel_units = '",
            cfg$channel_units %||% "all", "' / channel_include ignored.")
  df[df$level != "channel" | df$group_value %in% keep_ch, , drop = FALSE]
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
  x <- sub("_rh$", " (right)", x)
  x <- sub("_contra$", " (contralateral)", x)
  x <- sub("_ipsi$", " (ipsilateral)", x)
  sub("^SM_furman$", "Sensorimotor (C3/Cz/C4)", x)
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

  # Missing age / sex filled exactly as for the training subjects (training
  # median age, sex "U"), if the model was fitted with missing_covariates =
  # "impute"; otherwise those rows are dropped below, as in the screen.
  df <- apply_covariate_fill(df, mod$covariate_fill)

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
# drop: model variables to remove (default: every part of the metric).
refit_baseline <- function(mod, metric_z, cfg, drop = NULL) {
  labs <- attr(stats::terms(mod$formula), "term.labels")
  mv <- if (is.null(drop)) model_metric_vars(mod, metric_z) else drop
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

# =============================================================================
# Within-person held-out testing (V5)
# =============================================================================
# Label of the smooth whose (single) variable is `var`; NA if none.
smooth_label_for <- function(mod, var) {
  hit <- vapply(mod$smooth, function(s) identical(s$term, var), logical(1))
  if (any(hit)) mod$smooth[[which(hit)[1]]]$label else NA_character_
}

# Overall direction of one smooth on the model's own rows: correlation of the
# variable with its partial effect (+1 = increasing, -1 = decreasing).
term_direction <- function(mod, var) {
  lab <- smooth_label_for(mod, var)
  if (is.na(lab) || !var %in% names(mod$model)) return(NA_real_)
  f <- predict(mod, type = "terms")[, lab]
  if (stats::sd(f) == 0) return(0)
  suppressWarnings(stats::cor(mod$model[[var]], f))
}

# Subject-centred trial errors. Each held-out subject's own mean is removed
# from the observed ratings AND from each model's predictions, so only the
# ORDERING of trials within a person is scored (a new person's overall pain
# level, which the within part cannot know, cancels). Both models get the
# same treatment. Returns per-subject sums for boot_subjects().
within_person_errors <- function(ho) {
  ho %>%
    group_by(global_subjid) %>%
    mutate(yc = pain_rating - mean(pain_rating),
           mc = pred_metric - mean(pred_metric),
           bc = pred_base - mean(pred_base)) %>%
    summarise(num = sum((yc - mc)^2 - (yc - bc)^2), den = n(),
              sse_m = sum((yc - mc)^2), sse_b = sum((yc - bc)^2), .groups = "drop")
}

# The screen model refit on the held-out subjects alone: an independent-
# sample test of one term. Terms whose variables do not vary in the held-out
# data (one sex level, one cap size, a single experiment, ...) are dropped.
refit_holdout <- function(mod, ho, var, cfg) {
  dat <- as.data.frame(ho)
  for (v in c("global_subjid", "experiment_id", "sex", "cap_size"))
    if (v %in% names(dat)) dat[[v]] <- factor(dat[[v]])
  varies <- function(v) {
    x <- dat[[v]]
    if (is.null(x)) return(FALSE)
    if (is.factor(x) || is.character(x)) return(length(unique(x[!is.na(x)])) > 1)
    isTRUE(stats::sd(x, na.rm = TRUE) > 0)
  }
  labs <- attr(stats::terms(mod$formula), "term.labels")
  keep <- labs[vapply(labs, function(l) all(vapply(all.vars(str2lang(l)), varies, logical(1))),
                      logical(1))]
  if (!any(vapply(keep, function(l) var %in% all.vars(str2lang(l)), logical(1))))
    stop("tested term not estimable in the held-out data")
  fit <- bam(stats::reformulate(keep, response = "pain_rating"), data = dat, method = "fREML",
             select = TRUE, discrete = FALSE, nthreads = if (is.null(cfg$nthreads)) 1 else cfg$nthreads)
  st <- summary(fit)$s.table
  lab <- smooth_label_for(fit, var)
  list(p = unname(st[lab, "p-value"]), edf = unname(st[lab, "edf"]),
       direction = term_direction(fit, var), n_subj = dplyr::n_distinct(dat$global_subjid),
       dropped = paste(setdiff(labs, keep), collapse = "; "))
}