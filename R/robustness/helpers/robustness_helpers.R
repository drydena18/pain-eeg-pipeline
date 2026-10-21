# =============================================================================
# robustness_helpers.R
# -----------------------------------------------------------------------------
# Shared machinery for robustness_core.R. Builds on helpers/hit_helpers.R
# (prepare_holdout, predict_population, boot_subjects, model_metric_vars)
# and helpers/gamm_combo_core.R (smooth_or_linear), so every robustness test
# prepares, refits and predicts exactly as holdout_replicate.R does.
#
#   refit_formula()       the screen model's formula, rebuilt for a data subset
#   refit_pair()          metric model + nested baseline refit on a subset
#   to_fit_frame()        id columns -> factors (for bam)
#   to_pred_frame()       id columns -> character (for predict_population)
#   score_predictions()   held-out error, metric vs baseline, trial + subject
#   permute_primary()     subject-level (or within-subject) shuffle of the
#                         metric's primary part
#   stratified_folds()    whole-subject folds stratified by cap size
#   init_robustness_table() / write_robustness_rows()
#
# V1 (Oct 2026): replaces cv_helpers.R / cv_split.R (drafts that never ran).
# =============================================================================

library(DBI)
library(dplyr)
library(tibble)
library(mgcv)

`%||%` <- function(a, b) if (is.null(a)) b else a

ID_COLS <- c("global_subjid", "experiment_id", "sex", "cap_size")

to_fit_frame <- function(df) {
  df <- as.data.frame(df)
  for (v in intersect(ID_COLS, names(df))) df[[v]] <- factor(as.character(df[[v]]))
  droplevels(df)
}

to_pred_frame <- function(df) {
  df <- as.data.frame(df)
  for (v in intersect(ID_COLS, names(df))) df[[v]] <- as.character(df[[v]])
  df
}

# -----------------------------------------------------------------------------
# refit_formula()
# The screen model's terms, rebuilt for `data`:
#   - smooths of one variable get k re-derived from this subset
#     (smooth_or_linear: a subset can have fewer unique laser levels);
#   - random-effect smooths and factor terms are dropped when the subset has
#     only one level (one experiment left, one cap size in a fold, ...);
#   - terms containing any of drop_vars are removed (nested baseline).
# te() tensor terms are kept as written.
# -----------------------------------------------------------------------------
refit_formula <- function(mod, data, k_max = 10, drop_vars = character(0)) {
  labs <- attr(stats::terms(mod$formula), "term.labels")
  out <- character(0)
  for (l in labs) {
    e <- str2lang(l)
    vars <- all.vars(e)
    if (length(drop_vars) > 0 && any(vars %in% drop_vars)) next
    fn <- if (is.call(e)) as.character(e[[1]]) else ""
    if (fn == "s" && grepl("bs\\s*=\\s*['\"]re['\"]", l)) {
      if (n_distinct(data[[vars[1]]]) > 1) out <- c(out, l)
      next
    }
    if (fn == "s" && length(vars) == 1) {
      out <- c(out, smooth_or_linear(vars, data[[vars]], k_max = k_max))
      next
    }
    if (fn == "" && length(vars) == 1 && is.factor(data[[vars]])) {
      if (nlevels(droplevels(data[[vars]])) > 1) out <- c(out, l)
      next
    }
    out <- c(out, l)
  }
  if (length(out) == 0) out <- "1"
  stats::reformulate(out, response = "pain_rating")
}

fit_bam <- function(f, data, cfg) {
  bam(f, data = data, method = "fREML", select = TRUE,
      discrete = isTRUE(cfg$bam_discrete), nthreads = cfg$nthreads)
}

# Metric model and nested baseline (metric terms removed), both refit on the
# same subset, so they differ only by the metric. The refits inherit the
# screen model's $zscale / $gamm_meta (same column scaling).
# drop_vars: extra terms removed from BOTH models (e.g. "cap_size" when the
# test experiment's cap size is absent from the refit data).
refit_pair <- function(mod, data, metric_z, cfg, baseline = TRUE, drop_vars = character(0)) {
  data <- to_fit_frame(data)
  mv <- model_metric_vars(mod, metric_z)
  m <- fit_bam(refit_formula(mod, data, cfg$k_max, drop_vars = drop_vars), data, cfg)
  m$zscale <- mod$zscale; m$gamm_meta <- mod$gamm_meta
  b <- if (baseline) fit_bam(refit_formula(mod, data, cfg$k_max, drop_vars = c(mv, drop_vars)), data, cfg) else NULL
  list(metric = m, base = b)
}

# -----------------------------------------------------------------------------
# score_predictions()
# Population-level predictions of the metric and baseline models for test
# subjects. Change in MSE = metric - baseline (negative = metric better),
# trial level (per-subject sums, so the bootstrap resamples subjects) and
# subject level (each subject's mean rating). 95% CI and one-sided p for
# "no improvement" from resampling test subjects (boot_subjects).
# -----------------------------------------------------------------------------
score_predictions <- function(test, pred_metric, pred_base, n_boot = 2000) {
  d <- tibble(global_subjid = as.character(test$global_subjid), y = test$pain_rating,
              pm = pred_metric, pb = pred_base) %>%
    filter(!is.na(y), !is.na(pm), !is.na(pb))
  if (nrow(d) == 0) return(NULL)
  ts <- d %>% group_by(global_subjid) %>%
    summarise(num = sum((y - pm)^2 - (y - pb)^2), den = n(), .groups = "drop")
  ss <- d %>% group_by(global_subjid) %>%
    summarise(obs = mean(y), pm = mean(pm), pb = mean(pb), .groups = "drop") %>%
    mutate(dd = (obs - pm)^2 - (obs - pb)^2)
  bt <- boot_subjects(ts$num, ts$den, n_boot)
  bs <- boot_subjects(ss$dd, rep(1, nrow(ss)), n_boot)
  mse_tb <- mean((d$y - d$pb)^2); mse_tm <- mean((d$y - d$pm)^2)
  mse_sb <- mean((ss$obs - ss$pb)^2); mse_sm <- mean((ss$obs - ss$pm)^2)
  safe_cor <- function(a, b) if (length(a) > 2 && stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_
  tibble(
    n_test_subj = nrow(ss), n_test_trials = nrow(d),
    trial_mse_base = mse_tb, trial_mse_metric = mse_tm,
    trial_pct = 100 * bt[["est"]] / mse_tb,
    trial_pct_lo = 100 * bt[["lo"]] / mse_tb, trial_pct_hi = 100 * bt[["hi"]] / mse_tb,
    trial_p = bt[["p"]],
    subj_mse_base = mse_sb, subj_mse_metric = mse_sm,
    subj_pct = 100 * bs[["est"]] / mse_sb,
    subj_pct_lo = 100 * bs[["lo"]] / mse_sb, subj_pct_hi = 100 * bs[["hi"]] / mse_sb,
    subj_p = bs[["p"]],
    subj_r_base = safe_cor(ss$obs, ss$pb), subj_r_metric = safe_cor(ss$obs, ss$pm)
  )
}

# Primary change in MSE (%) for a vector of metric predictions; baseline
# fixed. Used inside the permutation loop (no bootstrap).
primary_pct <- function(test, pred_metric, pred_base, primary = "subject") {
  d <- tibble(s = as.character(test$global_subjid), y = test$pain_rating,
              pm = pred_metric, pb = pred_base) %>% filter(!is.na(y), !is.na(pm), !is.na(pb))
  if (primary == "trial") return(100 * (mean((d$y - d$pm)^2) / mean((d$y - d$pb)^2) - 1))
  ss <- d %>% group_by(s) %>% summarise(o = mean(y), m = mean(pm), b = mean(pb), .groups = "drop")
  100 * (mean((ss$o - ss$m)^2) / mean((ss$o - ss$b)^2) - 1)
}

# -----------------------------------------------------------------------------
# permute_primary()
# One null draw for the metric's primary part on held-out rows `ho`:
#   between (or pooled): each subject's mean of the metric is replaced by
#     another subject's mean (whole subjects, never trials; within
#     experiment when within_exp = TRUE, because metric scales differ by
#     hardware); the subject's own trial deviations are kept.
#   within: trial deviations are shuffled within each subject.
# -----------------------------------------------------------------------------
permute_primary <- function(ho, metric_z, decomposition, primary_part, within_exp = TRUE) {
  subj <- as.character(ho$global_subjid)
  raw <- sub("_z$", "", metric_z)
  if (identical(decomposition, "between_within") && identical(primary_part, "within")) {
    w <- paste0(raw, "_within_z")
    ho[[w]] <- stats::ave(ho[[w]], subj, FUN = function(v) v[sample.int(length(v))])
    return(ho)
  }
  b <- stats::ave(ho[[metric_z]], subj, FUN = function(v) mean(v, na.rm = TRUE))
  key <- unique(data.frame(s = subj, e = as.character(ho$experiment_id), stringsAsFactors = FALSE))
  donor <- key$s
  groups <- if (within_exp) split(seq_len(nrow(key)), key$e) else list(seq_len(nrow(key)))
  for (g in groups) donor[g] <- key$s[g][sample.int(length(g))]
  names(donor) <- key$s
  b_by_s <- tapply(b, subj, function(v) v[1])
  b_new <- unname(b_by_s[donor[subj]])
  if (identical(decomposition, "between_within")) {
    ho[[paste0(raw, "_between_z")]] <- b_new
  } else {
    ho[[metric_z]] <- b_new + (ho[[metric_z]] - b)
  }
  ho
}

# Whole-subject folds, stratified by cap size (each cap size spread over all
# folds). Returns a named integer vector: fold per global_subjid.
stratified_folds <- function(subj_cap, k, seed) {
  set.seed(seed)
  subj_cap <- distinct(subj_cap)
  fold <- integer(nrow(subj_cap))
  for (cs in unique(subj_cap$cap_size)) {
    idx <- which(subj_cap$cap_size == cs)
    idx <- idx[sample.int(length(idx))]
    fold[idx] <- rep_len(seq_len(k), length(idx))
  }
  setNames(fold, subj_cap$global_subjid)
}

# -----------------------------------------------------------------------------
# Results table
# One row per hit x test x variant. Prediction tests fill the trial_* /
# subj_* columns; stat / stat_lo / stat_hi / p_value hold the test's headline
# number (see robustness_core.R header for what they mean per test).
# -----------------------------------------------------------------------------
ROBUSTNESS_TABLE_SQL <- "
  CREATE TABLE IF NOT EXISTS robustness (
    fitted_at        VARCHAR,
    feature_set      VARCHAR,
    decomposition    VARCHAR,
    level            VARCHAR,
    group_value      VARCHAR,
    metric           VARCHAR,
    test             VARCHAR,
    variant          VARCHAR,
    status           VARCHAR,
    n_train_subj     INTEGER,
    n_test_subj      INTEGER,
    n_test_trials    INTEGER,
    trial_mse_base   DOUBLE,
    trial_mse_metric DOUBLE,
    trial_pct        DOUBLE,
    trial_pct_lo     DOUBLE,
    trial_pct_hi     DOUBLE,
    trial_p          DOUBLE,
    subj_mse_base    DOUBLE,
    subj_mse_metric  DOUBLE,
    subj_pct         DOUBLE,
    subj_pct_lo      DOUBLE,
    subj_pct_hi      DOUBLE,
    subj_p           DOUBLE,
    subj_r_base      DOUBLE,
    subj_r_metric    DOUBLE,
    stat             DOUBLE,
    stat_lo          DOUBLE,
    stat_hi          DOUBLE,
    p_value          DOUBLE,
    note             VARCHAR,
    PRIMARY KEY (feature_set, decomposition, level, group_value, metric, test, variant)
  )"

init_robustness_table <- function(con) {
  ensure_derived_table(con, "robustness", c("decomposition", "subj_pct_lo"), ROBUSTNESS_TABLE_SQL)
}

write_robustness_rows <- function(con, rows) {
  if (is.null(rows) || nrow(rows) == 0) return(invisible(NULL))
  cols <- dbListFields(con, "robustness")
  rows <- as.data.frame(rows)
  for (cn in setdiff(cols, names(rows))) rows[[cn]] <- NA
  rows <- rows[, cols]
  keys <- unique(rows[, c("feature_set", "decomposition", "level", "group_value", "metric", "test")])
  for (i in seq_len(nrow(keys))) {
    invisible(dbExecute(con,
      "DELETE FROM robustness WHERE feature_set = ? AND decomposition = ? AND level = ?
         AND group_value = ? AND metric = ? AND test = ?",
      params = unname(as.list(keys[i, ]))))
  }
  invisible(dbAppendTable(con, "robustness", rows))
}