# =============================================================================
# gamm_combo_core.R
# -----------------------------------------------------------------------------
# Shared, combo-capable GAMM fitting engine.
#
# Used by gamm_screen.R (which source_core.R and channel_core.R call). A
# "model" is defined by a character vector of metric column names (length 1
# today, length 2-3 once you're ready for interaction screening). Swapping
# from single-metric to multi-metric combinations is a one-argument change
# (see generate_metric_combos()), not a rewrite.
#
# Public functions:
#   build_combo_formula()      construct a formula from a base + metric vector
#   generate_metric_combos()   enumerate which metric combinations to fit
#   fit_gamm_combo()           fit one bam(), return model + a result-row tibble
#   init_results_db()          open/create the DuckDB results table
#   write_result_row()         incremental upsert of one fit's result row
#   close_results_db()         clean disconnect
#
# Design notes:
#   - safe_k(): caps smooth k at (n_unique - 1), floor 3, ceiling supplied
#   - safe_bam(): tryCatch wrapper, returns NULL on failure rather than
#     halting the whole screen
#   - select = TRUE + fREML: double-penalty automatic smooth shrinkage
#   - age_z entered linearly, not as a smooth (between-subject, N too small
#     for a reliable nonlinear estimate at N ~ 30-40 per experiment)
#   - fixed/random effects included conditionally (nlevels > 1 / n_distinct
#     > 1) to avoid "contrasts not defined" / unidentifiable RE errors
#
# V2 (Oct 2026):
#   - Smooth-term stats are matched by the variables inside each term label
#     (exact), not by fixed-string substring. None of the 44 v3.1 names is a
#     substring of another today (checked), but the old code relied on that
#     staying true: one future name such as "slow_z" next to "erd_slow_z"
#     would silently report the wrong term's p-value in a 2-metric combo.
#   - gamm_fitted gains feature_set and n_subjects; feature_set is part of
#     the primary key. init_results_db() refuses to append to a pre-v3.1
#     gamm_fitted (no feature_set column) and asks for it to be archived.
#   - safe_bam(): nthreads and discrete are configurable; the warning names
#     the unit generically (ROI or channel).
#   - metric_decomposition = "between_within": each metric enters as its
#     subject mean (_between_z) + the trial deviation from it (_within_z).
#     edf / f_stat / p_value report the PRIMARY part (primary_part, default
#     between: the pain-sensitivity question); edf_other / f_other / p_other
#     the other part. Columns decomposition and primary_part added.
#   - Each saved model carries $zscale (mean / SD of every z-scored column)
#     and $gamm_meta (level, unit, metrics, feature_set).
#   - smooth_or_linear(): a term with < 3 unique values enters linearly.
#     safe_k() floors k at 3, so a two-level covariate (laser_power in a
#     two-level experiment) or a 0/1 metric made every model at that unit
#     fail with "fewer unique covariate combinations than ... degrees of
#     freedom".
# =============================================================================

library(dplyr)
library(stringr)
library(purrr)
library(mgcv)
library(tibble)
library(duckdb)
library(DBI)

# =============================================================================
# safe_k / safe_bam
# =============================================================================
safe_k <- function(x, k_max = 10, k_min = 3) {
  n_unique <- length(unique(na.omit(x)))
  max(k_min, min(k_max, n_unique - 1L))
}

# A smooth needs at least k_min (3) unique values. A covariate or metric with
# fewer (e.g. a two-level laser_power, a 0/1 metric) enters linearly instead
# of failing the whole model.
smooth_or_linear <- function(var, x, k_max = 10, k_min = 3) {
  n_unique <- length(unique(na.omit(x)))
  if (n_unique < k_min) return(var)
  sprintf("s(%s, k = %d)", var, as.integer(safe_k(x, k_max = k_max, k_min = k_min)))
}

safe_bam <- function(formula, data, nm, unit = NA_character_, nthreads = 4L,
                     discrete = TRUE) {
  tryCatch(
    bam(formula  = formula,
        data     = data,
        method   = "fREML",
        select   = TRUE,
        discrete = discrete,
        nthreads = nthreads),
    error = function(e) {
      message("  [WARN] Model ", nm, " (unit ", unit, ") failed: ", conditionMessage(e))
      NULL
    }
  )
}

# =============================================================================
# metric_parts()
# -----------------------------------------------------------------------------
# The model columns a metric enters as.
#   "pooled"          : the metric itself, "x_z"
#   "between_within"  : "x_between_z" = the subject's mean of x_z (who has
#                       high / low values) and "x_within_z" = x_z minus that
#                       mean (this trial vs the subject's own average).
# With a per-subject random intercept in the model, a pooled s(x_z) is
# estimated almost entirely from within-subject variation, so a metric that
# differs BETWEEN people (a trait-like biomarker of pain sensitivity) is
# missed. The split estimates the two separately.
# =============================================================================
metric_parts <- function(metric, decomposition = "pooled") {
  if (identical(decomposition, "between_within")) {
    raw <- sub("_z$", "", metric)
    c(between = paste0(raw, "_between_z"), within = paste0(raw, "_within_z"))
  } else {
    c(pooled = metric)
  }
}

# Adds <raw>_between_z (subject mean of metric_z over that subject's
# non-missing trials) and <raw>_within_z (metric_z - subject mean).
# Used identically for training data (screen) and held-out data, where a new
# person's mean comes from their own trials.
add_between_within <- function(df, metric_z, subj = "global_subjid") {
  raw <- sub("_z$", "", metric_z)
  b <- stats::ave(df[[metric_z]], as.character(df[[subj]]),
                  FUN = function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
  df[[paste0(raw, "_between_z")]] <- b
  df[[paste0(raw, "_within_z")]]  <- df[[metric_z]] - b
  df
}

# =============================================================================
# build_combo_formula()
# -----------------------------------------------------------------------------
# base_formula : pain_rating ~ ... (covariates + REs already resolved by the
#                caller via safe_k / level guards)
# metrics      : character vector of z-scored column names, e.g.
#                  c("pre_paf_cog_hz_z")                         single
#                  c("pre_paf_cog_hz_z", "erd_asym_z")           additive pair
#                  c("pre_pow_slow_alpha_z", "pre_pow_fast_alpha_z")
#                                       tensor pair if listed in tensor_pairs
# data         : used only to compute safe_k() per metric
# tensor_pairs : list of length-2 character vectors fit as te(a, b)
# k_max        : smooth k ceiling per metric
# =============================================================================
build_combo_formula <- function(base_formula, metrics, data,
                                tensor_pairs = list(), k_max = 10,
                                decomposition = "pooled") {
  if (length(metrics) == 0L) return(base_formula)

  is_tensor_pair <- length(metrics) == 2L &&
    any(map_lgl(tensor_pairs, function(p) setequal(p, metrics)))

  parts <- lapply(metrics, metric_parts, decomposition = decomposition)
  part_names <- names(parts[[1]])

  if (is_tensor_pair) {
    # One tensor per part: te(a_between, b_between) + te(a_within, b_within)
    term <- paste(vapply(part_names, function(pn) {
      v <- c(parts[[1]][[pn]], parts[[2]][[pn]])
      k_vals <- vapply(v, function(x) as.integer(safe_k(data[[x]], k_max = k_max)), integer(1))
      sprintf("te(%s, %s, k = c(%d, %d))", v[1], v[2], k_vals[1], k_vals[2])
    }, character(1)), collapse = " + ")
  } else {
    vars <- unlist(parts, use.names = FALSE)
    terms <- vapply(vars, function(v) smooth_or_linear(v, data[[v]], k_max = k_max),
                    character(1))
    term <- paste(terms, collapse = " + ")
  }

  update(base_formula, as.formula(paste(". ~ . +", term)))
}

# =============================================================================
# generate_metric_combos()
# -----------------------------------------------------------------------------
# Returns a tibble: model_name, metrics (list-column), n_metrics
# =============================================================================
generate_metric_combos <- function(metrics, sizes = 1L, include_baseline = TRUE) {
  combo_list <- list()

  if (include_baseline) {
    combo_list[["m00_baseline"]] <- character(0)
  }

  for (sz in sizes) {
    if (sz > length(metrics)) next
    combo_mat <- combn(metrics, sz)
    for (j in seq_len(ncol(combo_mat))) {
      combo <- combo_mat[, j]
      nm <- paste0("m_", sz, "metric_", paste(str_remove(combo, "_z$"), collapse = "_x_"))
      combo_list[[nm]] <- combo
    }
  }

  tibble(
    model_name = names(combo_list),
    metrics    = unname(combo_list),
    n_metrics  = map_int(combo_list, length)
  )
}

# =============================================================================
# smooth_term_variables()
# -----------------------------------------------------------------------------
# "s(pre_sf_ratio_z)" -> "pre_sf_ratio_z"; "te(a,b)" -> c("a", "b").
# mgcv labels smooths as <type>(<var>[,<var>...])[:<by>]; the by-variable
# part is dropped.
# =============================================================================
smooth_term_variables <- function(label) {
  inner <- sub("^[A-Za-z0-9_.]+\\((.*)\\)(:.*)?$", "\\1", label)
  trimws(strsplit(inner, ",", fixed = TRUE)[[1]])
}

# =============================================================================
# extract_combo_result_row()
# -----------------------------------------------------------------------------
# One summary row per fitted model: identity columns + per-metric smooth-term
# stats (edf, F, p-value) + whole-model fit stats. For combos with 2+
# metrics, stats are pulled for EACH metric and pipe-delimited, so a
# 2-metric row stays one row. In a tensor pair both metrics map to the same
# te() term and so report the same stats.
# =============================================================================
extract_combo_result_row <- function(mod, model_name, level, group_value, metrics,
                                     rds_path, feature_set, n_subjects,
                                     decomposition = "pooled", primary_part = "between") {
  sm <- summary(mod)
  s_tab <- sm$s.table

  term_vars <- if (is.null(s_tab)) list() else
    setNames(lapply(rownames(s_tab), smooth_term_variables), rownames(s_tab))

  p_tab <- sm$p.table

  get_term_stat <- function(metric, stat_col) {
    hit <- if (length(term_vars) == 0) character(0) else
      names(term_vars)[vapply(term_vars, function(v) metric %in% v, logical(1))]
    if (length(hit) > 0) return(as.numeric(s_tab[hit[1], stat_col]))
    # Metric entered linearly (too few unique values for a smooth): report
    # it from the parametric table on the same scale (edf 1, F = t^2).
    if (!is.null(p_tab) && metric %in% rownames(p_tab)) {
      t_val <- p_tab[metric, 3]
      return(switch(stat_col,
                    "edf"     = 1,
                    "F"       = as.numeric(t_val^2),
                    "p-value" = as.numeric(p_tab[metric, 4]),
                    NA_real_))
    }
    NA_real_
  }

  # Stats for one part of every metric in the combo, pipe-delimited
  part_stats <- function(part) {
    if (length(metrics) == 0L) return(list(edf = NA_character_, f = NA_character_, p = NA_character_))
    vars <- vapply(metrics, function(m) metric_parts(m, decomposition)[[part]], character(1))
    list(edf = paste(round(vapply(vars, get_term_stat, numeric(1), "edf"), 4), collapse = "|"),
         f   = paste(round(vapply(vars, get_term_stat, numeric(1), "F"),   4), collapse = "|"),
         p   = paste(signif(vapply(vars, get_term_stat, numeric(1), "p-value"), 4), collapse = "|"))
  }
  split <- identical(decomposition, "between_within")
  prim  <- part_stats(if (split) primary_part else "pooled")
  other <- if (split) part_stats(setdiff(c("between", "within"), primary_part)) else
    list(edf = NA_character_, f = NA_character_, p = NA_character_)

  # Computed before tibble(): inside tibble(), `metrics` would already refer
  # to the pipe-joined string column, and length() of that is always 1.
  n_metrics_val <- length(metrics)
  metrics_str   <- paste(metrics, collapse = "|")

  tibble(
    fitted_at   = as.character(Sys.time()),
    feature_set = feature_set,
    level       = level,
    group_value = group_value,
    model_name  = model_name,
    metrics     = metrics_str,
    n_metrics   = n_metrics_val,
    decomposition = decomposition,
    primary_part  = if (split) primary_part else "pooled",
    edf         = prim$edf,
    f_stat      = prim$f,
    p_value     = prim$p,
    edf_other   = other$edf,
    f_other     = other$f,
    p_other     = other$p,
    aic         = AIC(mod),
    bic         = BIC(mod),
    dev_expl    = sm$dev.expl,
    r_sq        = sm$r.sq,
    n_obs       = as.integer(nobs(mod)),
    n_subjects  = as.integer(n_subjects),
    rds_path    = rds_path
  )
}

# =============================================================================
# fit_gamm_combo()
# -----------------------------------------------------------------------------
# Builds the formula, fits via safe_bam(), saves the .rds and returns
# list(model = ..., result_row = tibble). The caller writes result_row
# immediately (incremental writes).
#
# level       : "source" or "channel"
# group_value : the ROI name (e.g. "Insula_lh") or channel label (e.g. "Cz")
# =============================================================================
fit_gamm_combo <- function(data, level, group_value, metrics, base_formula, model_name,
                           out_dir, tensor_pairs = list(), k_max = 10,
                           feature_set = NA_character_, nthreads = 4L,
                           discrete = TRUE, zscale = NULL,
                           decomposition = "pooled", primary_part = "between") {

  formula_use <- build_combo_formula(base_formula, metrics, data,
                                     tensor_pairs = tensor_pairs, k_max = k_max,
                                     decomposition = decomposition)

  # Rows with NA in any of this model's metric columns are dropped by bam();
  # count the subjects that actually contribute.
  vars <- unlist(lapply(metrics, metric_parts, decomposition = decomposition), use.names = FALSE)
  used <- if (length(vars) == 0) rep(TRUE, nrow(data)) else
    stats::complete.cases(data[, vars, drop = FALSE])
  n_subj <- n_distinct(data$global_subjid[used])

  mod <- safe_bam(formula_use, data, model_name, group_value,
                  nthreads = nthreads, discrete = discrete)
  if (is.null(mod)) return(NULL)

  # Stored on the model object so downstream scripts (plots, sex follow-up,
  # holdout) need only the .rds: original-unit scaling of every z-scored
  # column, and what the model is.
  if (!is.null(zscale)) {
    raw <- c("age", "laser_power", "trial_index", sub("_z$", "", metrics))
    mod$zscale <- zscale[intersect(raw, names(zscale))]
  }
  mod$gamm_meta <- list(level = level, group_value = group_value, model_name = model_name,
                        metrics = metrics, feature_set = feature_set,
                        decomposition = decomposition,
                        primary_part = if (identical(decomposition, "between_within")) primary_part else "pooled")

  rds_path <- file.path(out_dir, paste0(model_name, ".rds"))
  saveRDS(mod, rds_path)

  result_row <- extract_combo_result_row(mod, model_name, level, group_value, metrics,
                                         rds_path, feature_set, n_subj,
                                         decomposition = decomposition, primary_part = primary_part)

  list(model = mod, result_row = result_row)
}

# =============================================================================
# DuckDB results table — incremental write design
# -----------------------------------------------------------------------------
# init_results_db(): opens (or creates) results.duckdb and ensures
#   gamm_fitted exists. Stops if an existing gamm_fitted predates v3.1 (no
#   feature_set column): those rows were fit on differently-defined metrics
#   and must not be mixed with new ones.
#
# write_result_row(): appends ONE row immediately (delete-then-insert on the
#   primary key). Single-process only; if screens are ever parallelised
#   against the same file this needs a lock or a queue.
#
# close_results_db(): disconnects and shuts down the embedded instance.
# =============================================================================
init_results_db <- function(db_path) {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  con <- dbConnect(duckdb::duckdb(), dbdir = db_path)

  if (dbExistsTable(con, "gamm_fitted")) {
    cols <- dbListFields(con, "gamm_fitted")
    if (!all(c("feature_set", "decomposition") %in% cols)) {
      dbDisconnect(con, shutdown = TRUE)
      stop("results.duckdb holds a gamm_fitted table from an older version of the GAMM ",
           "stage (no feature_set / decomposition column).\nArchive it first, e.g. rename\n  ",
           db_path, "\nto results_pre_v3.1.duckdb, then re-run.")
    }
  }

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS gamm_fitted (
      fitted_at    VARCHAR,
      feature_set  VARCHAR,
      level        VARCHAR,
      group_value  VARCHAR,
      model_name   VARCHAR,
      metrics      VARCHAR,
      n_metrics    INTEGER,
      decomposition VARCHAR,
      primary_part VARCHAR,
      edf          VARCHAR,
      f_stat       VARCHAR,
      p_value      VARCHAR,
      edf_other    VARCHAR,
      f_other      VARCHAR,
      p_other      VARCHAR,
      aic          DOUBLE,
      bic          DOUBLE,
      dev_expl     DOUBLE,
      r_sq         DOUBLE,
      n_obs        INTEGER,
      n_subjects   INTEGER,
      rds_path     VARCHAR,
      PRIMARY KEY (feature_set, level, group_value, model_name)
    )
  ")
  con
}

write_result_row <- function(con, result_row) {
  if (is.null(result_row)) return(invisible(NULL))

  dbExecute(con,
    "DELETE FROM gamm_fitted
     WHERE feature_set = ? AND level = ? AND group_value = ? AND model_name = ?",
    params = list(result_row$feature_set, result_row$level,
                  result_row$group_value, result_row$model_name)
  )
  dbAppendTable(con, "gamm_fitted", as.data.frame(result_row))
  invisible(NULL)
}

close_results_db <- function(con) {
  dbDisconnect(con, shutdown = TRUE)
}