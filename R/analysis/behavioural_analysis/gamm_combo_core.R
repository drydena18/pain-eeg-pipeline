# =============================================================================
# gamm_combo_core.R
# =============================================================================
# Shared, combo-capable GAMM fitting engine.
#
# Used by both run_gamm_alpha_metrics_v2.R (channel_space) and run_gamm_source.R
# (source-space, looped over ROI). Neither caller script duplicates fitting
# logic, they handle data loading / QC / z-scoring / looping, and delegate
# every bam() call through this module.
#
# Core idea: a "model" is no longer a hardcoded formula. It's defined by a 
# character vector of metric column names. Swapping from single-metric to
# multi-metric combinations is a one-argument change (see generate_metric_combos()),
# not a rewrite.
#
# Public functions:
#   build_combo_formula()       construct a formula from a base + metric vector
#   generate_metric_combos()    enumerate which metric combination to fit
#   fit_gamm_combo()            fit one bam(), return model + a result-row tibble
#   init_results_db()           open/create the DuckDB results table
#   write_results_row()         incremental upsert of one fit's result row
#   close_results_db()          clean diconnect
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
    n_unique <- length(unique(x))
    max(k_min, min(k_max, n_unique - 1L))
}

safe_bam <- function(formula, data, nm, roi = NA_character_) {
  tryCatch(
    bam(formula  = formula,
        data     = data,
        method   = "fREML",
        select   = TRUE,
        discrete = TRUE,
        nthreads = 4L),
    error = function(e) {
      message("  [WARN] Model ", nm, " (ROI ", roi, ") failed: ", conditionMessage(e))
      NULL
    }
  )
}

# =============================================================================
# build_combo_formula()
# =============================================================================
# base_formula : an existing pain_rating - ... formula (covariates + REs
#                already resolved by the caller via safe_k/level-guards)
# metrics : character vector of z-scored column names, e.g.
#           c("CoG_pre_z")  -> single metric
#           c("CoG_pre_z", "delta_ERD_z")   -> additive pair
#           c("CoG_pre_z", "delta_ERD_z")   -> tensor pair, if also listed
#                                               tensor_pairs
# data : df used only to compute safe_k() per metric (not stored)
# tensor_pairs : optional list of length-2 chacter vectors. Any metric
#                pair present here is fit as te(a, b) instead of additive
#                s(a) + s(b).
# k_max : smooth k ceiling per metric (passed to safe_k)
#
# Returns a formula object ready for safe_bam()
# =============================================================================
build_combo_formula <- function(base_formula, metrics, data, tensor_pairs = list(), k_max = 10) {
    if (length(metrics) == 0L) return(base_formula)

    is_tensor_pair <- length(metrics) == 2L && any(map_lgl(tensor_pairs, function(p) setequal(p, metrics)))

    if (is_tensor_pair) {
        k_vals <- vapply(metrics, function(m) as.integer(safe_k(data[[m]], k_max = k_max)), integer(1))
        term <- sprintf("te(%s, %s, k = c(%d, %d))", metrics[1], metrics[2], k_vals[1], k_vals[2])
    } else {
        terms <- vapply(metrics, function(m) {
            k_m <- safe_k(data[[m]], k_max = k_max)
            sprintf("s(%s, k = %d)", m, k_m)
        }, character(1))
        term <- paste(terms, collapse = " + ")
    }

    update(base_formula, as.formula(paste(". ~ . +", term)))
}

# =============================================================================
# generate_metric_combos()
# =============================================================================
# metrics : character vector of all candidate metric column names available
#           for this dataset (caller filters to what actually exists/has
#           variance before passing it)
# sizes : integer vector of combination sizes to generate. Default 1L (single
#         metric). Flip to 1:2 or 1:3 and every downstream piece (formula building,
#         fitting, result table) handles it without modification.
# include_baseline : if TRUE, prepends a NULL-metrics entry (covariates
#                    only) as the comparison floor
#
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
            nm <- paste0("m_", sz, "metric_", paste(str_remove(combo, "_z$"), collapse = "_X"))
            combo_list[[nm]] <- combo
        }
    }

    tibble(
        model_name = names(combo_list),
        metrics = combo_list,
        n_metrics = map_int(combo_list, length)
    )
}

# =============================================================================
# extract_combo_result_row()
# =============================================================================
# Pulls one summary row per fitted model: identity columns + per-metric
# smooth-term stats (edf, F, p-value) + whole-model fit stats. For combos
# with 2+ metrics, smooth term stats are pulled for EACH metric in the 
# combo and pipe-delimited, so a 2-metric row stays one row (queryable),
# not two.
# =============================================================================
extract_combo_result_row <- function(mod, model_name, roi, metrics, rds_path) {
    sm <- summary(mod)
    s_tab <- sm$s.table

    get_term_stat <- function(metric, stat_col) {
        if (is.null(s_tab)) return(NA_real_)
        idx <- which(str_detect(rownames(s_tab), fixed(metric)))
        if (length(idx) == 0) return(NA_real_)
        as.numeric(s_tab[idx[1], stat_col])
    }

    if (length(metrics) == 0L) {
        edf_str <- f_str <- p_str <- NA_character_
    } else {
        edf_str <- paste(round(vapply(metrics, get_term_stat, numeric(1), "edf"), 4), collapse + "|")
        f_str <- paste(round(vapply(metrics, get_term_stat, numeric(1), "F"), 4), collapse = "|")
        p_str <- paste(signif(vapply(metrics, get_term_stat, numeric(1), "p-value"), 4), collapse = "|")
    }

    n_metrics_val <- length(metrics)
    metrics_str <- paste(metrics, collapse = "|")

    tibble(
        fitted_at = as.character(Sys.time()),
        roi = roi,
        model_name = model_name,
        metrics = metrics_str,
        n_metrics = n_metrics_val,
        edf = edf_str,
        f_stat = f_str,
        p_value = p_str,
        aic = AIC(mod),
        bic = BIC(mod),
        dev_expl = sm$dev.expl,
        r_sq = sm$r.sq,
        n_obs = nobs(mod),
        rds_path = rds_path
    )
}

# =============================================================================
# fit_gamm_combo()
# =============================================================================
# One-stop call: builds the formula, fits via safe_bam(), saves the .rds
# (skipped if NULL), and returns list(model = ..., result_row = tibble).
# Caller is responsible for calling write_result_row() with the result_row
# immediately after, kept separate so the caller controls transaction
# timing (incremental writes per design)
# =============================================================================
fit_gamm_combo <- function(data, roi, metrics, base_formula, model_name, out_dir, tensor_pairs = list(), k_max = 10) {
    formula_use <- build_combo_formula(base_formula, metrics, data, tensor_pairs = tensor_pairs, k_max = k_max)

    mod <- safe_bam(formula_use, data, model_name, roi)
    if (is.null(mod)) return(NULL)

    rds_path <- file.path(out_dir, paste0(model_name, ".rds"))
    saveRDS(mod, rds_path)

    result_row <- extract_combo_result_row(mod, model_name, roi, metrics, rds_path)

    list(model = mod, result_row = result_row)
}

# =============================================================================
# DuckDB results table - incremental write design
# =============================================================================
# init_results_db(): opens (or creates) the DuckDB file and ensures the
#   gamm_fitted table exists with the right schema. Return the connection:
#   caller keeps it open for the duration of the run and passes it to
#   write_result_row() after every successful fit.
#
# write_result_row(): appends ONE row immediately. DuckDB's single-file
#   writer handles this safely for single-prcess run.
#
# close_results_db(): disconnects cleanly, including shuting down the
#   embedded DuckDB instance (dbDisconnect(..., shutdown = TRUE)).
# =============================================================================
init_results_db <- function(db_path) {
    dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
    con <- dbConnect(duckdb::duckdb(), dbdir = db_path)

    dbExecute(con, "
    CREATE TABLE IF NOT EXISTS gamm_fitted (
      fitted_at   VARCHAR,
      roi         VARCHAR,
      model_name  VARCHAR,
      metrics     VARCHAR,
      n_metrics   INTEGER,
      edf         VARCHAR,
      f_stat      VARCHAR,
      p_value     VARCHAR,
      aic         DOUBLE,
      bic         DOUBLE,
      dev_expl    DOUBLE,
      r_sq        DOUBLE,
      n_obs       INTEGER,
      rds_path    VARCHAR,
      PRIMARY KEY (roi, model_name)
    )
  ")
  con
}

write_result_row <- function(con, result_row) {
  if (is.null(result_row)) return(invisible(NULL))
 
  # Upsert: delete any existing (roi, model_name) row, then insert.
  # Re-running a screen for the same ROI/model overwrites cleanly rather
  # than erroring on the primary key or silently duplicating.
  dbExecute(con,
    "DELETE FROM gamm_fitted WHERE roi = ? AND model_name = ?",
    params = list(result_row$roi, result_row$model_name)
  )
  dbWriteTable(con, "gamm_fitted", result_row, append = TRUE)
  invisible(NULL)
}

close_results_db <- function(con) {
    dbDisconnect(con, shutdown = TRUE)
}