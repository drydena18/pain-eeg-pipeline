# =============================================================================
# sex_diff_core.R
# =============================================================================
# Follow-up to the GAMM screen: does the metric -> pain curve differ between
# males and females? Run from R/gamm/ after source_core.R / channel_core.R.
#
# For every significant single-metric model (p < gamm_cfg$hits_p_sig, current
# feature_set, both levels), the screen model is refitted on EXACTLY the same
# rows (its own model frame) with sex as an ordered factor:
#
#   screen : ... + sex      + s(metric_z)
#   here   : ... + sex_o    + s(metric_z) + s(metric_z, by = sex_o)
#
# With an ordered factor, s(metric_z) is the curve for the reference sex and
# s(metric_z, by = sex_o) is the DIFFERENCE curve (other - reference). Its
# p-value tests "do the curve shapes differ?" directly; the sex_o term
# carries the overall offset. (An unordered `by = sex` would instead give
# two curves each tested against zero, which does not answer the question.)
# A metric that entered the screen linearly (< 3 unique values) gets
# metric_z:sex_o instead.
#
# With the between/within split (metric_decomposition = "between_within"),
# the difference curve is added for the PRIMARY part only (default: the
# subject-mean term, i.e. "does the metric -> pain-sensitivity curve differ
# by sex?"); the term is recorded in term_var.
#
# Reads:  results.duckdb  gamm_fitted; the screen .rds files
# Writes: results.duckdb  gamm_sex_diff (one row per hit)
#         <out_dir>/<level>/<unit>/sexdiff_<metric>.rds
#
# Multiple comparisons: diff_p is uncorrected; diff_p_fdr is
# Benjamini-Hochberg across all hits of this run (both levels together).
# =============================================================================

library(dplyr)
library(tibble)
library(mgcv)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/hit_helpers.R")

ref_sex <- gamm_cfg$sex_reference

results_con <- init_results_db(gamm_cfg$results_db_path)

# Derived entirely from gamm_fitted + the .rds files, so an older layout is
# simply rebuilt.
if (dbExistsTable(results_con, "gamm_sex_diff") &&
    !"term_var" %in% dbListFields(results_con, "gamm_sex_diff"))
  invisible(dbExecute(results_con, "DROP TABLE gamm_sex_diff"))

invisible(dbExecute(results_con, "
  CREATE TABLE IF NOT EXISTS gamm_sex_diff (
    fitted_at      VARCHAR,
    feature_set    VARCHAR,
    level          VARCHAR,
    group_value    VARCHAR,
    metric         VARCHAR,
    term_var       VARCHAR,
    status         VARCHAR,
    sex_reference  VARCHAR,
    sex_other      VARCHAR,
    n_subj_ref     INTEGER,
    n_subj_other   INTEGER,
    term_type      VARCHAR,
    diff_edf       DOUBLE,
    diff_f         DOUBLE,
    diff_p         DOUBLE,
    diff_p_fdr     DOUBLE,
    sex_offset     DOUBLE,
    sex_offset_p   DOUBLE,
    aic_screen     DOUBLE,
    aic_sex        DOUBLE,
    delta_aic      DOUBLE,
    rds_path       VARCHAR,
    PRIMARY KEY (feature_set, level, group_value, metric)
  )"))

hits <- load_hits(results_con, gamm_cfg$feature_set, p_sig = gamm_cfg$hits_p_sig)
if (nrow(hits) == 0) stop("No significant single-metric models in gamm_fitted (feature_set ",
                          gamm_cfg$feature_set, ", p < ", gamm_cfg$hits_p_sig, ").")
message("Sex-difference follow-up: ", nrow(hits), " significant model(s) (p < ",
        gamm_cfg$hits_p_sig, "); reference sex = ", ref_sex)

# -----------------------------------------------------------------------------
# One refit
# -----------------------------------------------------------------------------
refit_by_sex <- function(hit) {
  row <- tibble(fitted_at = as.character(Sys.time()), feature_set = gamm_cfg$feature_set,
                level = hit$level, group_value = hit$group_value, metric = hit$metrics,
                status = "ok", sex_reference = ref_sex, sex_other = NA_character_,
                n_subj_ref = NA_integer_, n_subj_other = NA_integer_, term_type = NA_character_,
                diff_edf = NA_real_, diff_f = NA_real_, diff_p = NA_real_, diff_p_fdr = NA_real_,
                sex_offset = NA_real_, sex_offset_p = NA_real_,
                aic_screen = NA_real_, aic_sex = NA_real_, delta_aic = NA_real_,
                rds_path = NA_character_, term_var = NA_character_)

  if (!hit$rds_ok) { row$status <- "rds_missing"; return(row) }
  mod <- readRDS(hit$rds_path)
  d <- mod$model
  # The term whose shape is compared between sexes: the primary part of a
  # split metric (subject mean by default), or the metric itself if pooled.
  m <- primary_var(mod, hit$metrics)
  row$term_var <- m

  if (!"sex" %in% names(d)) { row$status <- "no_sex_term"; return(row) }
  lv <- levels(droplevels(d$sex))
  if (length(lv) != 2 || !ref_sex %in% lv) {
    row$status <- "not_two_sexes"; return(row)
  }
  other <- setdiff(lv, ref_sex)
  row$sex_other <- other
  subj <- distinct(d, global_subjid, sex)
  row$n_subj_ref   <- sum(subj$sex == ref_sex)
  row$n_subj_other <- sum(subj$sex == other)
  if (min(row$n_subj_ref, row$n_subj_other) < gamm_cfg$sex_min_subjects) {
    row$status <- "too_few_subjects"; return(row)
  }

  d$sex_o <- factor(as.character(d$sex), levels = c(ref_sex, other), ordered = TRUE)
  contrasts(d$sex_o) <- contr.treatment(levels(d$sex_o))   # sex_o<other> = offset vs reference

  # Same terms as the screen model, with sex -> sex_o and the difference term added
  terms_screen <- attr(stats::terms(mod$formula), "term.labels")
  sm <- mod$smooth[vapply(mod$smooth, function(s) identical(s$term, m), logical(1))]
  if (length(sm) == 1) {
    k <- sm[[1]]$bs.dim
    diff_term <- sprintf("s(%s, by = sex_o, k = %d)", m, k)
    row$term_type <- "smooth"
  } else {
    diff_term <- sprintf("%s:sex_o", m)
    row$term_type <- "linear"
  }
  terms_new <- c(replace(terms_screen, terms_screen == "sex", "sex_o"), diff_term)
  f_sex <- stats::reformulate(terms_new, response = "pain_rating")

  fit <- tryCatch(
    bam(f_sex, data = d, method = "fREML", select = TRUE,
        discrete = !isFALSE(gamm_cfg$bam_discrete), nthreads = gamm_cfg$nthreads),
    error = function(e) { message("  [WARN] ", conditionMessage(e)); NULL })
  if (is.null(fit)) { row$status <- "fit_failed"; return(row) }

  s_sum <- summary(fit)
  off_name <- paste0("sex_o", other)
  if (off_name %in% rownames(s_sum$p.table)) {
    row$sex_offset   <- s_sum$p.table[off_name, 1]
    row$sex_offset_p <- s_sum$p.table[off_name, 4]
  }
  if (row$term_type == "smooth") {
    lab <- sprintf("s(%s):sex_o%s", m, other)
    if (lab %in% rownames(s_sum$s.table)) {
      row$diff_edf <- s_sum$s.table[lab, "edf"]
      row$diff_f   <- s_sum$s.table[lab, "F"]
      row$diff_p   <- s_sum$s.table[lab, "p-value"]
    }
  } else {
    lab <- grep(sprintf("^%s:sex_o|^sex_o.*:%s$", m, m), rownames(s_sum$p.table), value = TRUE)
    if (length(lab) == 1) {
      row$diff_edf <- 1
      row$diff_f   <- s_sum$p.table[lab, 3]^2
      row$diff_p   <- s_sum$p.table[lab, 4]
    }
  }
  row$aic_screen <- AIC(mod)
  row$aic_sex    <- AIC(fit)
  row$delta_aic  <- row$aic_sex - row$aic_screen

  fit$zscale <- mod$zscale
  fit$gamm_meta <- c(mod$gamm_meta, list(sex_reference = ref_sex, sex_other = other,
                                         screen_rds = hit$rds_path))
  out_dir <- dirname(hit$rds_path)
  row$rds_path <- file.path(out_dir, paste0("sexdiff_", sub("_z$", "", m), ".rds"))
  saveRDS(fit, row$rds_path)
  row
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
out <- vector("list", nrow(hits))
for (i in seq_len(nrow(hits))) {
  h <- hits[i, ]
  message(sprintf("  [%d/%d] %-8s %-12s %s", i, nrow(hits), h$level, h$group_value, h$metrics))
  out[[i]] <- refit_by_sex(h)
}
res <- bind_rows(out)
ok <- !is.na(res$diff_p)
res$diff_p_fdr[ok] <- p.adjust(res$diff_p[ok], method = "BH")

for (i in seq_len(nrow(res))) {
  invisible(dbExecute(results_con,
    "DELETE FROM gamm_sex_diff WHERE feature_set = ? AND level = ? AND group_value = ? AND metric = ?",
    params = list(res$feature_set[i], res$level[i], res$group_value[i], res$metric[i])))
}
invisible(dbAppendTable(results_con, "gamm_sex_diff", as.data.frame(res)))

message("\nStatus: ", paste(sprintf("%s = %d", names(table(res$status)), table(res$status)),
                            collapse = ", "))
sig <- res %>% filter(!is.na(diff_p), diff_p < 0.05) %>% arrange(diff_p)
message("Curves differing by sex (uncorrected p < 0.05): ", nrow(sig),
        " of ", sum(ok), "  |  BH-FDR < 0.05: ", sum(res$diff_p_fdr < 0.05, na.rm = TRUE))
for (i in seq_len(min(20, nrow(sig)))) {
  message(sprintf("    %-8s %-12s %-28s p %s  FDR %s  dAIC = %.1f",
                  sig$level[i], sig$group_value[i], sig$metric[i],
                  fmt_p(sig$diff_p[i]), fmt_p(sig$diff_p_fdr[i]), sig$delta_aic[i]))
}
close_results_db(results_con)
message("\nTable: ", gamm_cfg$results_db_path, " (gamm_sex_diff)")