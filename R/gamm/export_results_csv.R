# =============================================================================
# export_results_csv.R
# -----------------------------------------------------------------------------
# One tidy CSV of the GAMM results, for Excel: every single-metric model of
# the current feature_set / decomposition at the levels and units selected
# in gamm_defaults.R (analysis_levels, channel_units, channel_include),
# with its screen statistics, change in AIC vs the unit's baseline model,
# and the held-out replication result if holdout_replicate.R has run.
#
# Run from R/gamm/:   source("export_results_csv.R")
# Writes: <out_dir>/results_summary_<decomposition>.csv
#
# Columns: level, unit, unit_label, metric, metric_label, hit (screen
# p < hits_p_sig), p_between / edf_between / F_between (the tested part),
# p_within / edf_within / F_within, delta_aic (vs baseline, same rows; NA if
# the rows differ), dev_expl, r_sq, n_subjects, n_obs, holdout_status,
# holdout_pct_change (primary level, negative = better), holdout_p_fdr,
# replicated, rds_path.
#
# V1 (Nov 2026).
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(dplyr)
})

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/hit_helpers.R")

build_results_table <- function(con, cfg) {
  num <- function(x) suppressWarnings(as.numeric(x))
  decomp <- cfg$metric_decomposition
  fits <- dbGetQuery(con, sprintf("
    SELECT level, group_value, model_name, metrics, n_metrics, p_value, edf, f_stat,
           p_other, edf_other, f_other, dev_expl, r_sq, aic, n_obs, n_subjects, rds_path
    FROM gamm_fitted
    WHERE feature_set = ? AND decomposition = ? AND level IN (%s)",
    sql_in(con, get_analysis_levels(cfg))),
    params = list(cfg$feature_set, decomp)) %>% as_tibble()

  base <- fits %>% filter(model_name == "m00_baseline") %>%
    select(level, group_value, aic_base = aic, n_obs_base = n_obs)

  out <- fits %>%
    filter(n_metrics == 1) %>%
    filter_channel_units(cfg, what = "model") %>%
    left_join(base, by = c("level", "group_value")) %>%
    transmute(level, unit = group_value, unit_label = unit_label(group_value),
              metric = sub("_z$", "", metrics), metric_label = metric_label(metrics),
              hit = num(p_value) < cfg$hits_p_sig,
              p_between = num(p_value), edf_between = num(edf), F_between = num(f_stat),
              p_within = num(p_other), edf_within = num(edf_other), F_within = num(f_other),
              delta_aic = ifelse(!is.na(n_obs_base) & n_obs == n_obs_base, aic - aic_base, NA_real_),
              dev_expl, r_sq, n_subjects, n_obs, rds_path, metrics)

  # Held-out replication, if it has been run for this version
  if (dbExistsTable(con, "holdout_replication") &&
      "decomposition" %in% dbListFields(con, "holdout_replication")) {
    pct <- if (identical(cfg$holdout_primary, "trial")) "trial_pct" else "subj_pct"
    ho <- dbGetQuery(con, sprintf("
      SELECT level, group_value AS unit, metric AS metrics, status AS holdout_status,
             %s AS holdout_pct_change, primary_p_fdr AS holdout_p_fdr, replicated
      FROM holdout_replication WHERE feature_set = ? AND decomposition = ?", pct),
      params = list(cfg$feature_set, decomp)) %>% as_tibble()
    out <- out %>% left_join(ho, by = c("level", "unit", "metrics"))
  }

  out %>% select(-metrics) %>% relocate(rds_path, .after = last_col()) %>%
    arrange(level, unit, p_between)
}

# -- Run ----
con <- dbConnect(duckdb::duckdb(), gamm_cfg$results_db_path, read_only = TRUE)
tab <- build_results_table(con, gamm_cfg)
dbDisconnect(con, shutdown = TRUE)

out_file <- file.path(gamm_cfg$out_dir,
                      sprintf("results_summary_%s.csv", gamm_cfg$metric_decomposition))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, out_file, row.names = FALSE, na = "")

message(sprintf("%d models (%d hits at p < %.2g) -> %s", nrow(tab), sum(tab$hit, na.rm = TRUE),
                gamm_cfg$hits_p_sig, out_file))
print(head(as.data.frame(tab[order(tab$p_between),
                             c("unit", "metric", "p_between", "delta_aic", "dev_expl")]), 10),
      row.names = FALSE)