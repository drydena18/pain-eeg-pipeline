# =============================================================================
# diagnose_screen.R
# -----------------------------------------------------------------------------
# Why are there no (or few) hits? Read-only summary of the current screen
# (feature_set, decomposition, analysis_levels / channel_units from
# gamm_defaults.R) and of the data it was fitted on. Run from R/gamm/:
#   source("diagnose_screen.R")
#
#   1. p-value and edf distribution of the tested (between) and the other
#      (within) part. Under no effect at all ~5% of p are < 0.05; far fewer
#      means p-values are piled up near 1 (terms shrunk to zero by
#      select = TRUE, or absorbed by the subject / experiment terms).
#   2. Training subjects per experiment / cap size: power for a
#      between-subject effect depends on the number of SUBJECTS, not trials.
#   3. Alignment status (subjects silently excluded by the merge).
#   4. Full summary() of the model with the smallest between-subject p.
# =============================================================================

suppressPackageStartupMessages({ library(DBI); library(duckdb); library(dplyr); library(mgcv) })
source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/hit_helpers.R")
cfg <- gamm_cfg
num <- function(x) suppressWarnings(as.numeric(x))
qs  <- function(x) round(stats::quantile(x, c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE), 3)

# -- 1. Screen results ----
rcon <- dbConnect(duckdb::duckdb(), cfg$results_db_path, read_only = TRUE)
fits <- dbGetQuery(rcon, sprintf("
  SELECT level, group_value, metrics, p_value, edf, p_other, edf_other, n_subjects, n_obs,
         bam_discrete, rds_path
  FROM gamm_fitted
  WHERE feature_set = ? AND decomposition = ? AND n_metrics = 1 AND level IN (%s)",
  sql_in(rcon, get_analysis_levels(cfg))),
  params = list(cfg$feature_set, cfg$metric_decomposition)) %>%
  as_tibble() %>% filter_channel_units(cfg, "model") %>%
  mutate(p_b = num(p_value), edf_b = num(edf), p_w = num(p_other), edf_w = num(edf_other))
dbDisconnect(rcon, shutdown = TRUE)

cat("\n== 1. SCREEN (", nrow(fits), " single-metric models; units: ",
    paste(unique(fits$group_value), collapse = ", "), ") ==\n", sep = "")
cat("Fitted with discrete = TRUE (p-values invalid): ", sum(as.logical(fits$bam_discrete), na.rm = TRUE), "\n", sep = "")
cat("\nBETWEEN-subject part (tested):\n")
cat("  p < 0.05: ", sum(fits$p_b < .05, na.rm = TRUE), "   p < 0.10: ", sum(fits$p_b < .10, na.rm = TRUE),
    "   (expected by chance alone: ~", round(.05 * nrow(fits), 1), " / ~", round(.10 * nrow(fits), 1), ")\n", sep = "")
cat("  p > 0.90: ", sum(fits$p_b > .90, na.rm = TRUE), "   NA: ", sum(is.na(fits$p_b)), "\n", sep = "")
cat("  p quantiles  : "); print(qs(fits$p_b))
cat("  edf quantiles: "); print(qs(fits$edf_b))
cat("  edf < 0.01 (term shrunk to zero): ", sum(fits$edf_b < 0.01, na.rm = TRUE), "\n", sep = "")
cat("\nWITHIN-subject part:\n")
cat("  p < 0.05: ", sum(fits$p_w < .05, na.rm = TRUE), "\n", sep = "")
cat("  p quantiles  : "); print(qs(fits$p_w))
cat("\nSubjects / trials per model: subjects ", paste(range(fits$n_subjects), collapse = "-"),
    ", trials ", paste(range(fits$n_obs), collapse = "-"), "\n", sep = "")

# -- 2-3. Data ----
dcon <- connect_data_db(cfg$data_db_path, read_only = TRUE)
units <- unique(fits$group_value[fits$level == "channel"])
if (length(units) > 0) {
  cat("\n== 2. TRAINING SUBJECTS WITH CLUSTER / CHANNEL DATA ==\n")
  print(as.data.frame(dbGetQuery(dcon, sprintf("
    SELECT s.experiment_name, s.cap_size, sp.split_group, COUNT(DISTINCT s.global_subjid) AS n_subjects
    FROM subjects s
    JOIN subject_split sp USING (global_subjid)
    WHERE s.global_subjid IN (SELECT DISTINCT global_subjid FROM spectral_metrics
                              WHERE channel IN (%s))
    GROUP BY 1, 2, 3 ORDER BY 1, 3", sql_in(dcon, units)))), row.names = FALSE)
}
if (dbExistsTable(dcon, "alignment_qc")) {
  cat("\n== 3. MERGE ALIGNMENT STATUS (all subjects) ==\n")
  print(as.data.frame(dbGetQuery(dcon,
    "SELECT status, COUNT(*) AS n FROM alignment_qc GROUP BY status ORDER BY n DESC")), row.names = FALSE)
}
disconnect_db(dcon)

# -- 4. Best model ----
best <- fits %>% filter(!is.na(p_b), file.exists(rds_path)) %>% slice_min(p_b, n = 1, with_ties = FALSE)
if (nrow(best) == 1) {
  cat("\n== 4. MODEL WITH THE SMALLEST BETWEEN-SUBJECT p: ", best$group_value, " / ", best$metrics, " ==\n", sep = "")
  print(summary(readRDS(best$rds_path)))
}