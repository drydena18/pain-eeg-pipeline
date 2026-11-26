# =============================================================================
# report_within_effect.R                                              V1.0.1
# -----------------------------------------------------------------------------
# One-shot report of the within-person (trial-to-trial) results, replacing
# the line-by-line queries: screen summary, within-person hits by window and
# unit, the pre-stimulus table for the cluster units, and - for one TARGET
# unit x metric - positive control, direction, shape, effect size in rating
# points, a figure of the within-person smooth, and the held-out result if
# holdout_replicate.R (holdout_part = "within") has run.
#
# V1.0.1: the text report is now always written (V1.0.0 left it empty when
# sourced, because the log was never closed).
#
# Run from R/gamm/:
#   Rscript report_within_effect.R                       # S1_contra / pre_sf_logratio_z
#   Rscript report_within_effect.R S1_contra pre_sf_balance_z
#   or in R:  report_unit <- "SM_furman"; source("report_within_effect.R")
#
# Reads:  results.duckdb (gamm_fitted, holdout_replication; read-only), the .rds
# Writes: <out_dir>/within_report/
#           screen_within_hits_by_window.csv   hits per window x unit
#           pre_window_clusters.csv            every pre_ metric x cluster unit
#           target_across_units.csv            the target metric in each unit
#           within_<unit>_<metric>.png         partial effect of the within term
#           within_<unit>_<metric>.txt         this report (console text)
#
# Units follow gamm_cfg (analysis_levels must include "channel";
# channel_units / cluster_names / channel_include), so stored results from an
# older unit selection or a different number of subjects are not mixed in.
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(dplyr)
  library(tibble)
  library(mgcv)
})

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/hit_helpers.R")
cfg <- gamm_cfg

args <- commandArgs(trailingOnly = TRUE)
if (!exists("report_unit"))   report_unit   <- if (length(args) >= 1) args[1] else "S1_contra"
if (!exists("report_metric")) report_metric <- if (length(args) >= 2) args[2] else "pre_sf_logratio_z"
if (!exists("open_results_ro"))
  open_results_ro <- function(path) dbConnect(duckdb::duckdb(), path, read_only = TRUE)

out_dir <- file.path(cfg$out_dir, "within_report")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stem <- sprintf("within_%s_%s", report_unit, sub("_z$", "", report_metric))
log_file <- file.path(out_dir, paste0(stem, ".txt"))
while (sink.number() > 0) sink()          # clear a sink left open by an earlier failed run
sink(log_file, split = TRUE)

# The report body runs inside a function so the log is always closed (and
# written to disk), even if a step stops with an error. on.exit() at the top
# level of a sourced script does nothing - V1.0.0 left the file empty.
report_body <- function() {

num <- function(x) suppressWarnings(as.numeric(x))
say <- function(...) cat(sprintf(...), "\n", sep = "")
rule <- function(title) cat("\n== ", title, " ", strrep("=", max(3, 70 - nchar(title))), "\n", sep = "")
# Cluster units carry an underscore (S1_contra, SM_furman); electrodes do not
cluster_units <- function(u) u[grepl("_", u)]

# -----------------------------------------------------------------------------
# Load screen results (current feature set, decomposition, channel units)
# -----------------------------------------------------------------------------
if (!"channel" %in% get_analysis_levels(cfg)) stop("analysis_levels must include 'channel'.")
rcon <- open_results_ro(cfg$results_db_path)
fits <- dbGetQuery(rcon, "
  SELECT level, group_value, metrics, p_value, edf, p_other, edf_other, f_other,
         n_subjects, n_obs, rds_path
  FROM gamm_fitted
  WHERE feature_set = ? AND decomposition = ? AND n_metrics = 1 AND level = 'channel'",
  params = list(cfg$feature_set, cfg$metric_decomposition)) %>%
  as_tibble() %>%
  filter_channel_units(cfg, "model") %>%
  mutate(p_b = num(p_value), p_w = num(p_other), edf_w = num(edf_other), f_w = num(f_other),
         n_subj = num(n_subjects), n_obs = num(n_obs),
         window = sub("_.*", "", metrics))
ho_tab <- if (dbExistsTable(rcon, "holdout_replication") &&
              "tested_part" %in% dbListFields(rcon, "holdout_replication"))
  as_tibble(dbGetQuery(rcon, "SELECT * FROM holdout_replication
                              WHERE feature_set = ? AND decomposition = ? AND tested_part = 'within'",
                       params = list(cfg$feature_set, cfg$metric_decomposition))) else NULL
dbDisconnect(rcon, shutdown = TRUE)
if (nrow(fits) == 0) stop("No single-metric channel-level models for feature_set ", cfg$feature_set,
                          " / ", cfg$metric_decomposition, " within channel_units = '", cfg$channel_units, "'.")

say("Within-person report  |  feature set %s  |  %s  |  channel_units = %s  |  %s",
    cfg$feature_set, cfg$metric_decomposition, cfg$channel_units, format(Sys.time(), "%Y-%m-%d %H:%M"))
say("Target: %s / %s", report_unit, report_metric)

# -----------------------------------------------------------------------------
rule("1. SCREEN SUMMARY")
# -----------------------------------------------------------------------------
n_mod <- nrow(fits)
say("%d single-metric models, units: %s", n_mod, paste(sort(unique(fits$group_value)), collapse = ", "))
ss <- fits %>% group_by(group_value) %>%
  summarise(subjects = paste(unique(range(n_subj)), collapse = "-"),
            trials = paste(unique(range(n_obs)), collapse = "-"), .groups = "drop")
print(as.data.frame(ss), row.names = FALSE)
if (n_distinct(fits$n_subj) > 1)
  say("NOTE: units were fitted on different numbers of subjects; compare hit rates only between units with the same N.")
say("Between part: p < .05 in %d (chance ~%.1f), p < .10 in %d (chance ~%.1f)",
    sum(fits$p_b < .05, na.rm = TRUE), .05 * n_mod, sum(fits$p_b < .10, na.rm = TRUE), .10 * n_mod)
say("Within part : p < .05 in %d (chance ~%.1f)", sum(fits$p_w < .05, na.rm = TRUE), .05 * n_mod)

# -----------------------------------------------------------------------------
rule("2. WITHIN-PERSON HITS BY WINDOW x UNIT (p < .05)")
# -----------------------------------------------------------------------------
hw <- fits %>% group_by(window, group_value) %>%
  summarise(hits = sum(p_w < .05, na.rm = TRUE), models = n(), .groups = "drop") %>%
  mutate(cell = sprintf("%d/%d", hits, models))
wide <- hw %>% select(window, group_value, cell) %>%
  tidyr::pivot_wider(names_from = group_value, values_from = cell)
print(as.data.frame(wide), row.names = FALSE)
write.csv(hw, file.path(out_dir, "screen_within_hits_by_window.csv"), row.names = FALSE)
say("Only pre_ effects can PRECEDE the stimulus; post_/whole_/delta_ effects may reflect the response to it.")

# -----------------------------------------------------------------------------
rule("3. PRE-STIMULUS METRICS x CLUSTER UNITS (within-person p)")
# -----------------------------------------------------------------------------
pre <- fits %>% filter(window == "pre", group_value %in% cluster_units(unique(fits$group_value)))
if (nrow(pre) > 0) {
  pw <- pre %>% mutate(metric = sub("^pre_|_z$", "", metrics), metric = sub("_z$", "", metric),
                       cell = sprintf("%.3g%s", p_w, ifelse(p_w < .05, "*", ""))) %>%
    select(metric, group_value, cell) %>%
    tidyr::pivot_wider(names_from = group_value, values_from = cell) %>% arrange(metric)
  print(as.data.frame(pw), row.names = FALSE)
  write.csv(pre %>% select(unit = group_value, metrics, p_between = p_b, p_within = p_w, edf_within = edf_w,
                           n_subjects = n_subj, n_obs),
            file.path(out_dir, "pre_window_clusters.csv"), row.names = FALSE)
  say("rel_slow / rel_fast / sf_logratio / sf_balance / sf_ratio are near-monotonic transforms of the")
  say("same slow:fast composition - several hits there count as ONE effect.")
} else say("(no cluster units with pre_ metrics in this selection)")

# -----------------------------------------------------------------------------
rule(sprintf("4. TARGET METRIC ACROSS UNITS: %s", report_metric))
# -----------------------------------------------------------------------------
ta <- fits %>% filter(metrics == report_metric) %>%
  transmute(unit = group_value, p_within = p_w, edf_within = edf_w, F_within = f_w,
            p_between = p_b, n_subjects = n_subj, n_obs) %>% arrange(p_within)
print(as.data.frame(ta), row.names = FALSE, digits = 3)
write.csv(ta, file.path(out_dir, "target_across_units.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
rule(sprintf("5. TARGET MODEL: %s / %s", report_unit, report_metric))
# -----------------------------------------------------------------------------
tr <- fits %>% filter(group_value == report_unit, metrics == report_metric)
if (nrow(tr) == 0) stop("Target ", report_unit, " / ", report_metric, " not in the screen results above.")
if (nrow(tr) > 1) say("NOTE: %d stored rows for the target; using the first.", nrow(tr))
rds <- tr$rds_path[1]
if (is.na(rds) || !file.exists(rds)) stop("Model file not found: ", rds)
m <- readRDS(rds)

parts <- metric_parts(report_metric, "between_within")
st <- summary(m)$s.table
lab_of <- function(v) smooth_label_for(m, v)
row_of <- function(lab) if (!is.na(lab) && lab %in% rownames(st)) st[lab, ] else
  c(edf = NA, Ref.df = NA, F = NA, `p-value` = NA)
lw <- lab_of(parts[["within"]]); lb <- lab_of(parts[["between"]]); ll <- lab_of("laser_power_z")
rw <- row_of(lw); rb <- row_of(lb); rl <- row_of(ll)

n_s <- if ("global_subjid" %in% names(m$model)) n_distinct(m$model$global_subjid) else NA
say("Fitted on %s subjects, %s trials; experiments (id: subjects): %s", n_s, nrow(m$model),
    if ("experiment_id" %in% names(m$model))
      paste(sprintf("%s: %d", names(tapply(m$model$global_subjid, m$model$experiment_id, n_distinct)),
                    tapply(m$model$global_subjid, m$model$experiment_id, n_distinct)), collapse = ", ") else "?")
say("Positive control  laser power : F = %.1f, edf = %.2f, p %s", rl[["F"]], rl[["edf"]], fmt_p(rl[["p-value"]]))
say("Between part (who)            : edf = %.2f, p %s", rb[["edf"]], fmt_p(rb[["p-value"]]))
say("Within part (trial-to-trial)  : edf = %.2f, F = %.2f, p %s", rw[["edf"]], rw[["F"]], fmt_p(rw[["p-value"]]))
say("Deviance explained (whole model): %.1f %%, adjusted R2 = %.3f", 100 * summary(m)$dev.expl, summary(m)$r.sq)

if (is.na(lw)) stop("No within-person smooth for ", report_metric, " in the model.")
dir_w <- term_direction(m, parts[["within"]])
shape <- if (is.na(rw[["edf"]])) "?" else if (rw[["edf"]] < 0.5) "shrunk to ~zero" else
  if (rw[["edf"]] < 1.5) "approximately linear" else "curved (edf >= 1.5)"
tt <- predict(m, type = "terms")[, lw]
q <- stats::quantile(tt, c(.05, .25, .5, .75, .95))
spread <- unname(q[["95%"]] - q[["5%"]])
res_sd <- sqrt(m$sig2)
say("Direction: %s (correlation of metric with its partial effect = %+.2f)",
    if (abs(dir_w) < 0.3) "no clear monotonic direction" else if (dir_w > 0) "INCREASING (higher value -> more pain)"
    else "DECREASING (higher value -> less pain)", dir_w)
say("Shape: %s", shape)
say("Partial effect across trials (rating points): 5%% %+.3f | 25%% %+.3f | 50%% %+.3f | 75%% %+.3f | 95%% %+.3f",
    q[[1]], q[[2]], q[[3]], q[[4]], q[[5]])
say("5th-95th percentile spread: %.3f points = %.2f residual SD (residual SD %.2f points)",
    spread, spread / res_sd, res_sd)
say("Rating scale assumed 0-10 (merge_cfg$rating_range).")

png(file.path(out_dir, paste0(stem, ".png")), width = 1600, height = 1200, res = 200)
plot(m, select = which(vapply(m$smooth, `[[`, "", "label") == lw), shade = TRUE, seWithMean = TRUE,
     rug = TRUE, xlab = sprintf("%s, within-person (z)", sub("_z$", "", report_metric)),
     ylab = "Partial effect on pain rating (0-10)",
     main = sprintf("%s: %s (within-person)  p %s", unit_label(report_unit),
                    metric_label(sub("_z$", "", report_metric)), fmt_p(rw[["p-value"]])))
abline(h = 0, lty = 2)
invisible(dev.off())
say("Figure: %s", file.path(out_dir, paste0(stem, ".png")))

# -----------------------------------------------------------------------------
rule("6. HELD-OUT TEST (holdout_part = 'within')")
# -----------------------------------------------------------------------------
hrow <- if (is.null(ho_tab)) NULL else ho_tab %>% filter(group_value == report_unit, metric == report_metric)
if (is.null(hrow) || nrow(hrow) == 0) {
  say("Not run yet for this target. Set holdout_part = \"within\" (and holdout_targets) in gamm_defaults.R,")
  say("then run holdout_replicate.R.")
} else {
  h <- hrow[1, ]
  say("Status: %s | pre-registered: %s | criterion: %s", h$status, h$pre_registered, h$primary_level)
  say("Held-out subjects: %s (%s trials)", h$n_holdout_subj, h$n_holdout_trials)
  say("Refit on held-out subjects only: p %s, edf = %.2f, direction %+.2f (training %+.2f)",
      fmt_p(h$refit_p), h$refit_edf, h$refit_direction, h$train_direction)
  say("Subject-centred prediction: %+.3f %% error [%+.3f, %+.3f], p %s",
      h$within_pct, 100 * h$within_dmse_lo / h$within_mse_base, 100 * h$within_dmse_hi / h$within_mse_base,
      fmt_p(h$within_p))
  say("REPLICATED: %s", if (isTRUE(as.logical(h$replicated))) "YES" else "no")
}

# -----------------------------------------------------------------------------
rule("7. ABSTRACT NUMBERS (check before use)")
# -----------------------------------------------------------------------------
pw_pow <- fits %>% filter(group_value == report_unit, grepl("^pre_pow_", metrics)) %>% pull(p_w)
say("Between-person: %d/%d models p < .05 (smallest p = %.3f).",
    sum(fits$p_b < .05, na.rm = TRUE), n_mod, min(fits$p_b, na.rm = TRUE))
say("Within-person %s %s: p %s, %s, ~%.2f rating points across trials (5th-95th percentile).",
    report_unit, report_metric, fmt_p(rw[["p-value"]]),
    if (dir_w > 0) "higher value -> more pain" else "higher value -> less pain", spread)
if (length(pw_pow)) say("Absolute pre-stimulus power in %s: within p %s.", report_unit,
                        paste(sprintf("%.2f", sort(pw_pow)), collapse = ", "))
say("Fitted on %s subjects.", n_s)
}
tryCatch(report_body(), finally = { while (sink.number() > 0) sink() })
message("Report written: ", log_file)