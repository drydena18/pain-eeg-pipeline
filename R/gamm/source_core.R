# =============================================================================
# source_core.R
# =============================================================================
# Source-space GAMM combo screen: every screened alpha / ERD metric against
# every ROI in source_metrics (18 ROIs in feature set v3.1: S1, S2, ACC,
# Insula, dlPFC, M1, each bilateral + _lh + _rh; restrict with
# gamm_cfg$roi_include). Run from R/gamm/.
#
# Reads:  data.duckdb     (subjects, trials, source_metrics, subject_split,
#                          alignment_qc - read-only)
# Writes: results.duckdb  (gamm_fitted, incremental per-model writes)
#         <out_dir>/source/<roi>/<model_name>.rds
#
# combo_sizes = 1L (single-metric models). Flip gamm_defaults.R's
# combo_sizes to c(1L, 2L) for pairwise screening - no code changes needed.
#
# Multiple-comparisons stance: screen uncorrected on training data; held-out
# replication is the correction (see gamm_defaults.R: apply_fdr_correction).
#
# No figures are made here: only .rds model objects and gamm_fitted rows.
# A plotting script that reads gamm_fitted for hits and renders them does
# not exist yet. (The old header named "export_significant_rois.R"; no such
# file was ever written.)
#
# V2 (Oct 2026): the per-unit loop now lives in helpers/gamm_screen.R and
# is shared with channel_core.R. Fixes carried by that move: `fir <-` typo
# (no result was ever written), missing ")" in s(trial_index_z, ...),
# connect_results_db() replaced by init_results_db() (creates gamm_fitted),
# explicit metric selection (p5_flag / phase / LEP are not screened).
# =============================================================================

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/gamm_screen.R")

run_gamm_screen("source", gamm_cfg)