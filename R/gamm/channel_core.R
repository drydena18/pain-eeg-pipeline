# =============================================================================
# channel_core.R
# =============================================================================
# Channel-space GAMM combo screen: every screened alpha / ERD metric against
# every channel in spectral_metrics (restrict with
# gamm_cfg$channel_include). One model set per channel - channels are NOT
# pooled into a single model with s(channel, bs = "re"). Run from R/gamm/.
#
# Reads:  data.duckdb     (subjects, trials, spectral_metrics, subject_split,
#                          alignment_qc - read-only)
# Writes: results.duckdb  (gamm_fitted, incremental per-model writes)
#         <out_dir>/channel/<channel>/<model_name>.rds
#
# V2 (Oct 2026): the per-unit loop now lives in helpers/gamm_screen.R and
# is shared with source_core.R (the two copies had drifted; see that file).
# =============================================================================

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("helpers/gamm_combo_core.R")
source("helpers/gamm_screen.R")

run_gamm_screen("channel", gamm_cfg)