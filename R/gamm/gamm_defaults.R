# =============================================================================
# gamm_defaults.R
# =============================================================================
# Configuration for the GAMM stage (source_core.R + channel_core.R).
# Mirrors the *_detault.m / *_default.py conventions used by preproc/spectral/
# source stages - one place to change paths, thresholds, and screening scope
# without touching the core scripts.
# =============================================================================

gamm_cfg <- list(

    # -- Database Paths ----
    data_db_path = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/data.duckdb",
    results_db_path = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/results.duckdb",

    # -- Output (model .rds objects still written to disk, per ROI/channel) ----
    out_dir = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/gamm_outputs",

    # -- QC Thresholds ----
    min_fooof_r2 = 0.80,
    min_trials_grp = 20L, # minimum trials after QC, per ROI or per channel

    # -- Combo Screening Scope ----
    # sizes = 1L -> single-metric models only
    # sizes = c(1L, 2L) -> add pairwise interaction screening
    # sizes = c(1L, 2L, 3L) -> add three-way combos
    combo_sizes = 1L,

    # Metric pairs to fit as te(a, b) tensor smooths instead of additive
    # s(a) + s(b), when both appear together in a 2-metric combo. Pairs not
    # listed here default to additive.
    tensor_pairs = list(
        c("pow_slow_z", "pow_fast_z")
    ),

    # -- Smooth Complexity Ceiling ----
    k_max = 10L,

    # -- Multiple Comparisons Stance ----
    # No correction applied at the screening stage (see source_core.R header
    # comment) - replication in the held-out CV split is the de facto
    # correction. This flag exists so that decision is visible in config,
    # not buried in code, in case it's revisited later.
    apply_fdr_correction = FALSE
)