# =============================================================================
# gamm_defaults.R
# =============================================================================
# Configuration for the GAMM stage (source_core.R + channel_core.R).
# Mirrors the *_default.m / *_default.py conventions used by preproc/spectral/
# source stages - one place to change paths, thresholds, and screening scope
# without touching the core scripts.
#
# V2 (Oct 2026), feature set v3.1:
#   - Metrics are now selected by name pattern (metric_include) instead of
#     "every column the pivot returns". The v3.1 tables hold 58 trial-level
#     columns per unit; only the alpha families and ERD are screened.
#     p5_flag (a QC flag), phase (circular) and LEP (a different signal)
#     are deliberately excluded - see metric_include below.
#   - tensor_pairs renamed to the v3.1 metric names.
#   - min_fooof_r2 now refers to the SUBJECT-level FOOOF fit (FOOOF is fit
#     once per subject x unit on the trial-mean PSD); default off, because
#     no screened metric is FOOOF-derived.
#   - feature_set is stamped on every gamm_fitted row.
#   - Follow-up, held-out and figure settings for sex_diff_core.R,
#     holdout_replicate.R and plot_gamm_hits.R.
# =============================================================================

gamm_cfg <- list(

    # -- Database Paths ----
    data_db_path    = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/data.duckdb",
    results_db_path = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/results.duckdb",

    # -- Output (model .rds objects still written to disk, per ROI/channel) ----
    out_dir = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/gamm_outputs",

    # Must match merge_cfg$feature_set. Stamped on every gamm_fitted row so
    # results from different feature definitions can never be mixed.
    feature_set = "v3.1",

    # -- Metric Selection ----
    # A metric is screened if its name matches ANY include pattern and NO
    # exclude pattern. 10 families x 4 windows + 4 ERD = 44 metrics.
    #
    #   whole_/pre_/post_/delta_ x {pow_slow_alpha, pow_fast_alpha,
    #     pow_alpha_total, paf_cog_hz, rel_slow_alpha, rel_fast_alpha,
    #     sf_ratio, sf_logratio, sf_balance, psi_cog}
    #   erd_slow, erd_fast, erd_pow_alpha_total, erd_asym
    #
    # Not screened:
    #   p5_flag           QC flag (0/1), used as a row filter, not a predictor
    #   *phase*           circular; s(phase) on a linear axis is wrong. Use
    #                     sin_phase/cos_phase in a dedicated model if needed.
    #   n2_* / p2_* / n2p2_*  LEP amplitudes/latencies - separate question
    metric_include = c(
        "^(whole|pre|post|delta)_(pow_slow_alpha|pow_fast_alpha|pow_alpha_total|paf_cog_hz|rel_slow_alpha|rel_fast_alpha|sf_ratio|sf_logratio|sf_balance|psi_cog)$",
        "^erd_(slow|fast|pow_alpha_total|asym)$"
    ),
    metric_exclude = character(0),

    # -- Between / within split ----
    # "between_within" (default): each metric enters as
    #     s(<metric>_between_z) : the subject's mean over their trials
    #                             (who has high / low values)
    #   + s(<metric>_within_z)  : this trial's deviation from that mean
    # The per-subject random intercept absorbs differences between people,
    # so a single pooled s(<metric>_z) is estimated almost entirely from
    # within-person variation and misses trait-like biomarkers. In a
    # synthetic test a true between-subject effect scored p = 0.67 pooled
    # and p < 2e-16 split.
    # "pooled": the previous behaviour, one s(<metric>_z).
    metric_decomposition = "between_within",
    # Which part's p-value goes into gamm_fitted (and so defines a "hit").
    # "between" answers "does a person's typical value predict their pain
    # sensitivity?"; "within" answers "does a trial's deviation predict
    # that trial's pain?". The other part is stored in *_other columns.
    primary_part = "between",

    # -- Unit Selection ----
    # NULL = every ROI / channel present in the database. The source tables
    # hold 18 ROIs: 6 bilateral + 6 _lh + 6 _rh. Screening all three is
    # three partly-redundant tests per region; restrict here if wanted, e.g.
    # roi_include = c("S1", "S2", "ACC", "Insula", "dlPFC", "M1")   (bilateral)
    # roi_include = "_(lh|rh)$"                                      (regex)
    roi_include     = NULL,
    channel_include = NULL,

    # -- QC ----
    # p5_flag = 1 marks trials whose pre-stimulus slow or fast alpha power is
    # below that unit's 5th percentile (ERD denominator near zero). TRUE
    # drops those trials before fitting (all metrics, not just ERD).
    exclude_p5_flag = FALSE,

    # Subject-level FOOOF fit quality. If numeric, a subject's data for a
    # unit is dropped when that unit's fooof_r2 (source_subject_metrics /
    # spectral_subject_metrics) is below the threshold. NULL = no filter.
    min_fooof_r2 = NULL,

    min_trials_grp = 20L, # minimum trials after QC, per ROI or per channel

    # -- Combo Screening Scope ----
    # sizes = 1L -> single-metric models only
    # sizes = c(1L, 2L) -> add pairwise interaction screening
    #   (44 metrics -> 946 pairs per unit; consider narrowing metric_include)
    combo_sizes = 1L,

    # Metric pairs to fit as te(a, b) tensor smooths instead of additive
    # s(a) + s(b), when both appear together in a 2-metric combo. Pairs not
    # listed here default to additive. Names are the z-scored column names.
    tensor_pairs = list(
        c("pre_pow_slow_alpha_z", "pre_pow_fast_alpha_z")
    ),

    # -- Smooth Complexity Ceiling ----
    k_max = 10L,

    # -- bam() options ----
    nthreads = 4L,
    # discrete = TRUE is much faster. Tested: with very few subjects
    # (~12 subjects / 2 experiments) a discrete fit that also contains a
    # subject-level term (age_z, s(experiment_id)) returned meaningless
    # smooth p-values (F = 695, p = 0.76); at 120 subjects / 9 experiments
    # discrete and non-discrete fits agreed. Set FALSE for small pilot runs.
    bam_discrete = TRUE,

    # -- Follow-up on significant models (sex_diff_core.R, plot_gamm_hits.R) ----
    # A "hit" = a single-metric model of this feature_set with p < hits_p_sig.
    hits_p_sig = 0.05,
    # Reference level for the sex-difference refit; the difference curve is
    # (other - reference). Must match the coding in subjects.sex ("M" / "F").
    sex_reference = "F",
    sex_min_subjects = 10L, # per sex, else the refit is skipped

    # -- Held-out replication (holdout_replicate.R) ----
    # Primary criterion: does the metric improve prediction for people the
    # model has never seen?
    #   "subject": error in each held-out subject's mean rating (who is more
    #              pain-sensitive), after laser power, age, sex, cap size and
    #              experiment
    #   "trial"  : error in individual trials (trial-to-trial pain)
    # Fix this before looking at the results.
    holdout_primary      = "subject",
    holdout_alpha        = 0.05,   # BH-adjusted, across all hits
    holdout_min_pct      = 1.0,    # and at least this % lower error than baseline
    holdout_n_boot       = 2000L,  # subject-resampling bootstrap draws
    holdout_seed         = 42L,
    holdout_min_subjects = 5L,

    # -- Figures (plot_gamm_hits.R) ----
    fig_dir    = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/gamm_outputs/figures",
    fig_width  = 6.5,   # inches; sex figures are 2 panels at 1.8 x this width
    fig_height = 4.2,
    fig_dpi    = 300,
    # Validated categorical pair (colour-blind safe, >= 3:1 on white)
    col_main   = "#2a78d6",
    col_sex    = c(F = "#eb6834", M = "#2a78d6"),
    col_diff   = "#52514e",

    # -- Multiple Comparisons Stance ----
    # No correction applied at the screening stage (see source_core.R header
    # comment) - replication in the held-out CV split is the de facto
    # correction. This flag exists so that decision is visible in config,
    # not buried in code, in case it's revisited later.
    apply_fdr_correction = FALSE
)