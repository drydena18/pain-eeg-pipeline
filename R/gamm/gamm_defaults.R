# =============================================================================
# gamm_defaults.R
# =============================================================================
# Configuration for the GAMM stage (source_core.R + channel_core.R) and the
# scripts that read its results (sex_diff_core.R, plot_gamm_hits.R,
# holdout_replicate.R, compare_levels.R, classical_tests.R,
# ../validate/power_analysis.R, ../robustness/robustness_core.R).
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
#
# V3 (Oct 2026):
#   - bam_discrete now FALSE. With discrete = TRUE, mgcv (tested 1.9-1)
#     returns invalid p-values for subject-level smooths: the between-subject
#     term of a true trait effect scored p = 0.999 (non-discrete: p = 7e-7;
#     same edf, same F), and null between terms also scored p ~ 1, at 66,
#     135 and 360 subjects. The between-subject p-value defines a hit, so the
#     screen produced no valid hits. Estimates and predictions are unaffected;
#     only the smooth-term tests are. gamm_fitted now records bam_discrete.
#   - robustness, comparison and classical-test settings added (bottom).
# =============================================================================

gamm_cfg <- list(

    # -- Database Paths ----
    data_db_path    = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/data.duckdb",
    results_db_path = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/results.duckdb",

    # -- Output (model .rds objects still written to disk, per ROI/channel) ----
    # between_within models: <out_dir>/<level>/<unit>/<model_name>.rds
    # pooled models:         <out_dir>/<level>/<unit>/pooled/<model_name>.rds
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
    # gamm_fitted and the follow-up tables are keyed on the decomposition,
    # so both can be stored side by side; every follow-up script reads the
    # decomposition set here.
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
    # Keep FALSE (see V3 note in the header). discrete = TRUE is several
    # times faster but its smooth-term p-values are invalid for terms that
    # are constant within a subject (the _between_z part, the subject random
    # effect), at any N tested. Only consider TRUE with
    # metric_decomposition = "pooled" and primary_part irrelevant, and check
    # first with check_discrete_between_p.R on your mgcv version.
    bam_discrete = FALSE,

    # -- Follow-up on significant models (sex_diff_core.R, plot_gamm_hits.R) ----
    # A "hit" = a single-metric model of this feature_set and decomposition
    # with p < hits_p_sig.
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

    # -- Robustness (../robustness/robustness_core.R) ----
    # Every test refits the screen model and its nested baseline (metric term
    # removed) and scores population-level predictions for people the refit
    # never saw, exactly as holdout_replicate.R does.
    #   loeo         leave one EXPERIMENT out: train on the other experiments'
    #                training subjects, predict every subject of the left-out
    #                experiment (new hardware / paradigm / lab)
    #   kfold        grouped k-fold CV inside the TRAINING set (whole
    #                subjects per fold, folds stratified by cap size)
    #   permutation  screen model held fixed; the metric's primary part is
    #                shuffled at SUBJECT level (within experiment) on the
    #                held-out set -> null for the held-out improvement
    #   split_half   refit on random halves of the training subjects
    #                (stratified by cap size); correlation of the two fitted
    #                metric curves = stability of the functional form
    #   fooof_sensitivity  refit after dropping subjects whose subject-level
    #                fooof_r2 at that unit is below each threshold
    robustness = list(
        tests            = c("loeo", "kfold", "permutation", "split_half", "fooof_sensitivity"),
        levels           = c("source", "channel"),
        # TRUE: only hits that replicated in holdout_replication (run
        # holdout_replicate.R first). Each hit costs ~2 refits per fold /
        # half / threshold, so running every screening hit is slow at full N.
        only_replicated  = TRUE,
        kfold_k          = 10L,
        n_perm           = 1000L,
        perm_within_experiment = TRUE,
        n_split_half     = 20L,
        fooof_thresholds = c(0.75, 0.85),
        min_test_subjects = 5L,
        n_boot           = 2000L,
        seed             = 42L
    ),

    # -- Channel vs source comparison (compare_levels.R) ----
    compare = list(
        # Paired held-out comparison: for each metric with a hit at both
        # levels, the best ROI and the best channel (by held-out subject
        # error) predict the same held-out subjects. Best-of-64 channels vs
        # best-of-18 ROIs favours channels; the report says so.
        n_boot = 2000L,
        grid_n = 100L
    ),

    # -- Classical tests (classical_tests.R) ----
    classical = list(
        split          = "all",        # "all" | "train" | "holdout"
        rois           = c("S1", "S2", "ACC", "Insula", "dlPFC", "M1"),   # bilateral ROIs
        rmanova_metrics = c("pre_sf_balance", "pre_paf_cog_hz", "erd_asym"),
        phase_metrics  = c("slow_phase", "slow_phase_post"),
        alpha          = 0.05,
        out_dir        = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/gamm_outputs/classical_tests"
    ),

    # -- Multiple Comparisons Stance ----
    # No correction applied at the screening stage (see source_core.R header
    # comment) - replication in the held-out CV split is the de facto
    # correction. This flag exists so that decision is visible in config,
    # not buried in code, in case it's revisited later.
    apply_fdr_correction = FALSE
)