"""
src_prestim.py - Pre-stimulus alpha metrics

All metrics are computed on the pre-stimulus source-space time course
(shape: n_epochs x n_rois x n_times_pre).

Metrics Implemented
--------------------
Power family (per trial x ROI):
    pow_slow        slow alpha [8, 10] Hz bandpower
    pow_fast        fast alpha [10, 12] Hz bandpower
    pow_alpha       total alpha [8, 12] Hz bandpower

Sub-band balance index [Metric 1 - BI_pre]
    BI_pre          (slow - fast) / (slow + fast + ε) ∈ [-1, +1]

Log ratio [Metric 2 - LR_pre] 
    LR_pre          ln(slow * ε) - ln(fast + ε) ∈ (-∞, +∞)

Centre of Gravity [Metric 3 - CoG_pre]
    CoG_pre         Spectral CoG over full alpha [8, 12] Hz (Hz)

Interaction term [ψ_cog]
    psi_cog         BI_pre x (CoG_pre - 10)

Hilbert instantaneous phase at stimulus onset [Metric 4]
    slow_phase      Phase of the slow-alpha (8-10 Hz) analytic signal at t = 0
    sin_phase       sin(slow_phase)     - linear GAMM regressor
    cos_phase       cos(slow_phase)     - linear GAMM regressor

ERD (event-related desynchronization) [Metric 5 - ΔERD]
    Computed in src_poststim.py because it requires post-stimulus window

Temporal variability index [Metric 6 - TVI_alpha]
    TVI_alpha       nMSSD of the per-trial BI_pre sequence.
                    This is a per-subject scalar, written to the GA CSV

References
-----------
    Furman et al. (2019, 2020); Tu et al. (2016); Nickel et al. (2022);
    Busch et al. (2009); Mathewson et al. (2009); Li et al. (2018)
"""

from __future__ import annotations

import numpy as np
from scipy.signal import hilbert

from src_spectral import src_psd_welch, src_bandpower, src_cog, src_bandpass_filter

# -- Noise-floor epsilon ------
# Use a small fixed constant rather than an arbitrary 1e-10
# The exact value matters less than being consistent across metrics
_EPS = 1e-12

# ====================================================================
# PER-TRIAL PRE-STIMULUS METRICS
# ====================================================================
def src_compute_prestim_metrics(
        tc_pre: np.ndarray,
        tc_full: np.ndarray,
        times_full: np.ndarray,
        sfreq: float,
        alpha: tuple[float, float],
        slow: tuple[float, float],
        fast: tuple[float, float],
        fmin: float,
        fmax: float,
        psd_window_sec: float = 0.5,
) -> list[dict]:
    """
    Compute all per-trial pre-stimulus alpha metrics for every (trial, ROI)

    The Hilbert phase at stimulus onset is extracted from the *full* epoch time
    course (tc_full) rather than the cropped pre-stimulus window (tc_pre) to
    avoid edge effects from the Hilbert transform. The pahse is then read at
    the sample nearest to t = 0 (stimulus onset).

    Args:
        tc_pre      : (n_epochs, n_rois, n_times_pre)  — pre-stim window only.
        tc_full     : (n_epochs, n_rois, n_times)       — full epoch.
        times_full  : (n_times,) — full epoch time axis in seconds.
        sfreq       : Sampling frequency in Hz.
        alpha       : (lo, hi) bounds for total alpha band (Hz).
        slow        : (lo, hi) bounds for slow alpha sub-band (Hz).
        fast        : (lo, hi) bounds for fast alpha sub-band (Hz).
        fmin, fmax  : PSD frequency range for Welch estimation.
        psd_window_sec : Welch segment length (s). Default 0.5 s is a good
                         balance between frequency resolution and stability for
                         typical pre-stim windows of 500–1000 ms.

    Returns:
        List of dicts, one per (trial, ROI). Each dict has fields:
            trial, roi, pow_slow, pow_fast, pow_alpha,
            BI_pre, LR_pre, CoG_pre, psi_cog,
            slow_phase, sin_phase, cos_phase
    """
    n_epochs, n_rois, _ = tc_pre.shape

    # Index of t = 0 in the full epoch (stimulus onset)
    t0_idx = int(mp.argmin(np.abs(times_full)))

    rows: list[dict] = []

    for ei in range(n_epochs):
        for ri in range(n_rois):
            x_pre = tc_pre[ei, ri, :]
            x_full = tc_full[ei, ri, :]

            # -- Power metrics -----
            freqs, psd = src_psd_welch(x_pre, sfreq, fmin, fmax, psd_window_sec)
            pow_slow = src_bandpower(freqs, psd, slow[0], slow[1])
            pow_fast = src_bandpower(freqs, psd, fast[0], fast[1])
            pow_alpha = src_bandpower(freqs, psd, alpha[0], alpha[1])

            # -- BI_pre -----
            bi_pre = (pow_slow - pow_fast) / (pow_slow + pow_fast + _EPS)

            # -- LR_pre -----
            lr_pre = float(np.log(pow_slow + _EPS) - np.log(pow_fast + _EPS))

            # -- CoG_pre -----
            cog_pre = src_cog(freqs, psd, alpha[0], alpha[1])

            # -- ψ_cog -----
            # Centred around 10 Hz (slow/fast boundary) so the main effects of
            # BI_pre and CoG_pre remain interpretable at the mean
            psi_cog = bi_pre * (cog_pre - 10.0) if not np.isnan(cog_pre) else float("nan")

            # -- Slow-alpha Hilbert phase at t = 0 -----
            # Filter the full epoch (avoids end-of-window edge effects),
            # then extract the instantaneous phase at the stimulus onset sample.
            try:
                x_filt = src_bandpass_filter(x_full, sfreq, slow[0], slow[1])
                analytic = hilbert(x_filt)
                slow_phase = float(np.angle(analytic(t0_idx)))
            except Exception:
                slow_phase = float("nan")

            rows.append({
                "trial": ei + 1,
                "roi_idx": ri,
                "pow_slow": pow_slow,
                "pow_fast": pow_fast,
                "pow_alpha": pow_alpha,
                "BI_pre": bi_pre,
                "LR_pre": lr_pre,
                "CoG_pre": cog_pre,
                "psi_cog": psi_cog,
                "slow_phase": slow_phase,
                "sin_phase": float(np.sin(slow_phase)) if not np.isnan(slow_phase) else float("nan"),
                "cos_phase": float(np.cos(slow_phase)) if not np.isnan(slow_phase) else float("nan"),
            })

    return rows

# ====================================================================
# GRAND-AVERAGE PRE-STIMULUS METRICS
# ====================================================================
def src_compute_ga_prestim_metrics(
        tc_pre: np.ndarray,
        sfreq: float,
        alpha: tuple[float, float],
        slow: tuple[float, float],
        fast: tuple[float, float],
        fmin: float,
        fmax: float,
        psd_window_sec: float = 0.5,
) -> tuple[list[dict], dict]:
    """
    Compute GA pre-stimulus metrics for each ROI

    The GA time course per ROI is the mean across trials (in the time domain),
    following the same convention as spectral_core.m

    Returns:
        ga_rows : List of dicts (one per ROI) for the GA CSV
        psd_by_roi : Dict mapping roi_idx -> (freqs, psd) for FOOOF and plotting
    """
    n_epochs, n_rois, _ = tc_pre.shape
    ga_rows: list[dict] = []
    psd_by_roi: dict = {}

    for ri in range(n_rois):
        x_ga = np.mean(tc_pre[:, ri, :], axis = 0)
        freqs, psd = src_psd_welch(x_ga, sfreq, fmin, fmax, psd_window_sec)
        psd_by_roi[ri] = (freqs, psd)

        pow_slow = src_bandpower(freqs, psd, slow[0], slow[1])
        pow_fast = src_bandpower(freqs, psd, fast[0], fast[1])
        pow_alpha = src_bandpower(freqs, psd, alpha[0], alpha[1])
        bi_pre = (pow_slow - pow_fast) / (pow_slow + pow_fast + _EPS)
        lr_pre = float(np.log(pow_slow + _EPS) - np.log(pow_fast + _EPS))
        cog_pre = src_cog(freqs, psd, alpha[0], alpha[1])
        psi_cog = bi_pre * (cog_pre - 10.0) if not np.isnan(cog_pre) else float("nan")

        ga_rows.append({
            "roi_idx": ri,
            "pow_slow": pow_slow,
            "pow_fast": pow_fast,
            "pow_alpha": pow_alpha,
            "BI_pre": bi_pre,
            "LR_pre": lr_pre,
            "CoG_pre": cog_pre,
            "psi_cog": psi_cog,
        })

    return ga_rows, psd_by_roi

# ====================================================================
# TEMPORAL VARIABILITY INDEX (TVI_alpha / nMSSD)
# ====================================================================
def src_compute_tvi_alpha(bi_pre_sequence: np.ndarray) -> float:
    """
    Compute the temporal variability index (TVI_alpha) from the per-trial
    BI_pre sequence for a single subject x ROI

    TVI_alpha is the normalized mean square successive difference (nMMSD):

        MSSD = mean( (b[k+1] - b[k])² ) for k = 1 ... K - 1
        Var = variance( b )
        TVI_alpha = MSSD / Var

    This isolates the *temporal autocorrelation structure* of the pre-stimulus
    state rather than its total amplitude variability
    
    Range: [0, 2]. Near 0 = rigid (slowly varying); near 2 = maximally
    alternating (each trial flips sign relative to the previous)

    Args:
        bi_pre_sequence : 1-D array of pre-trial BI_pre values, shape (K,)
    
    Returns:
        TVI_alpha as a float, or NaN if K < 3 (too few trials to be meaningful)
    """
    b = np.asarray(bi_pre_sequence, dtype = float)
    b = b[~np.isnan(b)]
    K = len(b)
    if K < 3:
        return float("nan")
    
    mssd = float(np.mean(np.diff(b) ** 2))
    var = float(np.var(b, ddof = 0))

    if var < _EPS:
        return float("nan")
    
    return mssd / var