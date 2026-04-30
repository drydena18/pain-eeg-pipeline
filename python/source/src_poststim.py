"""
src_poststim.py - Post-stimulus alpha metrics

Metrics Implemented
--------------------
ERD (event-related desynchronization) family [Metric 5 - sub-band ERD asymmetry]
    ERD_slow        Fractional slow-alpha power change: (post - pre) / pre
    ERD_fast        Fractional fast-alpha power change: (post - pre) / pre
    delta_ERD       ERD_slow - ERD_fast
                    < 0 : fast-ERD dominates (nociceptive gating)
                    > 0 : slow-ERD dominates (diffuse arousal / saliency)

Post-stimulus total alpha power
    pow_slow_post   Slow alpha [8, 10] Hz power i post-stim windows
    pow_fast_post   Fast alpha [10, 12] Hz power in post-stim windows
    pow_alpha_post  Total alpha [8, 12] Hz power in post-stim windows

Hilbert instantaneous phase (post-stimulus)
    slow_phase_post Instantaneous phase of the 8-10 Hz analytic signal
                    at the sample nearest to poststim_ref_t (default: 0.2 s,
                    corresponding to early post-stimulus alpha suppression)
    sin_phase_post  sin(slow_phase_post) - linear GAMM regressor
    cos_phase_post  cos(slow_phase_post) - linear GAMM regressor

Inter-trail phase coherence (ITC):
    Computed at the GA level (one value per ROI, not per trial) because ITC
    (also called phase-locking value) is inherently a population statistic.
    Written to the GA CSV.

Notes:
-------
ERD denominators: if pre-stimulus power n either band falls below the
5th percentile of the session distribution on any trial, ERD for that trial
is set to NaN rather than propagating a near-zero denominator. The threshold
is computed once per subject x ROI and passed in as pre_percetile_threshold

References
-----------
    Iannetti & Mouraux (2010); Ohara et al.. (2004); Furman et al. (2020)
"""

from __future__ import annotations

import numpy as np
from scipy.signal import hilbert

from src_spectral import src_psd_welch, src_bandpower, src_bandpass_filter

_EPS = 1e-12

# ====================================================================
# PER-TRIAL POST-STIMULUS METRICS
# ====================================================================
def src_compute_poststim_metrics(
        tc_pre: np.ndarray,
        tc_post: np.ndarray,
        tc_full: np.ndarray,
        times_full: np.ndarray,
        sfreq: float,
        alpha: tuple[float, float],
        slow: tuple[float, float],
        fast: tuple[float, float],
        fmin: float,
        fmax: float,
        poststim_ref_t: float = 0.2,
        psd_window_sec: float = 0.5,
) -> list[dict]:
    """
    Compute per-trial post-stimulus metrics for every (trial, ROI)

    Pre-stimulus power is re-computed here from tc_pre rather than imported
    from src_prestim.py to keep the two modules independent. The values are
    identical since both use src_psd_welch with the same parameters

    Args:
        tc_pre          : (n_epochs, n_rois, n_times_pre)  — pre-stim window only.
        tc_post         : (n_epochs, n_rois, n_times_post) — post-stim window only.
        tc_full         : (n_epochs, n_rois, n_times)      — full epoch.
        times_full      : (n_times,)                       — time vector for tc_full
        sfreq           : Sampling frequency in Hz
        alpha           : (fmin_alpha, fmax_alpha)         — total alpha band
        slow            : (fmin_slow, fmax_slow)           — slow alpha band
        fast            : (fmin_fast, fmax_fast)           — fast alpha band
        fmin, fmax      : Frequency range for PSD estimation
        poststim_ref_t  : Time at which post-stim phase is sampled (default: 0.2 s)
        psd_window_sec  : Desired Welch segment length in secods (default 0.5 s)

    Returns:
        List of dicts, one per (trial, ROI), with fields:
            trial, roi_idx,
            pow_slow_post, pow_fast_post, pow_alpha_post,
            ERD_slow, ERD-fast, delta_ERD,
            slow_phase_post, sin_phase_post, cos_phase_post
    """
    n_epochs, n_rois, _ = tc_pre.shape

    # Index of poststim_ref_t in the full epoch
    ref_idx = int(np.argmin(np.abs(times_full - poststim_ref_t)))

    # -- Session-level 5th-percentile thresholds for ERD denominator guard -----
    # Compute per ROI across all epochs and both sub-bands
    slow_pre_all = np.full((n_epochs, n_rois), np.nan)
    fast_pre_all = np.full((n_epochs, n_rois), np.nan)

    for ei in range(n_epochs):
        for ri in range(n_rois):
            f, p = src_psd_welch(tc_pre[ei, ri, :], sfreq, fmin, fmax, psd_window_sec)
            slow_pre_all[ei, ri] = src_bandpower(f, p, slow[0], slow[1])
            fast_pre_all[ei, ri] = src_bandpower(f, p, fast[0], fast[1])

    slow_thresh = np.nanpercentile(slow_pre_all, 5, axis = 0)
    fast_thresh = np.nanpercentile(fast_pre_all, 5, axis = 0)

    rows: list[dict] = []

    for ei in range(n_epochs):
        for ri in range(n_rois):
            # -- Pre-stim power for ERD denominator -----
            pow_slow_pre = slow_pre_all[ei, ri]
            pow_fast_pre = fast_pre_all[ei, ri]

            # -- Post-stim power -----
            f_post, p_post = src_psd_welch(tc_post[ei, ri, :], sfreq, fmin, fmax, psd_window_sec)
            pow_slow_post = src_bandpower(f_post, p_post, slow[0], slow[1])
            pow_fast_post = src_bandpower(f_post, p_post, fast[0], fast[1])
            pow_alpha_post = src_bandpower(f_post, p_post, alpha[0], alpha[1])

            # -- ERD (fractional change; NaN if denominator near noise floor) -----
            if pow_slow_pre > slow_thresh[ri] and pow_slow_pre > _EPS:
                erd_slow = float((pow_slow_post - pow_slow_pre) / pow_slow_pre)
            else:
                erd_slow = float("nan")

            if pow_fast_pre > fast_thresh[ri] and pow_fast_pre > _EPS:
                erd_fast = float((pow_fast_post - pow_fast_pre) / pow_fast_pre)
            else:
                erd_fast = float("nan")

            if not (np.isnan(erd_slow) or np.isnan(erd_fast)):
                delta_erd = float(erd_slow - erd_fast)
            else:
                delta_erd = float("nan")

            # -- Post-stimulus Hilbert phase at poststim_ref_t -----
            try:
                x_filt = src_bandpass_filter(tc_full[ei, ri, :], sfreq, slow[0], slow[1])
                analytic = hilbert(x_filt)
                phase_post = float(np.angle(analytic[ref_idx]))
            except Exception:
                phase_post = float("nan")

            rows.append({
                "trial": ei + 1,
                "roi_idx": ri,
                "pow_slow_post": pow_slow_post,
                "pow_fast_post": pow_fast_post,
                "pow_alpha_post": pow_alpha_post,
                "ERD_slow": erd_slow,
                "ERD_fast": erd_fast,
                "delta_ERD": delta_erd,
                "slow_phase_post": phase_post,
                "sin_phase_post": float(np.sin(phase_post)) if not np.isnan(phase_post) else float("nan"),
                "cos_phase_post": float(np.cos(phase_post)) if not np.isnan(phase_post) else float("nan"),
            })

    return rows

# ====================================================================
# GRAND-AVERAGE POST-STIM: INTER-TRIAL PHASE COHERENCE (ITC)
# ====================================================================
def src_compute_itc(
        tc_full: np.ndarray,
        times_full: np.ndarray,
        sfreq: float,
        slow: tuple[float, float],
        post_tmin: float,
        post_tmax: float,
) -> list[dict]:
    """
    Compute inter-trial phase coherence (ITC) of the slow-alpha band over the
    post-stimulus window.

    ITC (also called the phase-locking value, PLV) measures how consistently
    the slow-alpha oscilaltion is phase-reset by the laser stimulus across
    trials. A high ITC at a given time point means the phase is tightly
    clustered across trials; ITC = 1 is perfect coherence, ITC = 0 is uniform 
    phase distribution.

        ITC(t) = | (1 / K) Σ_k exp(i * φ_k(t)) |
    
    where φ_k(t) is the instantaneous phase on trial k at time t.

    Args:
        tc_full         : (n_epochs, n_rois, n_times) - full epoch source time course
        times_full      : (n_times,) in seconds
        sfreq           : Sampling frequency in Hz
        slow            : Slow alpha band (lo, hi) in Hz 
        post_tmin       : Start of post-stim window to report ITC over
        post_tmax       : End of post-stim window

    Returns:
        List of dicts, one per ROI:
            roi_idx, itc_mean, itc_peak, itc_peak_latency_ms
    """
    n_epochs, n_rois, _ = tc_full.shape
    mask_post = (times_full >= post_tmin) & (times_full <= post_tmax)
    times_post = times_full[mask_post]

    rows: list[dict] = []

    for ri in range(n_rois):
        # (n_epochs, n_times_post) complex analytic signal
        phases = np.full((n_epochs, int(np.sum(mask_post))), np.nan, dtype = complex)

        for ei in range(n_epochs):
            try:
                x_filt = src_bandpass_filter(tc_full[ei, ri, :], sfreq, slow[0], slow[1])
                analytic = hilbert(x_filt)
                # Unit-norm complex phasors for ITC
                phases[ei, :] = np.exp(1j * np.angle(analytic[mask_post]))
            except Exception:
                pass # leave as NaN; valid trials still contribute

        # ITC = | mean phasor | across trials
        itc = np.abs(np.nanmean(phases, axis = 0))

        itc_mean = float(np.nanmean(itc))
        peak_idx = int(np.nanargmax(itc))
        itc_peak = float(itc[peak_idx])
        itc_peak_latency_ms = float(times_post[peak_idx] * 1000.00)

        rows.append({
            "roi_idx": ri,
            "itc_mean": itc_mean,
            "itc_peak": itc_peak,
            "itc_peak_latency_ms": itc_peak_latency_ms,
        })

    return rows