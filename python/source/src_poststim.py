"""
src_poststim.py - Post-stimulus Hilbert phase + inter-trial phase coherence.

V2.0.0 changes vs V1.x:
    - REMOVED pow_slow_post/pow_fast_post/pow_alpha_post computation. These
      are now computed generically as part of the post_ window family by
      src_alpha_features.src_compute_window_alpha_features /
      src_compute_ga_window_alpha_features, called from source_core.py.
    - REMOVED ERD_slow/ERD_fast/delta_ERD computation, including the old
      NaN-on-low-power denominator guard. This is now src_erd.py
      (src_compute_erd_metrics), which always returns a numeric ERD value
      using an epsilon floor and flags low-power trials separately via
      p5_flag, matching spec_compute_interaction_metrics.m. See src_erd.py
      for why the noise-floor/percentile stats are computed per-ROI here
      rather than pooled across ROIs the way MATLAB pools across channels.
    - This file now owns ONLY:
        - Hilbert instantaneous phase at poststim_ref_t (slow_phase_post,
          sin_phase_post, cos_phase_post) — unchanged from V1.x.
        - Inter-trial phase coherence (ITC) — unchanged from V1.x.

Metrics Implemented
--------------------
Hilbert instantaneous phase (post-stimulus):
    slow_phase_post Instantaneous phase of the 8-10 Hz analytic signal
                    at the sample nearest to poststim_ref_t (default: 0.2 s,
                    corresponding to early post-stimulus alpha suppression).
    sin_phase_post  sin(slow_phase_post)  - linear GAMM regressor
    cos_phase_post  cos(slow_phase_post)  - linear GAMM regressor

Inter-trial phase coherence (ITC):
    Computed at the GA level (one value per ROI, not per trial) because ITC
    (also called phase-locking value) is inherently a population statistic.
    Written to the GA CSV.

References
──────────
    Iannetti & Mouraux (2010); Ohara et al. (2004); Furman et al. (2020).
"""

from __future__ import annotations

import numpy as np
from scipy.signal import hilbert

from src_spectral import src_bandpass_filter

_EPS = 1e-12


# =============================================================================
# PER-TRIAL POST-STIMULUS PHASE
# =============================================================================
def src_compute_poststim_phase(
    tc_full: np.ndarray,
    times_full: np.ndarray,
    sfreq: float,
    slow: tuple[float, float],
    poststim_ref_t: float = 0.2,
) -> list[dict]:
    """
    Compute the Hilbert instantaneous phase of the slow-alpha (8-10 Hz)
    signal at poststim_ref_t for every (trial, ROI).

    Args:
        tc_full         : (n_epochs, n_rois, n_times) — full epoch.
        times_full      : (n_times,) — full epoch time axis in seconds.
        sfreq           : Sampling frequency in Hz.
        slow            : (lo, hi) slow-alpha sub-band bounds in Hz.
        poststim_ref_t  : Time (s) at which the post-stim phase is read.

    Returns:
        List of dicts, one per (trial, ROI):
            {trial, roi_idx, slow_phase_post, sin_phase_post, cos_phase_post}
    """
    n_epochs, n_rois, _ = tc_full.shape
    ref_idx = int(np.argmin(np.abs(times_full - poststim_ref_t)))

    rows: list[dict] = []
    for ei in range(n_epochs):
        for ri in range(n_rois):
            try:
                x_filt = src_bandpass_filter(tc_full[ei, ri, :], sfreq, slow[0], slow[1])
                analytic = hilbert(x_filt)
                phase_post = float(np.angle(analytic[ref_idx]))
            except Exception:
                phase_post = float("nan")

            rows.append({
                "trial":           ei + 1,
                "roi_idx":         ri,
                "slow_phase_post": phase_post,
                "sin_phase_post":  float(np.sin(phase_post)) if not np.isnan(phase_post) else float("nan"),
                "cos_phase_post":  float(np.cos(phase_post)) if not np.isnan(phase_post) else float("nan"),
            })

    return rows


# =============================================================================
# GRAND-AVERAGE POST-STIMULUS PHASE
# =============================================================================
def src_compute_ga_poststim_phase(
    tc_full_ga: np.ndarray,
    times_full: np.ndarray,
    sfreq: float,
    slow: tuple[float, float],
    poststim_ref_t: float = 0.2,
) -> list[dict]:
    """
    Same as src_compute_poststim_phase but for the GA (trial-averaged) time
    course.

    Args:
        tc_full_ga : (1, n_rois, n_times) — GA mean time course, keepdims=True.
        Other args same as src_compute_poststim_phase.

    Returns:
        List of dicts, one per ROI: {roi_idx, slow_phase_post, sin_phase_post, cos_phase_post}
    """
    _, n_rois, _ = tc_full_ga.shape
    ref_idx = int(np.argmin(np.abs(times_full - poststim_ref_t)))

    rows: list[dict] = []
    for ri in range(n_rois):
        try:
            x_filt = src_bandpass_filter(tc_full_ga[0, ri, :], sfreq, slow[0], slow[1])
            analytic = hilbert(x_filt)
            phase_post = float(np.angle(analytic[ref_idx]))
        except Exception:
            phase_post = float("nan")

        rows.append({
            "roi_idx":         ri,
            "slow_phase_post": phase_post,
            "sin_phase_post":  float(np.sin(phase_post)) if not np.isnan(phase_post) else float("nan"),
            "cos_phase_post":  float(np.cos(phase_post)) if not np.isnan(phase_post) else float("nan"),
        })

    return rows


# =============================================================================
# GRAND-AVERAGE POST-STIM: INTER-TRIAL PHASE COHERENCE (ITC)
# =============================================================================
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
    post-stimulus window. Unchanged from V1.x.

    ITC (also called the phase-locking value, PLV) measures how consistently
    the slow-alpha oscillation is phase-reset by the laser stimulus across
    trials. A high ITC at a given time point means the phase is tightly
    clustered across trials; ITC = 1 is perfect coherence, ITC = 0 is uniform
    phase distribution.

        ITC(t) = | (1/K) sum_k  exp(i . phi_k(t)) |

    where phi_k(t) is the instantaneous phase on trial k at time t.

    Args:
        tc_full    : (n_epochs, n_rois, n_times) — full epoch source time courses.
        times_full : (n_times,) in seconds.
        sfreq      : Sampling frequency in Hz.
        slow       : Slow alpha band (lo, hi) in Hz.
        post_tmin  : Start of post-stim window to report ITC over.
        post_tmax  : End of post-stim window.

    Returns:
        List of dicts, one per ROI:
            roi_idx, itc_mean, itc_peak, itc_peak_latency_ms
    """
    n_epochs, n_rois, _ = tc_full.shape
    mask_post = (times_full >= post_tmin) & (times_full <= post_tmax)
    times_post = times_full[mask_post]

    rows: list[dict] = []

    for ri in range(n_rois):
        phases = np.full((n_epochs, int(np.sum(mask_post))), np.nan, dtype=complex)

        for ei in range(n_epochs):
            try:
                x_filt   = src_bandpass_filter(tc_full[ei, ri, :], sfreq, slow[0], slow[1])
                analytic = hilbert(x_filt)
                phases[ei, :] = np.exp(1j * np.angle(analytic[mask_post]))
            except Exception:
                pass

        itc = np.abs(np.nanmean(phases, axis=0))

        itc_mean = float(np.nanmean(itc))
        peak_idx = int(np.nanargmax(itc))
        itc_peak = float(itc[peak_idx])
        itc_peak_lat_ms = float(times_post[peak_idx] * 1000.0)

        rows.append({
            "roi_idx":             ri,
            "itc_mean":            itc_mean,
            "itc_peak":            itc_peak,
            "itc_peak_latency_ms": itc_peak_lat_ms,
        })

    return rows