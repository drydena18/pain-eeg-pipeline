"""
src_poststim.py - Post-stimulus Hilbert phase + inter-trial phase coherence.
V 3.0.0

MATLAB twin: spec_slow_phase_features.m.

Metrics
-------
    slow_phase_post   Slow-alpha analytic phase at the sample nearest
                      phase_ref_t (default 0.2 s)
    sin_phase_post    sin(slow_phase_post)
    cos_phase_post    cos(slow_phase_post)
    itc_mean, itc_peak, itc_peak_latency_ms
                      Inter-trial phase coherence |mean_k exp(i phi_k(t))|
                      over the post-stimulus window (subject x ROI, GA CSV)

V3.0.0 changes vs V2.0.0:
    - Shared FIR band-pass (same as MATLAB), vectorised over trials x ROIs.
    - Post window selected with src_window_mask (quarter-sample tolerance).
"""

from __future__ import annotations

import numpy as np

from src_prestim import _phase_at, _phase_rows
from src_spectral import src_analytic_signal, src_bandpass_filter, src_window_mask


# =============================================================================
# PER-TRIAL POST-STIMULUS PHASE
# =============================================================================
def src_compute_poststim_phase(
    tc_full: np.ndarray,
    times_full: np.ndarray,
    sfreq: float,
    slow: tuple[float, float],
    poststim_ref_t: float = 0.2,
    trans_bw: float = 1.5,
) -> list[dict]:
    """One row per (trial, ROI): {trial, roi_idx, slow_phase_post, sin_phase_post, cos_phase_post}."""
    phase = _phase_at(tc_full, times_full, sfreq, slow, poststim_ref_t, trans_bw)
    return _phase_rows(phase, "_post", with_trial = True)


# =============================================================================
# GRAND-AVERAGE POST-STIMULUS PHASE
# =============================================================================
def src_compute_ga_poststim_phase(
    tc_full_ga: np.ndarray,
    times_full: np.ndarray,
    sfreq: float,
    slow: tuple[float, float],
    poststim_ref_t: float = 0.2,
    trans_bw: float = 1.5,
) -> list[dict]:
    """tc_full_ga: (1, n_rois, n_times). One row per ROI."""
    phase = _phase_at(tc_full_ga, times_full, sfreq, slow, poststim_ref_t, trans_bw)[0]
    return _phase_rows(phase, "_post", with_trial = False)


# =============================================================================
# INTER-TRIAL PHASE COHERENCE (subject x ROI)
# =============================================================================
def src_compute_itc(
    tc_full: np.ndarray,
    times_full: np.ndarray,
    sfreq: float,
    slow: tuple[float, float],
    post_tmin: float,
    post_tmax: float,
    trans_bw: float = 1.5,
) -> list[dict]:
    """
    ITC(t) = |(1/K) sum_k exp(i phi_k(t))| over the post-stimulus window.

    Returns one dict per ROI: roi_idx, itc_mean, itc_peak, itc_peak_latency_ms
    (first maximum).
    """
    mask = src_window_mask(times_full, post_tmin, post_tmax, sfreq)
    times_post = np.asarray(times_full)[mask]
    z = src_analytic_signal(src_bandpass_filter(tc_full, sfreq, slow[0], slow[1], trans_bw))
    unit = np.exp(1j * np.angle(z[..., mask]))           # (n_epochs, n_rois, n_post)
    itc = np.abs(np.mean(unit, axis = 0))                # (n_rois, n_post)

    rows: list[dict] = []
    for ri in range(itc.shape[0]):
        peak_idx = int(np.argmax(itc[ri]))
        rows.append({
            "roi_idx":             ri,
            "itc_mean":            float(np.mean(itc[ri])),
            "itc_peak":            float(itc[ri, peak_idx]),
            "itc_peak_latency_ms": float(times_post[peak_idx] * 1000.0),
        })
    return rows
