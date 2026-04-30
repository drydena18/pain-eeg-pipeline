"""
src_lep.py - Laser-evoked potential (LEP) feature extraction

LEP's are extracted from the source-space ROI time courses rather than scalp
channels. This gives anatomically labelled responses (e.g., S1 N2, ACC P2)
that are directly interpretable in the source domain.

Features extracted (per trial x ROI)
-------------------------------------
    n2_amp_uv       N2 peak amplitude (most negative value in N2 window, µV)
    n2_lat_ms       N2 peak latency in ms
    p2_amp_uv       P2 peak amplitude (most positive value in P2 window, µV)
    p2_lat_ms       P2 peak latency in ms
    n2p2_amp_uv     N2-P2 peak-to-peak amplitude (P2 - N2, always positive)
    n2_mean_uv      Mean amplitude in the N2 window (area-based measure)
    p2_mean_uv      Mean amplitude in the P2 window

Grand-average features (per ROI)
---------------------------------
Same metrics applied to the mean time course across trials.

Notes
------
    - LEP windows (N2, P2) are configurable via the JSON config to accomodate 
      different ROIs and paradigm-specific latency shifts
    - All amplitudes are in the native units of the source estimate (typically 
      nA*m or arbitrary sLORETA units). Label them relative to each other; do not
      compare raw amplitudes across ROIs or studies)
    - For per-trial extraction, the same peak-search logic is applied to the
      single-trial time course without smoothing, which is intentionally noisy,
      subsequent GAMM modelling uses the amplitude as a continuous predictor rather
      than trusting the peak itself
"""

from __future__ import annotations

import numpy as np

# ====================================================================
# PEAK SEARCH PRIMITIVES
# ====================================================================
def _find_peak(
        x: np.ndarray,
        times: np.ndarray,
        t_lo: float,
        t_hi: float,
        polarity: str,
) -> tuple[float, float]:
    """
    Find the peak amplitude and latency within a time window

    Args:
        x           : 1-D time course array
        times       : Corresponding time axis in seconds
        t_lo        : Window start in seconds
        t_hi        : Window end in seconds
        polarity    : "neg" to find minimum (N2), "pos" to find maximum (P2)

    Returns:
        (peak_amp, peak_lat_ms) - both NaN if the window is empty
    """
    mask = (times >= t_lo) & (times <= t_hi)
    if not np.any(mask):
        return float("nan"), float("nan")
    
    x_win = x[mask]
    t_win = times[mask]

    if polarity == "neg":
        idx = int(np.argmin(x_win))
    else:
        idx = int(np.argmax(x_win))

    return float(x_win[idx]), float(t_win[idx] * 1000.00) # amp, latency in ms

# ====================================================================
# PER-TRIAL LEP EXTRACTION
# ====================================================================
def src_compute_lep_trial(
        tc_post: np.ndarray,
        times_post: np.ndarray,
        n2_window: tuple[float, float],
        p2_window: tuple[float, float],
) -> list[dict]:
    """
    Extract per-trial LEP features for every (trial, ROI).

    Args:
        tc_post     : (n_epochs, n_rois, n_times_post) - post-stim time courses
        times_post  : (n_times_post,) in seconds (should start at ~0 or post-offset)
        n2_window   : (t_lo, t_hi) in seconds for the N2 search window
        p2_window   : (t_lo, t_hi) in seconds for the P2 search window

    Returns:
        List of dicts, one per (trial, ROI):
            trial, roi_idx,
            n2_amp, n2_lat_ms,
            p2_amp, p2_lat_ms,
            n2p2_amp,
            n2_mean, p2_mean
    """
    n_epochs, n_rois, _ = tc_post.shape
    rows: list[dict] = []

    for ei in range(n_epochs):
        for ri in range(n_rois):
            x = tc_post[ei, ri, :]

            n2_amp, n2_lat = _find_peak(x, times_post, n2_window[0], n2_window[1], "neg")
            p2_amp, p2_lat = _find_peak(x, times_post, p2_window[0], p2_window[1], "pos")

            n2p2 = (p2_amp - n2_amp) if not (np.isnan(p2_amp) or np.isnan(n2_amp)) \
                else float("nan")
            
            # Mean amplitude (area-under-the-curve)
            n2_mask = (times_post >= n2_window[0]) & (times_post <= n2_window[1])
            p2_mask = (times_post >= p2_window[0]) & (times_post <= p2_window[1])
            n2_mean = np.mean(x[n2_mask]) if np.any(n2_mask) else float("nan")
            p2_mean = np.mean(x[p2_mask]) if np.any(p2_mask) else float("nan")

            rows.append({
                "trial": ei + 1,
                "roi_idx": ri,
                "n2_amp": n2_amp,
                "n2_lat_ms": n2_lat,
                "p2_amp": p2_amp,
                "p2_lat_ms": p2_lat,
                "n2p2_amp": n2p2,
                "n2_mean": n2_mean,
                "p2_mean": p2_mean
            })

    return rows

# ====================================================================
# GRAND-AVERAGE LEP EXTRACTION
# ====================================================================
def src_compute_leg_ga(
        tc_post: np.ndarray,
        times_post: np.ndarray,
        n2_window: tuple[float, float],
        p2_window: tuple[float, float],
) -> list[dict]:
    """
    Extract LEP features from the grand-average (mean across trials) time course
    for each ROI.

    Returns:
        List of dicts, one per ROI, with the same fields as src_compute_lep_trial
        but without the 'trial' key, and with an added 'n_trials' field
    """
    n_epochs, n_rois, _ = tc_post.shape
    rows: list[dict] = []

    for ri in range(n_rois):
        x_ga = np.mean(tc_post[:, ri, :], axis = 0)

        n2_amp, n2_lat = _find_peak(x_ga, times_post, n2_window[0], n2_window[1], "neg")
        p2_amp, p2_lat = _find_peak(x_ga, times_post, p2_window[0], p2_window[1], "pos")

        n2p2 = (p2_amp - n2_amp) if not (np.isnan(p2_amp) or np.isnan(n2_amp)) \
            else float("nan")
        
        n2_mask = (times_post >= n2_window[0]) & (times_post <= n2_window[1])
        p2_mask = (times_post >= p2_window[0]) & (times_post <= p2_window[1])
        
        rows.append({
            "roi_idx": ri,
            "n_trials": n_epochs,
            "n2_amp": n2_amp,
            "n2_lat_ms": n2_lat,
            "p2_amp": p2_amp,
            "p2_lat_ms": p2_lat,
            "n2p2_amp": n2p2,
            "n2_mean": np.mean(x_ga[n2_mask]) if np.any(n2_mask) else float("nan"),
            "p2_mean": np.mean(x_ga[p2_mask]) if np.any(p2_mask) else float("nan"),
        })

    return rows