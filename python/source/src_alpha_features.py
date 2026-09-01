"""
src_alpha_features.py - Shared per-window alpha feature computation.

Called once per window (whole / pre / post) from source_core.py, replacing
the previously-duplicated BI/LR/CoG computation that lived separately inside
src_prestim.py and src_poststim.py.
"""

from __future__ import annotations

import numpy as np

from src_spectral import src_psd_welch, src_bandpower

_EPS0 = 1e-12


# ===================================================================
# SINGLE-PSSD FEATURE COMPUTATION
# ===================================================================
def src_compute_alpha_features(
        freqs: np.ndarray,
        psd: np.ndarray,
        alpha: tuple[float, float],
        slow: tuple[float, float],
        fast: tuple[float, float],
) -> dict:
    """
    Compute the 10 alpha metrics + psi_cog from a single PSD

    Args:
        freqs   : (n_freqs,) frequency axis, Hz
        psd     : (n_freqs,) power spectral density
        alpha, slow, fast : (lo, hi) band bounds in Hz
    """
    pow_slow_alpha = src_bandpower(freqs, psd, slow[0], slow[1])
    pow_fast_alpha = src_bandpower(freqs, psd, fast[0], fast[1])
    pow_alpha_total = src_bandpower(freqs, psd, alpha[0], alpha[1])

    idxA = (freqs >= alpha[0]) & (freqs <= alpha[1])
    if np.any(idxA):
        nums = float(np.trapezoid(psd[idxA] * freqs[idxA], freqs[idxA]))
    else:
        nums = 0.0
    den = (pow_alpha_total if not np.isnan(pow_alpha_total) else 0.0) + _EPS0
    paf_cog_hz = nums / den

    ps = pow_slow_alpha if not np.isnan(pow_slow_alpha) else 0.0
    pf = pow_fast_alpha if not np.isnan(pow_fast_alpha) else 0.0
    pa = pow_alpha_total if not np.isnan(pow_alpha_total) else 0.0

    sf_ratio = ps / (pf + _EPS0)
    sf_logratio = float(np.log(ps + _EPS0) - np.log(pf + _EPS0))
    sf_balance = (ps - pf) / (ps + pf + _EPS0)
    slow_alpha_frac = ps / (ps + pf + _EPS0)
    rel_slow_alpha = ps / (pa + _EPS0)
    rel_fast_alpha = pf / (pa + _EPS0)

    psi_cog = sf_balance * (paf_cog_hz - 10.0)

    return {
        "pow_slow_alpha": pow_slow_alpha,
        "pow_fast_alpha": pow_fast_alpha,
        "pow_alpha_total": pow_alpha_total,
        "paf_cog_hz": paf_cog_hz,
        "sf_ratio": sf_ratio,
        "sf_logratio": sf_logratio,
        "sf_balance": sf_balance,
        "slow_alpha_frac": slow_alpha_frac,
        "rel_slow_alpha": rel_slow_alpha,
        "rel_fast_alpha": rel_fast_alpha,
        "psi_cog": psi_cog,
    }
                                   

# ===================================================================
# PER-TRIAL x ROI, FOR ONE WINDOW (whole / pre / post - caller decides
# which tc_window array to pass, e.g. tc, tc_pre, or tc_post)
# ===================================================================
def src_compute_window_alpha_features(
        tc_window: np.ndarray,
        sfreq: float,
        alpha: tuple[float, float],
        slow: tuple[float, float],
        fast: tuple[float, float],
        fmin: float,
        fmax: float,
        psd_window_sec: float = 0.5,
) -> list[dict]:
    """
    Compute the 11-metric alpha feature set for every (trial, ROI) in a 
    single window's time course.

    Args:
        tc_window       : (n_epochs, n_rois, n_times_window)
        sfreq           : Sampling frequency in Hz
        alpha, slow, fast : (lo, hi) band bounds in Hz
        fmin, fmax      : PSD frequency range for Welch estimation
        psd_window_sec  : Welch segment length (s)

    Returns:
        List of dicts, one per (trial, ROI): {trial, roi_idx, <11 metrics>}
        (unprefixed - caller prefixes with whole_/pre_/post_ as needed)
    """
    n_epochs, n_rois, _ = tc_window.shape
    rows: list[dict] = []

    for ei in range(n_epochs):
        for ri in range(n_rois):
            freqs, psd = src_psd_welch(tc_window[ei, ri, :], sfreq, fmin, fmax, psd_window_sec)
            feat = src_compute_alpha_features(freqs, psd, alpha, slow, fast)
            feat["trial"] = ei + 1
            feat["roi_idx"] = ri
            rows.append(feat)

    return rows


# ===================================================================
# GRAND-AVERAGE, FOR ONE WINDOW
# ===================================================================
def src_compute_ga_window_alpha_features(
        tc_window: np.ndarray,
        sfreq: float,
        alpha: tuple[float, float],
        slow: tuple[float, float],
        fast: tuple[float, float],
        fmin: float,
        fmax: float,
        psd_window_sec: float = 0.5,
) -> tuple[list[dict], dict]:
    """
    Compute GA alpha features for one window, per ROI.

    Returns:
        ga_rows     : List of dicts (one per ROI): {roi_idx, <11 metrics>}
        psd_by_roi  : Dict mapping roi_idx -> (freqs, psd), for reuse by
                      FOOOF / plotting so the PSD isn't recomputed a second time.
    """
    n_epochs, n_rois, _ = tc_window.shape
    ga_rows: list[dict] = []
    psd_by_roi: dict = {}

    for ri in range(n_rois):
        x_ga = np.mean(tc_window[:, ri, :], axis = 0)
        freqs, psd = src_psd_welch(x_ga, sfreq, fmin, fmax, psd_window_sec)
        psd_by_roi[ri] = (freqs, psd)

        feat = src_compute_alpha_features(freqs, psd, alpha, slow, fast)
        feat["roi_idx"] = ri
        ga_rows.append(feat)

    return ga_rows, psd_by_roi