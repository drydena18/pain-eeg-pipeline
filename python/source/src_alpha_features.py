"""
src_alpha_features.py - Shared per-window alpha feature computation.
V 3.1.0

Called once per window (whole / pre / post) from source_core.py. The
caller prefixes the unprefixed output with whole_ / pre_ / post_.

MATLAB twin: spec_window_alpha_features.m + spec_metrics_from_powers.m.
The definitions below are the parity contract for every spatial unit
(channel in MATLAB, ROI here):

    Trial x unit
        pow_*        mean of the filter-Hilbert power time course in the window
        paf_cog_hz   centre of gravity of the window's Welch PSD (8-12 Hz)
        ratios       from the three powers (see _metrics_from_powers)

    Subject x unit ("GA" here, "chan_summary" in MATLAB)
        pow_*        mean over trials of the trial powers
        paf_cog_hz   CoG of the mean-over-trials PSD
        ratios       from the mean powers (ratio of means)

V3.1.0 changes vs V3.0.0:
    - slow_alpha_frac removed: it equals (1 + sf_balance) / 2 exactly, so it
      duplicated every balance-index result. 10 feature families remain.

V3.0.0 changes vs V2.0.0:
    - Window masks use src_window_mask (quarter-sample tolerance), shared
      with MATLAB so both select the same samples.
    - Welch PSDs are computed for all trials x ROIs in one vectorised call.
    - Power scale is now the band-limited variance (see src_spectral V3.0.0).
"""

from __future__ import annotations

import numpy as np

from src_spectral import src_band_power_tc, src_cog_from_psd, src_psd_welch, src_window_mask

_EPS0 = 1e-12

_BAND_KEYS = ("slow", "fast", "alpha")

METRIC_NAMES = (
    "pow_slow_alpha", "pow_fast_alpha", "pow_alpha_total", "paf_cog_hz",
    "sf_ratio", "sf_logratio", "sf_balance",
    "rel_slow_alpha", "rel_fast_alpha", "psi_cog",
)


# ===================================================================
# BAND-POWER TIME COURSES (once per subject, full epoch)
# ===================================================================
def src_compute_band_power_tcs(
        tc: np.ndarray,
        sfreq: float,
        alpha: tuple[float, float],
        slow: tuple[float, float],
        fast: tuple[float, float],
        trans_bw: float = 1.5,
) -> dict:
    """
    Filter-Hilbert power time courses on the FULL epoch.

    Returns:
        Dict {"slow", "fast", "alpha"} -> (n_epochs, n_rois, n_times) power
    """
    bands = {"slow": slow, "fast": fast, "alpha": alpha}
    return {k: src_band_power_tc(tc, sfreq, lo, hi, trans_bw) for k, (lo, hi) in bands.items()}


# ===================================================================
# METRICS FROM POWERS (vectorised; works on scalars or arrays)
# ===================================================================
def _metrics_from_powers(ps, pf, pa, cog) -> dict:
    """The 10 alpha feature families (incl. psi_cog). NaN inputs propagate."""
    ps, pf, pa, cog = (np.asarray(v, dtype = float) for v in (ps, pf, pa, cog))

    sf_balance = (ps - pf) / (ps + pf + _EPS0)
    return {
        "pow_slow_alpha":  ps,
        "pow_fast_alpha":  pf,
        "pow_alpha_total": pa,
        "paf_cog_hz":      cog,
        "sf_ratio":        ps / (pf + _EPS0),
        "sf_logratio":     np.log(ps + _EPS0) - np.log(pf + _EPS0),
        "sf_balance":      sf_balance,
        "rel_slow_alpha":  ps / (pa + _EPS0),
        "rel_fast_alpha":  pf / (pa + _EPS0),
        "psi_cog":         sf_balance * (cog - 10.0),
    }


# ===================================================================
# ONE WINDOW: PER-TRIAL x ROI ROWS + GA ROWS + GA PSDs
# ===================================================================
def src_compute_window_alpha_features(
        power_tcs: dict,
        tc: np.ndarray,
        times: np.ndarray,
        tmin: float,
        tmax: float,
        sfreq: float,
        alpha: tuple[float, float],
        fmin: float,
        fmax: float,
        psd_window_sec: float = 2.0,
        df_target: float = 0.25,
        overlap: float = 0.5,
) -> tuple[list[dict], list[dict], dict]:
    """
    Compute the 10-family alpha feature set for one window [tmin, tmax] s.

    Args:
        power_tcs      : from src_compute_band_power_tcs (FULL epoch)
        tc             : (n_epochs, n_rois, n_times) FULL-epoch time courses
        times          : (n_times,) full-epoch time axis in seconds
        tmin, tmax     : window bounds in seconds (inclusive)
        sfreq          : Sampling frequency in Hz
        alpha          : (lo, hi) total alpha band in Hz, for the CoG
        fmin, fmax     : Welch PSD frequency range
        psd_window_sec : Welch segment length (clipped to the window length)
        df_target      : Welch frequency-grid spacing after zero-padding
        overlap        : Welch segment overlap fraction

    Returns:
        trial_rows    : one dict per (trial, ROI): {trial, roi_idx, <10 metrics>}
        ga_rows       : one dict per ROI: {roi_idx, <10 metrics>}
        ga_psd_by_roi : roi_idx -> (freqs, mean-over-trials PSD)
    """
    mask = src_window_mask(times, tmin, tmax, sfreq)
    if not np.any(mask):
        raise ValueError(
            f"Window [{tmin:.3f}, {tmax:.3f}] s has no samples in epoch range "
            f"[{times[0]:.3f}, {times[-1]:.3f}] s."
        )

    # Window-mean band power, (n_epochs, n_rois) per band
    pw = {k: power_tcs[k][:, :, mask].mean(axis = -1) for k in _BAND_KEYS}

    # Per-trial Welch PSDs (n_epochs, n_rois, n_freqs) and their CoG
    freqs, psd = src_psd_welch(tc[:, :, mask], sfreq, fmin, fmax, psd_window_sec, overlap, df_target)
    cog = src_cog_from_psd(freqs, psd, alpha)          # (n_epochs, n_rois)

    n_epochs, n_rois = cog.shape

    feat = _metrics_from_powers(pw["slow"], pw["fast"], pw["alpha"], cog)
    trial_rows: list[dict] = []
    for ei in range(n_epochs):
        for ri in range(n_rois):
            row = {k: float(feat[k][ei, ri]) for k in METRIC_NAMES}
            row["trial"] = ei + 1
            row["roi_idx"] = ri
            trial_rows.append(row)

    # GA: mean-over-trials powers; CoG from the mean-over-trials PSD
    ga_psd = np.nanmean(psd, axis = 0)                 # (n_rois, n_freqs)
    ga_cog = src_cog_from_psd(freqs, ga_psd, alpha)    # (n_rois,)
    ga_feat = _metrics_from_powers(
        np.nanmean(pw["slow"], axis = 0),
        np.nanmean(pw["fast"], axis = 0),
        np.nanmean(pw["alpha"], axis = 0),
        ga_cog,
    )
    ga_rows: list[dict] = []
    ga_psd_by_roi: dict = {}
    for ri in range(n_rois):
        row = {k: float(ga_feat[k][ri]) for k in METRIC_NAMES}
        row["roi_idx"] = ri
        ga_rows.append(row)
        ga_psd_by_roi[ri] = (freqs, ga_psd[ri])

    return trial_rows, ga_rows, ga_psd_by_roi
