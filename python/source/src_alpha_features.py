"""
src_alpha_features.py - Shared per-window alpha feature computation.
V2.0.0

Called once per window (whole / pre / post) from source_core.py, replacing
the previously-duplicated BI/LR/CoG computation that lived separately inside
src_prestim.py and src_poststim.py.

V2.0.0 changes vs V1.0.0:
    - Slow, fast, and total alpha POWER now come from filter-Hilbert bandpower
        time courses (src_spectral.src_band_power._tc) computed ONCE on the full
        epoch and averaged inside each window. Previously, they were integrated
        from a 0.5s Welch PSD whose ~1.95 Hz bins put a single bin inside each 2 Hz
        sub-band, so pow_slow_alpha / pow_fast_alpha were exactly 0.0 and
        every metric derived from them was 0.0 or NaN
    - paf_cog_hz is still a spectral centre of gravity, from a Welch PSD
        (2s segments, clipped to the window; zero-padded to 0.25 Hz bins).
    - NaN now PROPAGATES through the ratio metrics instead of being coerced
        to 0.0, so an un-estimable band shows up as NaN, not as a plausible 0.
    - GRAND AVERAGE (GA) FIX: GA power is not the MEAN OF PER-TRIAL POWERS,
        and the GA PSD (used for GA paf_cog_hz and FOOOF) is the MEAN OF
        PER-TRIAL PSDs. V1.0.0 averaed the time courses across trials FIRST and
        then took the PSD, which measures only the phase-locked (evoked) part
        of the signal.
    - Trial rows and GA rows are returned by ONE call per window, so each
        window's PSDs are computed once.
"""

from __future__ import annotations

import numpy as np

from src_spectral import src_psd_welch, src_band_power_tc

_EPS0 = 1e-12

# Band keys used throughout this module
_BAND_KEYS = ("slow", "fast", "alpha")


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
    Filter-Hilbert power time courses for thr slow, fast, and total alpha
    bands on the FULL epoch.

    Args:
        tc                  : (n_epochs, n_rois, n_times) full-epoch source time courses
        sfreq               : Sampling frequency in Hz
        alpha, slow, fast   : (lo, hi) band bounds in Hz
        trans_bw            : FIR transition bandwidth in Hz (see src_band_power_tc)

    Returns:
        Dict {"slow", "fast", "alpha"} -> (n_epochs, n_rois, n_times) power
    """
    bands = {"slow": slow, "fast": fast, "alpha": alpha}
    return {k: src_band_power_tc(tc, sfreq, lo, hi, trans_bw) for k, (lo, hi) in bands.items()}


# ===================================================================
# METRICS FROM POWERS (vectorized; works on scalars or arrays)
# ===================================================================
def _metrics_from_powers(ps, pf, pa, cog) -> dict:
    """
    The 10 alpha metrics + psi_cog from slow / fast / total power and CoG
    NaN inputs propagate to NaN outputs.
    """
    ps, pf, pa, cog = (np.asarray(v, dtype = float) for v in (ps, pf, pa, cog))

    sf_ratio        = ps / (pf + _EPS0)
    sf_logratio     = np.log(ps + _EPS0) - np.log(pf + _EPS0)
    sf_balance      = (ps - pf) / (ps + pf + _EPS0)
    slow_alpha_frac = ps / (ps + pf + _EPS0)
    rel_slow_alpha  = ps / (pa + _EPS0)
    rel_fast_alpha  = pf / (pa + _EPS0)
    psi_cog         = sf_balance / (cog + 10.0)

    return {
        "pow_slow_alpha":   ps,
        "pow_fast_alpha":   pf,
        "pow_total_alpha":  pa,
        "paf_cog_hz":       cog,
        "sf_ratio":         sf_ratio,
        "sf_logratio":      sf_logratio,
        "sf_balance":       sf_balance,
        "slow_alpha_frac":  slow_alpha_frac,
        "rel_slow_alpha":   rel_slow_alpha,
        "rel_fast_alpha":   rel_fast_alpha,
        "psi_cog":          psi_cog,
    }


def _cog_from_psd(freqs: np.ndarray, psd: np.ndarray, alpha: tuple[float, float]) -> float:
    """
    Spectral centre of gravity over the alpha band; NaN if not estimable.
    """
    idx = (freqs >= alpha[0]) & (freqs <= alpha[1])
    if np.count_nonzero(idx) < 2:
        return float("nan")
    den = float(np.trapezoid(psd[idx], freqs[idx]))
    if not np.isfinite(den) or den <= 0:
        return float("nan")
    return float(np.trapezoid(psd[idx] * freqs[idx], freqs[idx])) / den


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
) -> tuple[list[dict], list[dict], dict]:
    """
    Compute the 11-metric alpha freature set for one window [tmin, tmax] s.

    Power metrics: mean of the filter-Hilbert power time courses inside the
    window. paf_cog_hz: centre of gravity of a Welch PSD of the cropped window.

    Args:
        power_tcs       : from src_compute_band_power_tcs (FULL epoch)
        tc              : (n_epochs, n_rois, n_times) FULL-epoch time courses
                            (cropped here for the Welch PSD)
        times           : (n_times,) full-epoch time axis in seconds
        tmin, tmax      : window bounds in seconds (inclusive)
        sfreq           : Sampling frequency in Hz
        alpha           : (lo, hi) total alpha band in Hz, for the CoG
        fmin, fmax      : Welch PSD frequency range
        psd_window_sec  : Welch segment length (clipped to the window length)
        df_target       : Welch frequency-grid spacing after zero-padding

    Returns:
        trial_rows      : list of dicts, one per (trial, ROI):
                            {trial, roi_idx, <11 metrics>} (unprefixed)
        ga_rows         : list of dicts, one per ROI: {roi_idx, <11 metrics>}
                            from mean-over-trials powers and the mean-over-trials PSD
        ga_psd_by_roi   : dict roi_idx -> (freqs, mean-over-trials PSD), for
                            FOOOF / plotting
    """
    mask = (times >= tmin) & (times <= tmax)
    if not np.any(mask):
        raise ValueError(
            f"Window [{tmin:.3f}, {tmax:.3f}] s has no samples in epoch range "
            f"[{times[0]:.3f}, {times[-1]:.3f}] s."
            )

    # Window-mean band power, (n_epochs, n_rois) per band
    pw = {k: power_tcs[k][:, :, mask].mean(axis = -1) for k in _BAND_KEYS}

    # Per-trial Welch PSDs for the CoG, and their trial mean for the GA
    tc_win = tc[:, :, mask]
    n_epochs, n_rois, _ = tc_win.shape
    cog = np.full((n_epochs, n_rois), np.nan)
    ga_psd_by_roi: dict = {}

    for ri in range(n_rois):
        psds = []
        freqs = None
        for ei in range(n_epochs):
            freqs, psd = src_psd_welch(tc_win[ei, ri, :], sfreq, fmin, fmax, psd_window_sec, overlap = 0.5, df_target = df_target)
            psds.append(psd)
            cog[ei, ri] = _cog_from_psd(freqs, psd, alpha)
        ga_psd_by_roi[ri] = (freqs, np.mean(np.vstack(psds), axis = 0))

    # Per-trial rows
    feat = _metrics_from_powers(pw["slow"], pw["fast"], pw["alpha"], cog)
    trial_rows = list[dict] = []
    for ei in range(n_epochs):
        for ri in range(n_rois):
            row = {k: float(v[ei, ri]) for k, v in feat.items()}
            row["trial"] = ei + 1
            row["roi_idx"] = ri
            trial_rows.append(row)

    # GA rows: mean-over-trials powers; CoG from the mean-over-trials PSD
    ga_rows: list[dict] = []
    for ri in range(n_rois):
        freqs, ga_psd = ga_psd_by_roi[ri]
        ga_feat = _metrics_from_powers(
            np.nanmean(pw["slow"][:, ri]),
            np.nanmean(pw["fast"][:, ri]),
            np.nanmean(pw["alpha"][:, ri]),
            _cog_from_psd(freqs, ga_psd, alpha),
        )
        row = {k: float(v) for k, v in ga_feat.items()}
        row["roi_idx"] = ri
        ga_rows.append(row)

    return trial_rows, ga_rows, ga_psd_by_roi