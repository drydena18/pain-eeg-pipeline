"""
src_lep.py - Laser-evoked potential (LEP) feature extraction.
V 3.0.0

MATLAB twin: spec_lep_features.m. Same windows, same peak rule, same column
names in both pipelines.

Features (per trial x unit, and per unit on the trial-mean waveform)
--------------------------------------------------------------------
    n2_amp     most negative value in the N2 window
    n2_lat_ms  its latency (first occurrence), ms
    p2_amp     most positive value in the P2 window
    p2_lat_ms  its latency (first occurrence), ms
    n2p2_amp   p2_amp - n2_amp
    n2_mean    mean amplitude in the N2 window
    p2_mean    mean amplitude in the P2 window
    n_trials   (GA rows only) number of trials averaged

Units: uV in channel space; arbitrary sLORETA units in source space. A
sign-aligned ROI time course's polarity depends on the dominant source
orientation, so source "N2"/"P2" are the most negative/positive deflections
in their windows, not necessarily the scalp components' polarity.

V3.0.0 changes vs V2.x:
    - Searches the FULL epoch with src_window_mask (quarter-sample
      tolerance) instead of a pre-cropped post window; results are identical
      while the windows sit inside the post window, and the masks now match
      MATLAB exactly.
    - Vectorised over trials x ROIs.
"""

from __future__ import annotations

import numpy as np

from src_spectral import src_window_mask

LEP_NAMES = ("n2_amp", "n2_lat_ms", "p2_amp", "p2_lat_ms", "n2p2_amp", "n2_mean", "p2_mean")


def _lep_arrays(x: np.ndarray, times: np.ndarray, sfreq: float,
                n2_window: tuple[float, float], p2_window: tuple[float, float]) -> dict:
    """x: (..., n_times). Returns dict of (...)-shaped arrays."""
    times = np.asarray(times, dtype = float)
    out: dict = {}
    for name, win, fn in (("n2", n2_window, np.argmin), ("p2", p2_window, np.argmax)):
        mask = src_window_mask(times, win[0], win[1], sfreq)
        if not np.any(mask):
            nan = np.full(x.shape[:-1], np.nan)
            out[f"{name}_amp"], out[f"{name}_lat_ms"], out[f"{name}_mean"] = nan, nan.copy(), nan.copy()
            continue
        xw = x[..., mask]
        tw = times[mask]
        idx = fn(xw, axis = -1)
        out[f"{name}_amp"] = np.take_along_axis(xw, idx[..., None], axis = -1)[..., 0]
        out[f"{name}_lat_ms"] = tw[idx] * 1000.0
        out[f"{name}_mean"] = xw.mean(axis = -1)
    out["n2p2_amp"] = out["p2_amp"] - out["n2_amp"]
    return out


# ====================================================================
# PER-TRIAL LEP EXTRACTION
# ====================================================================
def src_compute_lep_trial(
        tc: np.ndarray,
        times: np.ndarray,
        sfreq: float,
        n2_window: tuple[float, float],
        p2_window: tuple[float, float],
) -> list[dict]:
    """tc: (n_epochs, n_rois, n_times) full epoch. One row per (trial, ROI)."""
    a = _lep_arrays(tc, times, sfreq, n2_window, p2_window)
    n_epochs, n_rois = tc.shape[:2]
    rows: list[dict] = []
    for ei in range(n_epochs):
        for ri in range(n_rois):
            row = {"trial": ei + 1, "roi_idx": ri}
            row.update({k: float(a[k][ei, ri]) for k in LEP_NAMES})
            rows.append(row)
    return rows


# ====================================================================
# GRAND-AVERAGE LEP EXTRACTION (trial-mean waveform)
# ====================================================================
def src_compute_lep_ga(
        tc: np.ndarray,
        times: np.ndarray,
        sfreq: float,
        n2_window: tuple[float, float],
        p2_window: tuple[float, float],
) -> list[dict]:
    """One row per ROI from the mean-over-trials time course, plus n_trials."""
    n_epochs = tc.shape[0]
    a = _lep_arrays(tc.mean(axis = 0), times, sfreq, n2_window, p2_window)
    rows: list[dict] = []
    for ri in range(tc.shape[1]):
        row = {"roi_idx": ri, "n_trials": n_epochs}
        row.update({k: float(a[k][ri]) for k in LEP_NAMES})
        rows.append(row)
    return rows
