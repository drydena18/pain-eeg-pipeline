"""
src_prestim.py - Pre-stimulus Hilbert phase + temporal variability index.
V 3.0.0

MATLAB twins: spec_slow_phase_features.m (phase), spec_compute_tvi_alpha.m (TVI).

Metrics
-------
    slow_phase   Phase of the slow-alpha analytic signal at t = 0 (radians)
    sin_phase    sin(slow_phase)
    cos_phase    cos(slow_phase)
    TVI_alpha    MSSD / Var of the per-trial pre_sf_balance sequence
                 (subject x ROI scalar, GA CSV only)

The phase is read from the FULL-epoch filtered signal (not the cropped
pre-stim window) to keep filter/Hilbert edge effects away from t = 0.

V3.0.0 changes vs V2.0.0:
    - Band-pass uses the shared FIR (src_spectral.src_bandpass_filter with
      trans_bw), identical to MATLAB, and is vectorised over trials x ROIs.
    - Sample nearest t = 0 chosen identically in both languages (first
      minimum of |t|).
    - TVI_alpha variance now uses ddof = 1 (unbiased), matching MATLAB's
      var(x, 0). The MSSD / Var ratio ranges over [0, 4]; ~2 means
      trial-to-trial independence (the old docstring said [0, 2]).
"""

from __future__ import annotations

import numpy as np

from src_spectral import src_analytic_signal, src_bandpass_filter

_EPS = 1e-12


def _phase_at(tc: np.ndarray, times: np.ndarray, sfreq: float,
              slow: tuple[float, float], t_ref: float, trans_bw: float) -> np.ndarray:
    """Slow-band analytic phase at the sample nearest t_ref; (n_epochs, n_rois)."""
    idx = int(np.argmin(np.abs(np.asarray(times) - t_ref)))
    z = src_analytic_signal(src_bandpass_filter(tc, sfreq, slow[0], slow[1], trans_bw))
    return np.angle(z[..., idx])


def _phase_rows(phase: np.ndarray, suffix: str, with_trial: bool) -> list[dict]:
    rows: list[dict] = []
    if with_trial:
        n_epochs, n_rois = phase.shape
        for ei in range(n_epochs):
            for ri in range(n_rois):
                p = float(phase[ei, ri])
                rows.append({"trial": ei + 1, "roi_idx": ri,
                             f"slow_phase{suffix}": p,
                             f"sin_phase{suffix}": float(np.sin(p)),
                             f"cos_phase{suffix}": float(np.cos(p))})
    else:
        for ri, p in enumerate(np.asarray(phase).ravel()):
            p = float(p)
            rows.append({"roi_idx": ri,
                         f"slow_phase{suffix}": p,
                         f"sin_phase{suffix}": float(np.sin(p)),
                         f"cos_phase{suffix}": float(np.cos(p))})
    return rows


# ====================================================================
# PER-TRIAL PRE-STIMULUS PHASE
# ====================================================================
def src_compute_prestim_phase(
        tc_full: np.ndarray,
        times_full: np.ndarray,
        sfreq: float,
        slow: tuple[float, float],
        trans_bw: float = 1.5,
) -> list[dict]:
    """One row per (trial, ROI): {trial, roi_idx, slow_phase, sin_phase, cos_phase}."""
    phase = _phase_at(tc_full, times_full, sfreq, slow, 0.0, trans_bw)
    return _phase_rows(phase, "", with_trial = True)


# ====================================================================
# GRAND-AVERAGE PRE-STIMULUS PHASE (phase of the trial-mean time course)
# ====================================================================
def src_compute_ga_prestim_phase(
        tc_full_ga: np.ndarray,
        times_full: np.ndarray,
        sfreq: float,
        slow: tuple[float, float],
        trans_bw: float = 1.5,
) -> list[dict]:
    """tc_full_ga: (1, n_rois, n_times). One row per ROI."""
    phase = _phase_at(tc_full_ga, times_full, sfreq, slow, 0.0, trans_bw)[0]
    return _phase_rows(phase, "", with_trial = False)


# ====================================================================
# TEMPORAL VARIABILITY INDEX
# ====================================================================
def src_compute_tvi_alpha(bi_pre_sequence: np.ndarray) -> float:
    """
    TVI_alpha = MSSD / Var for one subject x ROI.

        MSSD = mean((b[k+1] - b[k])^2)
        Var  = unbiased sample variance (ddof = 1)

    NaN when fewer than 3 valid trials or Var ~ 0.
    """
    b = np.asarray(bi_pre_sequence, dtype = float)
    b = b[~np.isnan(b)]
    if len(b) < 3:
        return float("nan")
    mssd = float(np.mean(np.diff(b) ** 2))
    var = float(np.var(b, ddof = 1))
    if var < _EPS:
        return float("nan")
    return mssd / var
