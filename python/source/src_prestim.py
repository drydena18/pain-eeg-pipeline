"""
src_prestim.py - Pre-stimulus Hilbert phase + temporal variability index.

V2.0.0 changes vs V1.x:
    - REMOVED the alpha-feature computation (pow_slow/pow_fast/pow_alpha,
      BI_pre, LR_pre, CoG_pre, psi_cog). These are now computed generically
      for all three windows (whole/pre/post) by
      src_alpha_features.src_compute_window_alpha_features /
      src_compute_ga_window_alpha_features, called directly from
      source_core.py, then prefixed with src_add_prefix. This eliminates
      the previous duplication of the same BI/LR/CoG formulas between this
      file and src_poststim.py.
    - This file now owns ONLY:
        - Hilbert instantaneous phase at stimulus onset (slow_phase,
          sin_phase, cos_phase) — unchanged from V1.x, still derived from
          the FULL epoch (tc_full) to avoid Hilbert edge effects, not from
          the cropped pre-stim window.
        - src_compute_tvi_alpha — unchanged; still a generic nMSSD over any
          1-D sequence. The caller now passes the pre_sf_balance sequence
          (renamed from BI_pre to match the new whole_/pre_/post_ naming
          scheme) rather than a bi_pre sequence, but the function itself
          didn't need to change.

Metrics Implemented
--------------------
Hilbert instantaneous phase at stimulus onset [Metric 4]
    slow_phase      Phase of the slow-alpha (8-10 Hz) analytic signal at t = 0
    sin_phase       sin(slow_phase)     - linear GAMM regressor
    cos_phase       cos(slow_phase)     - linear GAMM regressor

Temporal variability index [Metric 6 - TVI_alpha]
    TVI_alpha       nMSSD of the per-trial pre_sf_balance sequence.
                    This is a per-subject scalar, written to the GA CSV.
                    Stays pre-stim-only by design (out of scope for the
                    whole/pre/post/delta generalization).

References
-----------
    Furman et al. (2019, 2020); Tu et al. (2016); Nickel et al. (2022);
    Busch et al. (2009); Mathewson et al. (2009); Li et al. (2018)
"""

from __future__ import annotations

import numpy as np
from scipy.signal import hilbert

from src_spectral import src_bandpass_filter

_EPS = 1e-12


# ====================================================================
# PER-TRIAL PRE-STIMULUS PHASE
# ====================================================================
def src_compute_prestim_phase(
        tc_full: np.ndarray,
        times_full: np.ndarray,
        sfreq: float,
        slow: tuple[float, float],
) -> list[dict]:
    """
    Compute the Hilbert instantaneous phase of the slow-alpha (8-10 Hz)
    signal at stimulus onset (t = 0) for every (trial, ROI).

    Extracted from the FULL epoch time course (tc_full) rather than the
    cropped pre-stimulus window, to avoid edge effects from the Hilbert
    transform. The phase is then read at the sample nearest to t = 0.

    Args:
        tc_full    : (n_epochs, n_rois, n_times) — full epoch.
        times_full : (n_times,) — full epoch time axis in seconds.
        sfreq      : Sampling frequency in Hz.
        slow       : (lo, hi) slow-alpha sub-band bounds in Hz.

    Returns:
        List of dicts, one per (trial, ROI):
            {trial, roi_idx, slow_phase, sin_phase, cos_phase}
    """
    n_epochs, n_rois, _ = tc_full.shape
    t0_idx = int(np.argmin(np.abs(times_full)))

    rows: list[dict] = []
    for ei in range(n_epochs):
        for ri in range(n_rois):
            x_full = tc_full[ei, ri, :]
            try:
                x_filt = src_bandpass_filter(x_full, sfreq, slow[0], slow[1])
                analytic = hilbert(x_filt)
                slow_phase = float(np.angle(analytic[t0_idx]))
            except Exception:
                slow_phase = float("nan")

            rows.append({
                "trial":      ei + 1,
                "roi_idx":    ri,
                "slow_phase": slow_phase,
                "sin_phase":  float(np.sin(slow_phase)) if not np.isnan(slow_phase) else float("nan"),
                "cos_phase":  float(np.cos(slow_phase)) if not np.isnan(slow_phase) else float("nan"),
            })

    return rows


# ====================================================================
# GRAND-AVERAGE PRE-STIMULUS PHASE
# ====================================================================
def src_compute_ga_prestim_phase(
        tc_full_ga: np.ndarray,
        times_full: np.ndarray,
        sfreq: float,
        slow: tuple[float, float],
) -> list[dict]:
    """
    Same as src_compute_prestim_phase but for the GA (trial-averaged) time
    course.

    Args:
        tc_full_ga : (1, n_rois, n_times) — GA mean time course, keepdims=True.
        Other args same as src_compute_prestim_phase.

    Returns:
        List of dicts, one per ROI (no 'trial' key): {roi_idx, slow_phase, sin_phase, cos_phase}
    """
    _, n_rois, _ = tc_full_ga.shape
    t0_idx = int(np.argmin(np.abs(times_full)))

    rows: list[dict] = []
    for ri in range(n_rois):
        x_full = tc_full_ga[0, ri, :]
        try:
            x_filt = src_bandpass_filter(x_full, sfreq, slow[0], slow[1])
            analytic = hilbert(x_filt)
            slow_phase = float(np.angle(analytic[t0_idx]))
        except Exception:
            slow_phase = float("nan")

        rows.append({
            "roi_idx":    ri,
            "slow_phase": slow_phase,
            "sin_phase":  float(np.sin(slow_phase)) if not np.isnan(slow_phase) else float("nan"),
            "cos_phase":  float(np.cos(slow_phase)) if not np.isnan(slow_phase) else float("nan"),
        })

    return rows


# ====================================================================
# TEMPORAL VARIABILITY INDEX (TVI_alpha / nMSSD)
# ====================================================================
def src_compute_tvi_alpha(bi_pre_sequence: np.ndarray) -> float:
    """
    Compute the temporal variability index (TVI_alpha) from a per-trial
    balance-index sequence for a single subject x ROI.

    Unchanged from V1.x. Caller now passes the pre_sf_balance sequence
    (previously BI_pre) — see source_core.py.

        MSSD = mean( (b[k+1] - b[k])^2 ) for k = 1 ... K - 1
        Var  = variance( b )
        TVI_alpha = MSSD / Var

    This isolates the *temporal autocorrelation structure* of the pre-stimulus
    state rather than its total amplitude variability.

    Range: [0, 2]. Near 0 = rigid (slowly varying); near 2 = maximally
    alternating (each trial flips sign relative to the previous).

    Args:
        bi_pre_sequence : 1-D array of per-trial balance-index values, shape (K,)

    Returns:
        TVI_alpha as a float, or NaN if K < 3 (too few trials to be meaningful).
    """
    b = np.asarray(bi_pre_sequence, dtype=float)
    b = b[~np.isnan(b)]
    K = len(b)
    if K < 3:
        return float("nan")

    mssd = float(np.mean(np.diff(b) ** 2))
    var = float(np.var(b, ddof=0))

    if var < _EPS:
        return float("nan")

    return mssd / var