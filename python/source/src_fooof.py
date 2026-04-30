"""
src_fooof.py - FOOOF (Fitting Oscillations & One-Over-F) helpers

Uses either the 'specparam' package (new name, v2+) or the legacy 'fooof'
package transparently. An ImportError at module load is deferred rather than
fatal: functions check FOOOF availability and raise clearly at call time

Covers:
    - Single PSD fit with aperiodic parameter and alpha peak extraction
    - Batch GA fit over all ROIs
"""

from __future__ import annotations

import numpy as np

# -- specparam / fooof compatibility shim -----
try:
    from specparam import SpectralModel as FOOOF
    _FOOOF_PKG = "specparam"
except ImportError:
    try:
        from fooof import FOOOF
        _FOOOF_PKG = "fooof"
    except ImportError:
        FOOOF = None
        _FOOOF_PKG = None

def fooof_available() -> bool:
    """
    Return True if either specparam or fooof is importable
    """
    return FOOOF is not None

def fooof_package_name() -> str:
    """
    Return the name of the installed package, or 'unavailable'
    """
    return _FOOOF_PKG if _FOOOF_PKG is not None else "unavailable"

# ====================================================================
# SINGLE PSD FIT
# ====================================================================
def src_fit_fooof(
        freqs: np.ndarray,
        psd: np.ndarray,
        fooof_cfg: dict,
) -> tuple[dict, object]:
    """
    Fit a FOOOF model to a single PSD vector and extract aperiodic parameters
    plus the dominant alpha-band peak

    Args:
        freqs     : Frequency axis, shape (n_freqs,).
        psd       : Power spectral density values, shape (n_freqs,).
        fooof_cfg : Dict with keys:
                        aperiodic_mode     : "fixed" or "knee"
                        peak_width_limits  : [min, max] in Hz
                        max_n_peaks        : int
                        min_peak_height    : float
                        peak_threshold     : float
                        freq_range         : [fmin, fmax] for fitting

    Returns:
        metrics : Dict with keys:
                      fooof_offset, fooof_exponent, fooof_knee  (aperiodic)
                      fooof_alpha_cf, fooof_alpha_pw, fooof_alpha_bw  (peak)
        fm      : Fitted FOOOF model object (for plotting).

    Raises:
        RuntimeError is neither specparam not fooof is installed
    """
    if FOOOF is None:
        return RuntimeError(
            "FOOOF fitting requested but neither 'specparam' not 'fooof' is installed.\n"
            "Install with: pip install specparam (or: pip install fooof)"
        )
    
    mode = fooof_cfg.get("aperiodic_mode", "fixed")

    fm = FOOOF(
        aperiodic_mode = mode,
        peak_width_limits = tuple(fooof_cfg.get("peak_width_limits", [1.0, 12.0])),
        max_n_peaks = int(fooof_cfg.get("max_n_peaks", 6)),
        min_peak_height = float(fooof_cfg.get("min_peak_height", 0.1)),
        peak_threshold = float(fooof_cfg.get("peak_threshold", 2.0)),
        verbose = False,
    )
    freq_range = fooof_cfg.get("freq_range", [1.0, 40.0])
    fm.fit(freqs, psd, freq_range)

    # -- Aperiodic parameters -----
    ap = fm.get_params("aperiodic_params")
    metrics: dict = {}

    if mode == "fixed":
        metrics["fooof_offset"] = float(ap[0])
        metrics["fooof_exponent"] = float(ap[1])
        metrics["fooof_knee"] = float("nan")
    else:
        metrics["fooof_offset"] = float(ap[0])
        metrics["fooof_knee"] = float(ap[1])
        metrics["fooof_exponent"] = float(ap[2])

    # -- Alpha peak extraction -----
    # Pick the strongest peak whose centre frequency (CF) is inside [8, 12] Hz.
    try:
        peaks = np.atleast_2d(fm.get_params("peak_params"))
        alpha_peaks = [
            (cf, pw, bw) for cf, pw, bw in peaks
            if 8.0 <= cf <= 12.0
        ]
        if alpha_peaks:
            alpha_peaks.sort(key = lambda t: t[1], reverse = True)
            cf, pw, bw = alpha_peaks[0]
            metrics["fooof_alpha_cf"] = float(cf)
            metrics["fooof_alpha_pw"] = float(pw)
            metrics["fooof_alpha_bw"] = float(bw)
        else:
            metrics["fooof_alpha_cf"] = float("nan")
            metrics["fooof_alpha_pw"] = float("nan")
            metrics["fooof_alpha_bw"] = float("nan")
    except Exception:
        metrics["fooof_alpha_cf"] = float("nan")
        metrics["fooof_alpha_pw"] = float("nan")
        metrics["fooof_alpha_bw"] = float("nan")

    return metrics, fm

# ====================================================================
# BATCH GA FIR OVER ALL ROIs
# ====================================================================
def src_compute_fooof_ga(
        psd_by_roi_idx: dict,
        n_rois: int,
        sub: int,
        fooof_cfg: dict,
) -> tuple[list[dict], dict]:
    """
    Run FOOOF on the grand-average PSD for each ROI

    Individual ROI failures are aught and logged as NaN rows rather than
    aborting the entire subject

    Args:
        psd_by_roi_idx  : Dict mapping roi_idx (int) -> (freqs, psd)
        n_rois          : Total number of ROIs (used to iterate in order)
        sub             : Subject ID (written into output rows)
        fooof_cfg       : FOOOF config dict (see src_fit_fooof)

    Returns:
        fooof_rows : Lost of dicts (one per ROI) for the FOOOF GA CSV
        fm_by_roi  : Dict mapping roi_idx -> fitted FOOOF model (for plotting)
    """
    _nan_metrics = {
        "fooof_offset": float("nan"),
        "fooof_exponent": float("nan"),
        "fooof_knee": float("nan"),
        "fooof_alpha_cf": float("nan"),
        "fooof_alpha_pw": float("nan"),
        "fooof_alpha_bw": float("nan"),
    }

    fooof_rows: list[dict] = []
    fm_by_roi: dict = {}

    for ri in range(n_rois):
        freqs, psd = psd_by_roi_idx[ri]
        try:
            metrics, fm = src_fit_fooof(freqs, psd, fooof_cfg)
            fm_by_roi[ri] = fm
        except Exception as e:
            print(f"[WARN] FOOOF failed for ROI index {ri}: {e}")
            metrics = dict(_nan_metrics)

        fooof_rows.append({"subject": sub, "roi_idx": ri, **metrics})

    return fooof_rows, fm_by_roi