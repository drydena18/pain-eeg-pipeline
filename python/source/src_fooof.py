"""
src_fooof.py - FOOOF (Fitting Oscillations & One-Over-F) helpers

Uses either the 'specparam' package (v2+, new name) or the legacy 'fooof'
package transparently. An ImportError at module load is deferred rather than
fatal: functions check FOOOF availability and raise clearly at call time.

specparam 2.0 API note
-----------------------
specparam 2.0 renamed the get_params component strings:
    old (fooof / specparam 1.x)  -> new (specparam 2.x)
    "aperiodic_params"           -> "aperiodic"
    "peak_params"                -> "periodic"

Both are tried in order so the same code works with either package version.

Covers:
    - Single PSD fit with aperiodic parameter and alpha peak extraction
    - Batch GA fit over all ROIs
"""

from __future__ import annotations

import numpy as np

# -- specparam / fooof compatibility shim ─────────────────────────────────────
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
    """Return True if either specparam or fooof is importable."""
    return FOOOF is not None

def fooof_package_name() -> str:
    """Return the name of the installed package, or 'unavailable'."""
    return _FOOOF_PKG if _FOOOF_PKG is not None else "unavailable"


# =============================================================================
# API-VERSION-SAFE GET_PARAMS HELPERS
# =============================================================================

def _get_aperiodic(fm) -> np.ndarray:
    """
    Return the aperiodic parameter array from a fitted FOOOF/SpectralModel.

    Tries specparam 2.x key ("aperiodic") first, then falls back to the
    fooof / specparam 1.x key ("aperiodic_params").
    """
    for key in ("aperiodic", "aperiodic_params"):
        try:
            result = fm.get_params(key)
            if result is not None:
                return np.atleast_1d(result)
        except (AttributeError, TypeError, KeyError):
            continue
    raise RuntimeError(
        "Could not retrieve aperiodic parameters from the fitted model. "
        "Neither 'aperiodic' nor 'aperiodic_params' succeeded."
    )


def _get_peaks(fm) -> np.ndarray:
    """
    Return the peak parameter array from a fitted FOOOF/SpectralModel.

    Each row is [CF, PW, BW]. Returns an empty (0, 3) array if no peaks
    were fitted or if retrieval fails.

    Tries specparam 2.x key ("periodic") first, then the legacy key
    ("peak_params").
    """
    for key in ("periodic", "peak_params"):
        try:
            result = fm.get_params(key)
            if result is not None:
                arr = np.atleast_2d(result)
                if arr.shape[1] == 3:
                    return arr
        except (AttributeError, TypeError, KeyError, IndexError):
            continue
    return np.empty((0, 3))  # no peaks


# =============================================================================
# SINGLE PSD FIT
# =============================================================================

def src_fit_fooof(
        freqs: np.ndarray,
        psd: np.ndarray,
        fooof_cfg: dict,
) -> tuple[dict, object]:
    """
    Fit a FOOOF/SpectralModel to a single PSD vector and extract aperiodic
    parameters plus the dominant alpha-band peak.

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
        fm      : Fitted model object (for plotting).

    Raises:
        RuntimeError if neither specparam nor fooof is installed.
    """
    if FOOOF is None:
        raise RuntimeError(
            "FOOOF fitting requested but neither 'specparam' nor 'fooof' is installed.\n"
            "Install with: pip install specparam"
        )

    mode = fooof_cfg.get("aperiodic_mode", "fixed")

    fm = FOOOF(
        aperiodic_mode    = mode,
        peak_width_limits = tuple(fooof_cfg.get("peak_width_limits", [1.0, 12.0])),
        max_n_peaks       = int(fooof_cfg.get("max_n_peaks", 6)),
        min_peak_height   = float(fooof_cfg.get("min_peak_height", 0.1)),
        peak_threshold    = float(fooof_cfg.get("peak_threshold", 2.0)),
        verbose           = False,
    )
    freq_range = fooof_cfg.get("freq_range", [1.0, 40.0])
    fm.fit(freqs, psd, freq_range)

    # -- Aperiodic parameters ─────────────────────────────────────────────────
    ap = _get_aperiodic(fm)
    metrics: dict = {}

    if mode == "fixed":
        # ap = [offset, exponent]
        metrics["fooof_offset"]   = float(ap[0])
        metrics["fooof_exponent"] = float(ap[1])
        metrics["fooof_knee"]     = float("nan")
    else:
        # ap = [offset, knee, exponent]
        metrics["fooof_offset"]   = float(ap[0])
        metrics["fooof_knee"]     = float(ap[1])
        metrics["fooof_exponent"] = float(ap[2])

    # -- Alpha peak extraction ─────────────────────────────────────────────────
    # Pick the strongest peak whose centre frequency (CF) falls in [8, 12] Hz.
    try:
        peaks = _get_peaks(fm)
        alpha_peaks = [
            (cf, pw, bw) for cf, pw, bw in peaks
            if 8.0 <= cf <= 12.0
        ]
        if alpha_peaks:
            alpha_peaks.sort(key=lambda t: t[1], reverse=True)
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


# =============================================================================
# BATCH GA FIT OVER ALL ROIs
# =============================================================================

def src_compute_fooof_ga(
        psd_by_roi_idx: dict,
        n_rois: int,
        sub: int,
        fooof_cfg: dict,
) -> tuple[list[dict], dict]:
    """
    Run FOOOF/SpectralModel on the grand-average PSD for each ROI.

    Individual ROI failures are caught and logged as NaN rows rather than
    aborting the entire subject.

    Args:
        psd_by_roi_idx  : Dict mapping roi_idx (int) -> (freqs, psd).
        n_rois          : Total number of ROIs (used to iterate in order).
        sub             : Subject ID (written into output rows).
        fooof_cfg       : FOOOF config dict (see src_fit_fooof).

    Returns:
        fooof_rows : List of dicts (one per ROI) for the FOOOF GA CSV.
        fm_by_roi  : Dict mapping roi_idx -> fitted model (for plotting).
    """
    _nan_metrics = {
        "fooof_offset":   float("nan"),
        "fooof_exponent": float("nan"),
        "fooof_knee":     float("nan"),
        "fooof_alpha_cf": float("nan"),
        "fooof_alpha_pw": float("nan"),
        "fooof_alpha_bw": float("nan"),
    }

    fooof_rows: list[dict] = []
    fm_by_roi:  dict = {}

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