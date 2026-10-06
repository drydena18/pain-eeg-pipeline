"""
src_fooof.py - FOOOF / specparam aperiodic decomposition.
V 3.0.0

ONE implementation used by both pipelines:
    - source:   src_compute_fooof_ga() on the subject x ROI trial-mean PSD
    - spectral: python/spectral/fooof_bridge.py imports src_fit_fooof() and
                fits the subject x channel trial-mean PSD sent from MATLAB

Works with specparam 2.x (SpectralModel), specparam 1.x and legacy fooof.

Outputs per fit
---------------
    fooof_offset, fooof_exponent, fooof_knee (NaN in 'fixed' mode)
    fooof_r2, fooof_error
    fooof_alpha_cf, fooof_alpha_pw, fooof_alpha_bw
        strongest peak (by power) with CF inside alpha_band; NaN if none

V3.0.0 changes vs V2.x:
    - Alpha band is a parameter (was hard-coded 8-12 Hz).
    - Adds fooof_r2 and fooof_error (needed for fit-quality filtering, and
      present in the MATLAB outputs).
    - src_fit_fooof_rows(): batch fit of an (n_rows, n_freqs) matrix, used by
      the MATLAB bridge.
"""

from __future__ import annotations

import numpy as np

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

FOOOF_NAMES = (
    "fooof_offset", "fooof_exponent", "fooof_knee", "fooof_r2", "fooof_error",
    "fooof_alpha_cf", "fooof_alpha_pw", "fooof_alpha_bw",
)


def fooof_available() -> bool:
    return FOOOF is not None


def fooof_package_name() -> str:
    return _FOOOF_PKG if _FOOOF_PKG is not None else "unavailable"


def _nan_metrics() -> dict:
    return {k: float("nan") for k in FOOOF_NAMES}


# =============================================================================
# VERSION-SAFE ACCESSORS
# =============================================================================
def _get_params(fm, keys: tuple[str, ...]):
    for key in keys:
        try:
            result = fm.get_params(key)
            if result is not None:
                return result
        except (AttributeError, TypeError, KeyError, IndexError, ValueError):
            continue
    return None


def _get_aperiodic(fm) -> np.ndarray:
    ap = _get_params(fm, ("aperiodic", "aperiodic_params"))
    if ap is None:
        raise RuntimeError("Could not retrieve aperiodic parameters from the fitted model.")
    return np.atleast_1d(ap)


def _get_peaks(fm) -> np.ndarray:
    pk = _get_params(fm, ("periodic", "peak_params"))
    if pk is None:
        return np.empty((0, 3))
    arr = np.atleast_2d(np.asarray(pk, dtype = float))
    return arr if arr.ndim == 2 and arr.shape[1] == 3 and np.all(np.isfinite(arr)) else np.empty((0, 3))


def _get_r2_error(fm) -> tuple[float, float]:
    r2 = getattr(fm, "r_squared_", None)
    err = getattr(fm, "error_", None)
    if r2 is not None and err is not None:
        return float(r2), float(err)
    try:                                      # specparam 2.x
        res = fm.results.metrics.results
        r2 = next((v for k, v in res.items() if "rsquared" in k), np.nan)
        err = next((v for k, v in res.items() if k.startswith("error")), np.nan)
        return float(r2), float(err)
    except Exception:
        return float("nan"), float("nan")


# =============================================================================
# SINGLE PSD FIT
# =============================================================================
def src_fit_fooof(freqs: np.ndarray, psd: np.ndarray, fooof_cfg: dict,
                  alpha_band: tuple[float, float] = (8.0, 12.0)) -> tuple[dict, object]:
    """
    Fit one PSD and return (metrics, fitted_model).

    fooof_cfg keys: aperiodic_mode, peak_width_limits, max_n_peaks,
                    min_peak_height, peak_threshold, freq_range
    """
    if FOOOF is None:
        raise RuntimeError("Neither 'specparam' nor 'fooof' is installed (pip install specparam).")

    mode = str(fooof_cfg.get("aperiodic_mode", "fixed"))
    fm = FOOOF(
        aperiodic_mode    = mode,
        peak_width_limits = tuple(fooof_cfg.get("peak_width_limits", [1.0, 12.0])),
        max_n_peaks       = int(fooof_cfg.get("max_n_peaks", 6)),
        min_peak_height   = float(fooof_cfg.get("min_peak_height", 0.1)),
        peak_threshold    = float(fooof_cfg.get("peak_threshold", 2.0)),
        verbose           = False,
    )
    fm.fit(np.asarray(freqs, dtype = float), np.asarray(psd, dtype = float),
           list(fooof_cfg.get("freq_range", [1.0, 40.0])))

    metrics = _nan_metrics()
    ap = _get_aperiodic(fm)
    metrics["fooof_offset"] = float(ap[0])
    if mode == "knee" and ap.size >= 3:
        metrics["fooof_knee"] = float(ap[1])
        metrics["fooof_exponent"] = float(ap[2])
    elif ap.size >= 2:
        metrics["fooof_exponent"] = float(ap[1])

    metrics["fooof_r2"], metrics["fooof_error"] = _get_r2_error(fm)

    peaks = _get_peaks(fm)
    in_band = peaks[(peaks[:, 0] >= alpha_band[0]) & (peaks[:, 0] <= alpha_band[1])] if peaks.size else peaks
    if in_band.size:
        cf, pw, bw = in_band[int(np.argmax(in_band[:, 1]))]
        metrics["fooof_alpha_cf"], metrics["fooof_alpha_pw"], metrics["fooof_alpha_bw"] = float(cf), float(pw), float(bw)

    return metrics, fm


def src_fit_fooof_rows(freqs: np.ndarray, psd_rows: np.ndarray, fooof_cfg: dict,
                       alpha_band: tuple[float, float] = (8.0, 12.0)) -> list[dict]:
    """Fit every row of an (n_rows, n_freqs) matrix; failures -> NaN + fail_reason."""
    out: list[dict] = []
    for i, psd in enumerate(np.atleast_2d(psd_rows)):
        try:
            metrics, _ = src_fit_fooof(freqs, psd, fooof_cfg, alpha_band)
        except Exception as e:
            metrics = _nan_metrics()
            metrics["fail_reason"] = str(e)
        metrics["row"] = i + 1
        out.append(metrics)
    return out


# =============================================================================
# BATCH SUBJECT x ROI FIT (source pipeline)
# =============================================================================
def src_compute_fooof_ga(psd_by_roi_idx: dict, n_rois: int, sub: int, fooof_cfg: dict,
                         alpha_band: tuple[float, float] = (8.0, 12.0)) -> tuple[list[dict], dict]:
    """
    Fit the subject x ROI trial-mean PSD for every ROI.

    Returns:
        fooof_rows : one dict per ROI {roi_idx, <FOOOF_NAMES>}
        fm_by_roi  : roi_idx -> fitted model (for plotting)
    """
    fooof_rows: list[dict] = []
    fm_by_roi: dict = {}
    for ri in range(n_rois):
        freqs, psd = psd_by_roi_idx[ri]
        try:
            metrics, fm = src_fit_fooof(freqs, psd, fooof_cfg, alpha_band)
            fm_by_roi[ri] = fm
        except Exception as e:
            print(f"[WARN] FOOOF failed for sub-{sub:03d} ROI index {ri}: {e}")
            metrics = _nan_metrics()
        fooof_rows.append({"roi_idx": ri, **metrics})
    return fooof_rows, fm_by_roi
