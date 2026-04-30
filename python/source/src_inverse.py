"""
src_inverse.py - MNE forward solution and sLORETA inverse pipeline.

Covers:
    - Building (or loading from cache) the forward solution
    - Computing the noise covariance from a baseline window
    - Assembling the inverse operator
    - Applying sLORETA to all epochs -> per-trial source time course
    - Applying sLORETA to the grand-average evoked response
    - Extracting ROI time courses from source estimates
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np
import mne

from src_io import src_logmsg

# ====================================================================
# FORWARD SOLUTION + NOISE COVARIANCE + INVERSE OPERATOR
# ====================================================================
def src_make_inverse_operator(
        epochs: mne.Epochs,
        bem_sol: str,
        trans: str,
        src_space: str,
        noise_tmin: float,
        noise_tmax: float,
        mindist_mm: float,
        loose: float,
        depth: float,
        fwd_cache_path: Optional[str],
        logf,
) -> tuple:
    """
    Build (or reload from cache) the forward solution, compute the noise
    covariance from the pre-stimulus baseline, and assemble the inverse operator.

    The forward solution is expensive but deterministic for a given electrode
    layout and source space, so it is written to <fwd_cache_path> and reloaded
    on subsequent runs instead of being recomputed.

    Args:
        epochs         : MNE Epochs object (already loaded).
        bem_sol        : Path to the fsaverage BEM solution .fif file.
        trans          : Path to the EEG->MRI transform .fif file.
        src_space      : Path to the fsaverage source space .fif file.
        noise_tmin     : Start of baseline window for noise covariance (s).
        noise_tmax     : End of baseline window for noise covariance (s).
        mindist_mm     : Minimum source-to-skull distance in mm (typically 5).
        loose          : Loose orientation constraint (0 = fixed, 1 = free).
        depth          : Depth weighting exponent (0 = none, 0.8 = standard).
        fwd_cache_path : Where to save/load the forward solution. Pass None to
                         disable caching.
        logf           : Log file handle from src_open_log().

    Returns:
        inv : MNE InverseOperator
        fwd : MNE ForwardSolution (kept for extract_label_time_course)
    """
    # -- Forward solution -----
    if fwd_cache_path and os.path.exists(fwd_cache_path):
        src_logmsg(logf, "[FWD] Loading cached forward solution: %s", fwd_cache_path)
        fwd = mne.read_forward_solution(fwd_cache_path, verbose = "ERROR")
    else:
        src_logmsg(logf, "[FWD] Computing forward solution (mindist = %.1f mm)...", mindist_mm)
        fwd = mne.make_forward_solution(
            info = epochs.info,
            trans = trans,
            src = src_space,
            bem = bem_sol,
            eeg = True,
            meg = False,
            mindist = mindist_mm,
            verbose = "ERROR",
        )
        if fwd_cache_path:
            os.makedirs(os.path.dirname(fwd_cache_path), exist_ok = True)
            mne.write_forward_solution(fwd_cache_path, fwd, overwrite = True, verbose = "ERROR")
            src_logmsg(logf, "[FWD] Cached to: %s", fwd_cache_path)
        
    # -- Noise covariance -----
    src_logmsg(logf, "[COV] Noise covariance from baseline (%.3f - %.3f s)...",
               noise_tmin, noise_tmax)
    noise_cov = mne.compute_covariance(
        epochs,
        tmin = noise_tmin,
        tmax = noise_tmax,
        method = "empirical",
        rank = None,
        verbose = "ERROR",
    )

    # -- Inverse Operator -----
    inv = mne.minimum_norm.make_inverse_operator(
        epochs.info,
        fwd,
        noise_cov,
        loose = loose,
        depth = depth,
        verbose = "ERROR",
    )
    src_logmsg(logf, "[INV] Operator ready (loose = %.2f, depth = %.2f).", loose, depth)
    return inv, fwd

# ====================================================================
# APPLY INVERSE
# ====================================================================
def src_apply_inverse_epochs(
        epochs: mne.Epochs,
        inv,
        fwd,
        labels: list,
        lambda2: float,
        pick_ori: Optional[str],
        roi_extract_mode: str,
        logf,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply sLORETA to every epoch and extract ROI (region of interest) time courses.

    Args:
        epochs          : Loaded MNE Epochs.
        inv             : Inverse operator from src_make_inverse_operator().
        fwd             : Forward solution from src_make_inverse_operator().
        labels          : List of mne.Label objects (ROIs to extract).
        lambda2         : Regularisation parameter (= 1/SNR²).
        pick_ori        : Orientation constraint: None, "normal", or "vector".
        roi_extract_mode: Extraction mode passed to extract_label_time_course
                          (e.g. "mean_flip", "mean", "pca_flip").
        logf            : Log file handle.

    Returns:
        tc : np.ndarray, shape (n_epochs, n_rois, n_times)
        times : np.ndarray, shape (n_times,) - epoch time axis in seconds
    """
    src_logmsg(logf, "[INV] Applying sLORETA to %s epochs...", len(epochs))
    stcs = mne.minimum_norm.apply_inverse_epochs(
        epochs,
        inv,
        lambda2 = lambda2,
        method = "sLORETA",
        pick_ori = pick_ori,
        return_generator = False,
        verbose = "ERROR",
    )
    src_logmsg(logf, "[ROI] Extracting %d ROI time courses (mode = '%s')...",
               len(labels), roi_extract_mode)
    tc = mne.extract_label_time_course(
        stcs,
        labels = labels,
        src = fwd["src"],
        mode = roi_extract_mode,
        verbose = "ERROR",
    )
    tc = np.ndarray(tc) # (n_epochs, n_rois, n_times)
    return tc, epochs.times

def src_apply_inverse_evoked(
        epochs: mne.Epochs,
        inv,
        lambda2: float,
        pick_ori: Optional[str],
        logf,
):
    """
    Apply sLORETA to the grand-average evoked response.

    Used for the brain-snapshot QC figure and identifying
    LEP peak latencies across subjects.

    Returns:
        stc_ga : MNE SourceTimeCourse for the grand average.
    """
    src_logmsg(logf, "[EVOKED] Applying sLORETA to grand-average evoked...")
    evoked = epochs.average()
    stc_ga = mne.minimum_norm.apply_inverse(
        evoked,
        inv,
        lambda2 = lambda2,
        method = "sLORETA",
        pick_ori = pick_ori,
        verbose = "ERROR",
    )
    return stc_ga