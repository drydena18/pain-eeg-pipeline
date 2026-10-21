"""
src_inverse.py - MNE forward solution and sLORETA inverse pipeline.

Covers:
    - Building (or loading from cache) the forward solution
    - Computing the noise covariance from a baseline window
    - Assembling the inverse operator
    - Applying sLORETA to all epochs -> per-trial source time course
    - Applying sLORETA to the grand-average evoked response
    - Extracting ROI time courses from source estimates

Noise Covariance Note (Oct 2026):
    Default is mne.make_ad_hoc_cov() (diagonal, MNE default EEG std 0.2 uV).
    The features of interest are ongoing pre-stimulus alpha, and an emptical
    covariance from the pre-stimulus baseline treats that alpha as noise:
    sLORETA whitens and noise-normalized with it, attenuating the sources that
    dominate the baseline. With no noise-only EEG recording available, a
    non-informative diagonal model is the standard choice for ongoing
    activity (Brainstorm noise-covariance tutorial; MNE developer guidance).
    For EEG-only data a uniform diagonal is a scaled identity: the std only
    rescales sLORETA output by a global 1/std, so spatial patterns and ratio
    metrics are unchanged and absolute power is in arbitraty but identical
    units for every subject/ Keep the default std fixed across all experiments
    so pre / post / ERD contrasts stay comparible across caps. "expirical"
    (baseline window) remains available as a sensitivity check.

Memory Note (Oct 2026):
    src_apply_inverse_epochs streams epochs through apply_inverse_epochs
    (return_generator = True) straight into extract_label_time_course, so only
    one vertex-level SourceEstimate exists at a time.
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
        noise_method: str = "ad_hoc",
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
        noise_method   : "ad_hoc" (default) -> mne.make_ad_hoc_cov(epochs.info),
                         diagonal, identical for every subject:
                         "empirical" -> mne.compute_covariance over
                         [noise_tmin, noise_tmax] (sensitivity check only)

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
    noise_method = str(noise_method).lower()
    if noise_method == "ad_hoc":
        src_logmsg(logf, "[COV] Ad hoc diagonal noise covariance (MNE default EEG std).")
        noise_cov = mne.make_ad_hoc_cov(epochs.info, verbose = "ERROR")
    elif noise_method == "empirical":
        src_logmsg(logf, "[COV] Empirical noise covariance from baseline (%.3f - %.3fs)...", noise_tmin, noise_tmax)
        noise_cov = mne.compute_covariance(
            epochs,
            tmin = noise_tmin,
            tmax = noise_tmax,
            method = "empirical",
            rank = None,
            verbose = "ERROR",
        )
    else:
        raise ValueError(f"noise_method = '{noise_method}' is not supported ('ad_hoc' or 'empirical').")

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

    Epochs are streamed: apply_inverse_epochs yields one SourceEstimate at a
    time and extract_label_time_course rreduces it to ROI time courses before
    the next epoch is inverted. Peak memory is one vertex-level estimate
    (~0.25 GB on ico-5, ~1 GV transient for loose orientation) plus the
    (n_epochs, n_rois, n_times) output, instead of all epochs at once.

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
    # -- Orientation Guard -----
    if pick_ori == "vector":
        raise ValueError(
            "src_apply_inverse_epochs: pick_ori = 'vector' gives (n_rois, 3, n_times) "
            "per epoch; the feature pipeline expects one signed time course per ROI. "
            "Use pick_ori = 'normal'."
        )
    if pick_ori is None:
        src_logmsg(logf, "[WARN] pick_ori = None returns non-negative vector magnitude; "
                         "band power and LEP polarity will be invalid. Use 'normal'.")

    src_logmsg(logf, "[INV] Applying sLORETA to %s epochs (streamed)...", len(epochs))
    stcs = mne.minimum_norm.apply_inverse_epochs(
        epochs,
        inv,
        lambda2 = lambda2,
        method = "sLORETA",
        pick_ori = pick_ori,
        return_generator = True,
        verbose = "ERROR",
    )
    src_logmsg(logf, "[ROI] Extracting %d ROI time courses (mode = '%s')...",
                   len(labels), roi_extract_mode)
    tc = mne.extract_label_time_course(
        stcs,
        labels = labels,
        src = fwd["src"],
        mode = roi_extract_mode,
        return_generator = False,
        verbose = "ERROR",
    )
    tc = np.stack(tc, axis = 0)

    if tc.shape[0] != len(epochs):
        raise RuntimeError(
            f"src_apply_inverse_epochs: got {tc.shape[0]} ROI time courses "
            f"for {len(epochs)} epochs."
        )
    src_logmsg(logf, "[ROI] Time courses ready: %d epochs x %d ROIs x %d samples (%.1f MB).",
                   tc.shape[0], tc.shape[1], tc.shape[2], tc.nbytes / 1e6)
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
        stc_ga : MNE SourceTimeCourse for the grand average
    """
    src_logmsg(logf, "[EVOKED] Applying sLORETA to grand-average evoked...")
    evoked = epochs.average()
    evoked = mne.minimum_norm.apply_inverse(
        evoked,
        inv,
        lambda2 = lambda2,
        method = "sLORETA",
        pick_ori = pick_ori,
        verbose = "ERROR",
    )
    return stc_ga