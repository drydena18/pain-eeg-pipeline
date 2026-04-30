"""
src_spectral.py - Shared spectral computation primitives.

These are low-level building blocks used by src_prestim.py, src_poststim.py,
and src_fooof.py. No domain-specific metric logic lives here.

Covers:
    - Welch PSD estimation on a 1-D time course
    - Band-power integation (trapezoidal rule)
    - Spectral CoG - used as PAF proxy
    - Bandpass filtering (zero-phase FIR) for Hilbert-based phase extraction
"""

from __future__ import annotations

import numpy as np
import mne

# ====================================================================
# WELCH PSD
# ====================================================================
def src_psd_welch(
        x: np.ndarray,
        sfreq: float,
        fmin: float,
        fmax: float,
        window_sec: float = 0.5,
        overlap: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute a Welch PSD estimate for a 1-D time series.

    The Welch window length is chosen as the minimum of window_sec and the 
    full signal duration, so short pre-stimulus windows (e.g., 500 ms) are
    handled gracefully without crashing.

    Args:
        x           : 1-D signal array, shape (n_times,)
        sfreq       : Sampling frequency in Hz
        fmin, fmax  : Frequency range to return
        window_sec  : Desired Welch segment length in secods (default 0.5 s)
        overlap     : Segment overlap fraction in [0, 1] (default 0.5)

    Returns:
        freqs : np.ndarray, shape (n_freqs,)
        psd   : np.ndarray, shape (n_freqs,)
    """
    n_times = len(x)
    n_win = min(int(window_sec * sfreq), n_times)
    n_win = max(n_win, 8) # minimum sanity bound
    n_step = max(1, int(n_win * (1.0 - overlap)))

    psd, freqs = mne.time_frequency.psd_array_welch(
        x[np.newaxis, :],
        sfreq = sfreq,
        fmin = fmin,
        fmax = fmax,
        n_fft = max(int(2 ** np.ceil(np.log2(n_win))), n_win),
        n_per_seg = n_win,
        n_overlap = n_win - n_step,
        average = "mean",
        verbose = "ERROR",
    )
    return freqs, psd.squeeze()

# ====================================================================
# BAND-POWER
# ====================================================================
def src_bandpower(
        freqs: np.ndarray,
        psd: np.ndarray,
        lo: float,
        hi: float,
) -> float:
    """
    Integrate PSD over the band [lo, hi] Hz using the trapezoidal rule.

    Returns NaN if no frequency bins fall within [lo, hi].
    """
    idx = (freqs >= lo) & (freqs <= hi)
    if not np.any(idx):
        return float("nan")
    return float(np.trapezoid(psd[idx], freqs[idx]))

# ====================================================================
# CENTRE OF GRAVITY (peak alpha frequency proxy)
# ====================================================================
def src_cog(
        freqs: np.ndarray,
        psd: np.ndarray,
        lo: float,
        hi: float,
) -> float:
    """
    Compute the spectral centre of gravoty (CoG) over the band [lo, hi] Hz.

    CoG is used as a robust PAF proxy that is well-defined even when the
    alpha peak is broad or multimodal.

    Returns NaN if the band is empty or total power is zero.
    """
    idx = (freqs >= lo) & (freqs <= hi)
    f, p = freqs[idx], psd[idx]
    denom = np.sum(p)
    if denom <= 0 or len(f) == 0:
        return float("nan")
    return float(np.sum(f * p) / denom)

# ====================================================================
# BANDPASS FILTER (for Hilbert-based phase)
# ====================================================================
def src_bandpass_filter(
        x: np.ndarray,
        sfreq: float,
        lo: float,
        hi: float,
) -> np.ndarray:
    """
    Apply a zero-phase FIR (finite impulse response) bandpass filter to a 1-D
    signal, suitable for subsequent Hilbert transform phase extraction.

    MNE's filter_data uses a Hamming-windowed FIT with automatic order
    selection to achieve the specified transition bandwidth, which avoids
    phase distortion problems that affect IIR (infinite impulse response) desgins.

    Args:
        x       : 1-D signal array, shape (n_times,)
        sfreq   : Sampling frequency in Hz
        lo      : Lower passband edge in Hz
        hi      : Upper passband edge in Hz

    Returns:
        x_filt : Filtered 1-D array, same shape as x
    """
    return mne.filter.filter_data(
        x[np.newaxis, :].astype(np.float64),
        sfreq = sfreq,
        l_freq = lo,
        h_freq = hi,
        method = "fir",
        fir_window = "hamming",
        verbose = "ERROR",
    ).squeeze()