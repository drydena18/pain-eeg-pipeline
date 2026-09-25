"""
src_spectral.py - Shared spectral computation primitives.
V2.0.0

These are low-level building blocks used by src_prestim.py, src_poststim.py,
and src_fooof.py. No domain-specific metric logic lives here.

V2.0.0 changes vs V1.0.0:
    - NEW src_band_power_tc: filter-Hilbert band-power time courses. Slow,
        fast, and total alpha power are now estimated by bandpass-filtering the
        FULL epoch, taking |analytic signal|^2, and averaging that power inside
        each analysis window. The frequency selectivity is set by the filer,
        which is identical for every window, instead of by the Welch segment
        length, which differed between the whole / pre / post windows and was
        too short (0.5s -> 1.95 Hz bins) to put more than one bin inside a
        2 Hz sub-band.
    - src_bandpower now returns NaN when fewer than 2 bins fall in the band.
        The trapezoidal rule over a single point has zero width and previously
        returned exactly 0.0, which silently zeroed every slow/fast metric.
    - src_psd_welch zero-pads to a target bin spacing (df_target, default
    0.25 Hz). This only interpolates the spectram (true resolution is still
    ~1 / segment length) but gives the centre-of-gravity and FOOOF fits
    enough bins to work with.

Covers:
    - Welch PSD estimation on a 1-D time course
    - Band-power integation (trapezoidal rule)
    - Spectral CoG - used as PAF proxy
    - Bandpass filtering (zero-phase FIR) for Hilbert-based phase extraction
"""

from __future__ import annotations

import numpy as np
import mne
from scipy.signal import hilbert

# ====================================================================
# WELCH PSD
# ====================================================================
def src_psd_welch(
        x: np.ndarray,
        sfreq: float,
        fmin: float,
        fmax: float,
        window_sec: 2.0,
        overlap: 0.5,
        df_target: float = 0.25,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute a Welch PSD estimate for a 1-D time series.

    The Welch segment length is min(window_sec, signal duration), so short
    windows (e.g., a 0.9s pre-stim window) become a single Hann-tapered
    segment rather than crashing. n_fft is zero-padded so the frequency grid
    spacing is at most df_target Hz.

    Args:
        x           : 1-D signal array, shape (n_times,)
        sfreq       : Sampling frequency in Hz
        fmin, fmax  : Frequency range to return
        window_sec  : Desired Welch segment length in secods (default 0.5 s)
        overlap     : Segment overlap fraction in [0, 1] (default 0.5)
        df_target   : Maximum frequency-grid spacing in Hz after zero-padding

    Returns:
        freqs : np.ndarray, shape (n_freqs,)
        psd   : np.ndarray, shape (n_freqs,)
    """
    n_times = len(x)
    n_win = min(int(window_sec * sfreq), n_times)
    n_win = max(n_win, 8) # minimum sanity bound
    n_step = max(1, int(n_win * (1.0 - overlap)))
    n_fft = max(int(2 ** np.ceil(np.log2(n_win))), int(np.ceil(sfreq / df_target)))

    psd, freqs = mne.time_frequency.psd_array_welch(
        x[np.newaxis, :],
        sfreq = sfreq,
        fmin = fmin,
        fmax = fmax,
        n_fft = n_fft,
        n_per_seg = n_win,
        n_overlap = n_win - n_step,
        average = "mean",
        verbose = "ERROR",
    )
    return freqs, psd.squeeze()

# ====================================================================
# BAND-POWER FROM A PSD
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
    if np.count_nonzero(idx) < 2:
        return float("nan")
    return float(np.trapezoid(psd[idx], freqs[idx]))

# ====================================================================
# FILTER-HILBERT BAND-POWER TIME COURSE
# ====================================================================
def src_band_power_tc(
        x: np.ndarray,
        sfreq: float,
        lo: float,
        hi: float,
        trans_bw: float = 1.5,
) -> np.ndarray:
    """
    Instantaneous band power |analytic signal|^2 of x in [lo, hi] Hz.

    The signal is bandpass filtered along the LAST axis with a zero-phase,
    Hamming-windowed FIT (firwin design) with explicit transition bandwidths,
    then Hilbert transformed. Filter the FULL epoch and crop afterwards:
    filtering a cropped window would put the filter's edge transitents inside
    the window you are measuring.

    Choice of trans_bw (Hz) trades frequency selectivity against temporal
    smearing. Filter length ~ 3.3 / trans_bw seconds:
        1.0 Hz -> 3.3 s (longer than a 3s epoch; relies heavily on padding)
        1.5 Hz -> 2.2 s (default; sharpest filter that fits inside the epoch)
        2.0 Hz -> 1.65 s
    With 8-10 / 10-12 Hz bands and trans_bw = 1.5, a pure 9 Hz sinusoid is
    ~95 % slow ad a pure 11 Hz sinusoid is ~95 % fast; power at exactly 10 Hz
    is split 50/50 by construction (it sits on the shared band edge).

    By Parseval, the time-average of |analytic|^2 over a window approximates
    the PSD integrated over the passband, so values are on the same scale as
    src_bandpower() on a PSD.

    Args:
        x           : array, shape (..., n_times)
        sfreq       : Sampling frequency in Hz
        lo, hi      : Passband edges in Hz
        trans_bw    : Transition bandwidth in Hz, applied to both edges

    Returns:
        power : array, same shape as x
    """
    shape = x.shape
    x2 = np.ascontiguousarray(x.reshape(-1, shape[-1]), dtype = np.float64)
    x_filt = mne.filter.filter_data(
        x2,
        sfreq = sfreq,
        l_freq = lo,
        h_freq = hi,
        l_trans_bandwidth = trans_bw,
        h_trans_bandwidth = trans_bw,
        method = "fir",
        fir_window = "hamming",
        fir_design = "firwin",
        phase = "zero",
        pad = "reflect_limited",
        verbose = "ERROR",
    )
    power = np.abs(hilbert(x_filt, axis = -1)) ** 2
    return power.reshape(shape)

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
    Apply a zero-phase FIR (finite impulse response) bandpass filter to a 1;D
    signal, suitable for subsequent Hilbert transform phase extraction.

    MNE's filter_data uses a Hamming-windowed FIR with automatic order
    selection to achieve the specified transition bandwidth, which avoids
    phase distortion problems that affect IIR (infinite impulse repsonse) designs.

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