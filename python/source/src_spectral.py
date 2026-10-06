"""
src_spectral.py - Shared spectral computation primitives.
V 3.0.0

These are the low-level estimators used by src_alpha_features.py,
src_erd.py, src_prestim.py, src_poststim.py and src_fooof.py. No
domain-specific metric logic lives here.

PARITY CONTRACT (V3.0.0)
------------------------
Every function in this file has a line-for-line MATLAB twin in
matlab/spectral/helpers/ so that channel-space (MATLAB) and source-space
(Python) features are computed by the SAME estimator:

    src_window_mask        <->  spec_window_mask.m
    src_fir_bandpass       <->  spec_fir_bandpass.m
    src_zero_phase_filter  <->  spec_zero_phase_filter.m
    src_analytic_signal    <->  spec_analytic_signal.m
    src_band_power_tc      <->  spec_band_power_tcs.m
    src_psd_welch          <->  spec_welch_psd.m
    src_cog_from_psd       <->  spec_cog_from_psd.m

The FIR filter and Welch PSD are written out explicitly (no mne.filter, no
mne.time_frequency) because library defaults (filter length rules, window
symmetry, detrending, padding) differ between MNE/SciPy and EEGLAB/MATLAB.
tests/parity_test.py checks the two implementations agree numerically.

V3.0.0 changes vs V2.0.0:
    - Band-pass filter is now an explicit windowed-sinc FIR:
        N = floor(3.3 * sfreq / trans_bw + 0.5), forced odd (Hamming rule)
        h = (lowpass(hi) - lowpass(lo)) * symmetric Hamming(N)
        normalised to unit gain at the band centre.
      The -6 dB points sit EXACTLY on the band edges (8 / 10 / 12 Hz), so the
      slow (8-10) and fast (10-12) bands meet at 10 Hz instead of overlapping
      by a full transition width as MNE's default design did.
    - Zero-phase application: reflect-pad by min(N, n_times - 1) samples,
      FFT convolution with the symmetric (linear-phase) kernel, crop.
    - Band power is now |analytic|^2 / 2, i.e. the variance of the
      band-limited signal. By Parseval this equals the one-sided PSD
      integrated over the passband (V2.0.0 was 2x larger and its docstring
      wrongly claimed the scales matched).
    - Welch PSD uses scipy.signal.welch with explicit settings: periodic
      Hamming window, no detrending, density scaling, zero-padding to
      df_target. Vectorised over leading axes.
    - Window masks are inclusive with a quarter-sample tolerance so that a
      boundary sample is never dropped by floating-point noise in the time
      axis (EEGLAB and MNE build the axis differently).
    - src_bandpower kept (NaN for < 2 bins); not used by the pipeline.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import fftconvolve, welch

# Hamming-window length factor: transition bandwidth ~ 3.3 * sfreq / N.
FIR_LENGTH_FACTOR = 3.3


# ====================================================================
# TIME-WINDOW MASK
# ====================================================================
def src_window_mask(times: np.ndarray, tmin: float, tmax: float, sfreq: float) -> np.ndarray:
    """Inclusive [tmin, tmax] mask with a quarter-sample tolerance."""
    tol = 0.25 / float(sfreq)
    times = np.asarray(times, dtype = float)
    return (times >= tmin - tol) & (times <= tmax + tol)


# ====================================================================
# FIR BAND-PASS DESIGN (shared definition, see module docstring)
# ====================================================================
def src_fir_bandpass(sfreq: float, lo: float, hi: float, trans_bw: float) -> np.ndarray:
    """
    Linear-phase windowed-sinc band-pass FIR.

    Args:
        sfreq    : Sampling frequency in Hz
        lo, hi   : -6 dB band edges in Hz (0 < lo < hi < sfreq / 2)
        trans_bw : Transition bandwidth in Hz, centred on each edge

    Returns:
        h : (N,) filter taps, N odd, unit gain at (lo + hi) / 2
    """
    if not (0.0 < lo < hi < sfreq / 2.0):
        raise ValueError(f"Bad band [{lo}, {hi}] Hz for sfreq = {sfreq} Hz.")
    n = int(np.floor(FIR_LENGTH_FACTOR * sfreq / trans_bw + 0.5))   # round half up, as MATLAB
    if n % 2 == 0:
        n += 1
    m = (n - 1) // 2
    k = np.arange(n) - m

    def _lowpass(fc: float) -> np.ndarray:
        w = 2.0 * fc / sfreq
        return w * np.sinc(w * k)          # np.sinc(x) = sin(pi x) / (pi x)

    win = 0.54 - 0.46 * np.cos(2.0 * np.pi * np.arange(n) / (n - 1))   # symmetric Hamming
    h = (_lowpass(hi) - _lowpass(lo)) * win

    fc = 0.5 * (lo + hi)
    gain = np.abs(np.sum(h * np.exp(-2j * np.pi * fc * np.arange(n) / sfreq)))
    return h / gain


# ====================================================================
# ZERO-PHASE APPLICATION ALONG THE LAST AXIS
# ====================================================================
def src_zero_phase_filter(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    """
    Apply a symmetric odd-length FIR along the last axis with zero phase.

    Reflect-pads (mirror, edge sample not repeated) by min(N, n_times - 1)
    samples on each side, convolves ('same', centred), then crops.
    """
    x = np.asarray(x, dtype = np.float64)
    n_t = x.shape[-1]
    pad = min(len(h), n_t - 1)
    pad_width = [(0, 0)] * (x.ndim - 1) + [(pad, pad)]
    xp = np.pad(x, pad_width, mode = "reflect")
    kernel = np.asarray(h, dtype = np.float64).reshape((1,) * (x.ndim - 1) + (-1,))
    y = fftconvolve(xp, kernel, mode = "same", axes = -1)
    return y[..., pad:pad + n_t]


# ====================================================================
# ANALYTIC SIGNAL (same definition as MATLAB hilbert / scipy hilbert)
# ====================================================================
def src_analytic_signal(x: np.ndarray) -> np.ndarray:
    """FFT-based analytic signal along the last axis."""
    x = np.asarray(x, dtype = np.float64)
    n = x.shape[-1]
    spec = np.fft.fft(x, axis = -1)
    g = np.zeros(n)
    g[0] = 1.0
    if n % 2 == 0:
        g[1:n // 2] = 2.0
        g[n // 2] = 1.0
    else:
        g[1:(n + 1) // 2] = 2.0
    return np.fft.ifft(spec * g, axis = -1)


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
    Instantaneous band power |analytic(x_band)|^2 / 2 along the last axis.

    Filter the FULL epoch and crop afterwards so the filter's edge
    transients fall outside the analysed windows. With trans_bw = 1.5 Hz at
    500 Hz the kernel is 1101 taps (2.2 s), which fits inside the 3 s epoch.

    The window mean of the returned power equals the variance of the
    band-limited signal, i.e. the one-sided PSD integrated over [lo, hi].

    Args:
        x        : array (..., n_times), e.g. (n_epochs, n_rois, n_times)
        sfreq    : Sampling frequency in Hz
        lo, hi   : -6 dB band edges in Hz
        trans_bw : Transition bandwidth in Hz

    Returns:
        power : array, same shape as x
    """
    h = src_fir_bandpass(sfreq, lo, hi, trans_bw)
    if len(h) > x.shape[-1]:
        raise ValueError(
            f"FIR length {len(h)} samples exceeds the epoch length {x.shape[-1]}; "
            "increase features.band_power.trans_bw_hz or lengthen the epoch."
        )
    z = src_analytic_signal(src_zero_phase_filter(x, h))
    return (np.abs(z) ** 2) / 2.0


# ====================================================================
# BAND-PASS FILTER (for Hilbert phase / ITC) - same FIR as band power
# ====================================================================
def src_bandpass_filter(
        x: np.ndarray,
        sfreq: float,
        lo: float,
        hi: float,
        trans_bw: float = 1.5,
) -> np.ndarray:
    """Zero-phase band-pass along the last axis using src_fir_bandpass."""
    return src_zero_phase_filter(x, src_fir_bandpass(sfreq, lo, hi, trans_bw))


# ====================================================================
# WELCH PSD
# ====================================================================
def src_psd_welch(
        x: np.ndarray,
        sfreq: float,
        fmin: float,
        fmax: float,
        window_sec: float = 2.0,
        overlap: float = 0.5,
        df_target: float = 0.25,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Welch PSD along the last axis.

    Segment length = min(window_sec * sfreq, n_times) samples (>= 8), so a
    0.9 s window becomes one Hamming-tapered segment. n_fft is zero-padded
    so the grid spacing is at most df_target Hz (true resolution is still
    ~1 / segment length). Periodic Hamming, no detrending, density scaling,
    one-sided, mean across segments.

    Returns:
        freqs : (n_freqs,)
        psd   : (..., n_freqs)
    """
    x = np.asarray(x, dtype = np.float64)
    n_times = x.shape[-1]
    n_win = max(min(int(np.floor(window_sec * sfreq + 1e-9)), n_times), 8)
    n_step = max(1, int(np.floor(n_win * (1.0 - overlap) + 1e-9)))
    n_fft = max(int(2 ** np.ceil(np.log2(n_win))), int(np.ceil(sfreq / df_target - 1e-9)))

    freqs, psd = welch(
        x,
        fs = sfreq,
        window = "hamming",
        nperseg = n_win,
        noverlap = n_win - n_step,
        nfft = n_fft,
        detrend = False,
        return_onesided = True,
        scaling = "density",
        axis = -1,
        average = "mean",
    )
    keep = (freqs >= fmin - 1e-9) & (freqs <= fmax + 1e-9)
    return freqs[keep], psd[..., keep]


# ====================================================================
# SPECTRAL CENTRE OF GRAVITY
# ====================================================================
def src_cog_from_psd(freqs: np.ndarray, psd: np.ndarray, band: tuple[float, float]) -> np.ndarray:
    """
    Power-weighted mean frequency over [band] along the last axis.
    NaN when fewer than 2 bins fall in the band or band power <= 0.
    """
    idx = (freqs >= band[0] - 1e-9) & (freqs <= band[1] + 1e-9)
    psd = np.asarray(psd, dtype = float)
    if np.count_nonzero(idx) < 2:
        return np.full(psd.shape[:-1], np.nan)
    fb = freqs[idx]
    pb = psd[..., idx]
    den = np.trapezoid(pb, fb, axis = -1)
    num = np.trapezoid(pb * fb, fb, axis = -1)
    with np.errstate(invalid = "ignore", divide = "ignore"):
        cog = num / den
    return np.where(np.isfinite(den) & (den > 0), cog, np.nan)


# ====================================================================
# BAND-POWER FROM A PSD (kept for reference; not used by the pipeline)
# ====================================================================
def src_bandpower(freqs: np.ndarray, psd: np.ndarray, lo: float, hi: float) -> float:
    """Trapezoidal PSD integral over [lo, hi]; NaN if fewer than 2 bins."""
    idx = (freqs >= lo) & (freqs <= hi)
    if np.count_nonzero(idx) < 2:
        return float("nan")
    return float(np.trapezoid(psd[idx], freqs[idx]))
