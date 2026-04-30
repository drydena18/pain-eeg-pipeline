"""
src_plot.py - Plotting helpers for the source localization pipeline

All functions:
    - Use the "Agg" matplotlib backend (no display required)
    - Accept a logf handle and log the saved path
    - Catch exceptions internally and warn rather than crash te subject loop
    - Produce 150 dpi PNGs

Figures Produced
-----------------
    src_plot_ga_psd         GA pre-stim PSD overlay for all ROIs
    src_plot_fooof          Single FOOOF model fit
    src_plot_erd            ERD_slow and ERD_fast time series across trials
    src_plot_lep_ga         Grand-average LEP waveform per ROI
    src_plot_phase_polar    Polar histogram of Hilbert phase at stimulus onset
    src_plot_brain          2-D lateral-view brain snapshot (sLORETA GA map)
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src_io import src_logmsg

# ====================================================================
# GA PRE-STIM PSD OVERLAY
# ====================================================================
def src_plot_ga_psd(
        path_png: str,
        psd_by_roi_idx: dict,
        roi_names: list[str],
        sub_str: str,
        fmin: float,
        fmax: float,
        logf,
):
    """
    Overlay grand-average pre-stimulus PSD curves for all ROIs on a single
    semi-log (dB) figure. One curve per ROI, coloured by a qualitative palette
    """
    try:
        fig, ax = plt.subplots(figsize = (10, 5))
        for ri, roi in enumerate(roi_names):
            freqs, psd = psd_by_roi_idx[ri]
            ax.plot(freqs, 10.0 * np.log10(psd * 1e-20), label = roi, linewidth = 1.3)

        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Power (dB)")
        ax.set_xlim(fmin, fmax)
        ax.axvspan(8, 10, color = "blue", alpha = 0.06, label = "Slow α")
        ax.axvspan(10, 12, color = "red", alpha = 0.06, label = "Fast α")
        ax.set_title(f"{sub_str} GA Pre-stimulus ROI PSD (sLORETA)")
        ax.legend(fontsize = 7, ncol = max(1, len(roi_names) // 8))
        fig.tight_layout()
        fig.savefig(path_png, dpi = 150)
        plt.close(fig)
        src_logmsg(logf, "[FIG] %s", path_png)
    except Exception as e:
        src_logmsg(logf, "[WARN] GA PSD plot failed: %s", str(e))

# ====================================================================
# FOOOF MODEL FIT
# ====================================================================
def src_plot_fooof(
        path_png: str,
        fm,
        title: str,
        logf,
):
    """
    Save a FOOOF model fit plot (aperiodic + peaks overlay).
    Works with both the legacy fooof and the new specparam API
    """
    try:
        # return_fig is supported by specparam; fooof uses a slightly different API
        try:
            fig = fm.plot(plot_peaks = "shade", add_legend = True, return_fig = True)
        except TypeError:
            fm.plot(plot_peaks = "shade", add_legend = True)
            fig = plt.gcf()

        fig.axes[0].set_title(title)
        fig.tight_layout()
        fig.savefig(path_png, dpi = 150)
        plt.close(fig)
        src_logmsg(logf, "[FIG] %s", path_png)
    except Exception as e:
        src_logmsg(logf, "[WARN] FOOOF plot failed (%s): %s", title, str(e))

# ====================================================================
# ERD TRIAL SERIES
# ====================================================================
def src_plot_erd(
        path_png: str,
        poststim_rows: list[dict],
        roi_names: list[str],
        sub_str: str,
        logf,
):
    """
    Plot the per-trial ERD_slow and ERD_fast series for each ROI in a grid.
    Each subplot shows one ROI's slow (blue) and fast (red) ERD across trials,
    making drift and outlier trials immediately visible
    """
    try:
        import pandas as pd
        df = pd.DataFrame(poststim_rows)
        df["roi"] = df["roi_idx"].map(lambda ri: roi_names[ri] if ri < len(roi_names) else str(ri))

        n_rois = len(roi_names)
        ncols = min(4, n_rois)
        nrows = int(np.ceil(n_rois / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize = (4 * ncols, 4 * nrows), squeeze = False)
        fig.suptitle(f"{sub_str} Per-trial ERD (sLORETA)", fontsize = 10)

        for ri, roi in enumerate(roi_names):
            ax = axes[ri // ncols][ri % ncols]
            sub_df = df[df["roi"] == roi].sort_values("trial")
            ax.plot(sub_df["trial"], sub_df["ERD_slow"], color = "steelblue", 
                    linewidth = 0.8, label = "slow")
            ax.plot(sub_df["trial"], sub_df["ERD_fast"], color = "tomato",
                    linewidth = 0.8, label = "fast")
            ax.axhline(0, color = "k", linewidth = 0.5, linestyle = "--")
            ax.set_title(roi, fontsize = 8)
            ax.set_xlabel("Trial", fontsize = 7)
            ax.set_ylabel("ERD", fontsize = 7)
            if ri == 0:
                ax.legend(fontsize = 6)

        # Hide unused subplots
        for ri in range(n_rois, nrows * ncols):
            axes[ri // ncols][ri % ncols].set_visible(False)

        fig.tight_layout()
        fig.savefig(path_png, dpi = 150)
        plt.close(fig)
        src_logmsg(logf, "[FID] %s", path_png)
    except Exception as e:
        src_logmsg(logf, "[WARN] ERD plot failed: %s", str(e))

# ====================================================================
# GRAND-AVERAGE LEP WAVEFORM
# ====================================================================
def src_plot_lep_ga(
        path_png: str,
        tc_post: np.ndarray,
        times_post: np.ndarray,
        roi_names: list[str],
        n2_window: tuple[float, float],
        p2_window: tuple[float, float],
        sub_str: str,
        logf,
):
    """
    Plot the grand-average LEP waveform for each ROI
    N2 and P2 search windows are shaded in blue and red respectively
    """
    try:
        n_rois = len(roi_names)
        ncols = min(4, n_rois)
        nrows = int(np.ceil(n_rois / ncols))
        t_ms = times_post * 1000.0

        fig, axes = plt.subplots(nrows, ncols, figsize = (4 * ncols, 3 * nrows), squeeze = False)

        for ri, roi in enumerate(roi_names):
            ax = axes[ri // ncols][ri % ncols]
            x_ga = np.mean(tc_post[:, ri, :], axis = 0)
            ax.plot(t_ms, x_ga, color = "black", linewidth = 0.5)
            ax.axvspan(n2_window[0] * 1000, n2_window[1] * 1000, 
                       color = "steelblue", alpha = 0.12, label = "N2 window")
            ax.axvspan(p2_window[0] * 1000, p2_window[1] * 1000,
                       color = "tomato", alpha = 0.12, label = "P2 window")
            ax.set_title(roi, fontsize = 8)
            ax.set_xlabel("Time (ms)", fontsize = 7)
            ax.set_ylabel("Source amp.", fontsize = 7)
            if ri == 0:
                ax.legend(fontsize = 6)

        for ri in range(n_rois, nrows * ncols):
            axes[ri // ncols][ri % ncols].set_visible(False)

        fig.tight_layout()
        fig.savefig(path_png, dpi = 150)
        plt.close(fig)
        src_logmsg(logf, "[FIG] %s", path_png)
    except Exception as e:
        src_logmsg(logf, "[WARN] LEP GA plot failed: %s", str(e))

# ====================================================================
# HILBERT PHASE POLAR HISTOGRAM
# ====================================================================
def src_plot_phase_polar(
        path_png: str,
        trial_rows: list[dict],
        roi_names: list[str],
        phase_col: str,
        sub_str: str,
        title_suffix: str,
        logf,
):
    """
    Plot a polar histogram of instantaneous Hilbert phase for each ROI

    The uniform distribution (grey ring) and the mean resultant vector (red
    arrow) are overlaid so that phase clustering (or lack thereof) is
    immediately visible. Works for both slow_phase (pre-stim) and
    slow_phase_post (post-stim)

    Args:
        trial_rows      : List of per-trial dicts containing 'roi_idx' and phase_col
        phase_col       : Which phase field to plot (e.g., "slow_phase", "slow_phase_post")
        title_suffix    : Short string appended to the suptitle (e.g., "pre-stim onset")
    """
    try:
        import pandas as pd
        df = pd.DataFrame(trial_rows)
        df["roi"] = df["roi_idx"].map(lambda ri: roi_names[ri] if ri < len(roi_names) else str(ri))

        if phase_col not in df.columns:
            src_logmsg(logf, "[WARN] Phase column '%s' not in trial rows - skipping polar plot.", phase_col)
            return
        
        n_rois = len(roi_names)
        ncols = min(4, n_rois)
        nrows = int(np.ceil(n_rois / ncols))

        fig = plt.figure(figsize = (4 * ncols, 3.5 * nrows))
        fig.suptitle(f"{sub_str} Phase distribution ({title_suffix})", fontsize = 10)

        for ri, roi in enumerate(roi_names):
            ax = fig.add_subplot(nrows, ncols, ri + 1, projection = "polar")
            phases = df.loc[df["roi"] == roi, phase_col].dropna().values

            if len(phases) < 4:
                ax.set_title(roi + "\n(insufficient data)", fontsize = 7)
                continue

            n_bins = 18
            counts, edges = np.histogram(phases, bins = n_bins, range = (-np.pi, np.pi))
            width = edges[1] - edges[0]
            centres = 0.5 * (edges[:-1] + edges[1:])
            ax.bar(centres, counts, width = width, color = "steelblue",
                   alpha = 0.7, edgecolor = "white", linewidth = 0.4)
            
            # Mean resultant vector (phase clustering)
            r_mean = np.abs(np.mean(np.exp(1j * phases)))
            theta_mu = np.angle(np.mean(np.exp(1j * phases)))
            ax.annotate("", xy = (theta_mu, r_mean * counts.max()),
                        xytext = (0, 0),
                        arrowprops = dict(arrowstyle = "->", color = "red", lw = 1.5))
            
            ax.set_title(f"{roi}\nR = {r_mean:.2f}", fontsize = 7)
            ax.set_yticks([])

        fig.tight_layout()
        fig.savefig(path_png, dpi = 150)
        plt.close(fig)
        src_logmsg(logf, "[FIG] %s", path_png)
    except Exception as e:
        src_logmsg(logf, "[WARN] Phase polar plot failed: %s", str(e))


# ====================================================================
# BRAIN SNAPSHOT (sLORETA GA MAP)
# ====================================================================
def src_plot_brain(
        path_png: str,
        stc_ga,
        subjects_dir: str,
        t_sec: float,
        sub_str: str,
        logf,
):
    """
    Save a 2-D lateral-view brain snapshot of the GA sLORETA activation map
    at a given time point using MNE's matplotlib backend

    Args:
        stc_ga          : MNE SourceTimeEstimate from src_apply_inverse_evoked()
        subjects_dir    : Directory containing the fsaverage/ folder
        t_sec           : Time in seconds at which to snapshot the map
        sub_str         : e.g., "sub-001" used in the figure title
    """
    try:
        brain = stc_ga.plot(
            subject = "fsaverage",
            subjects_dir = subjects_dir,
            initial_time = t_sec,
            hemi = "split",
            views = "lat",
            backend = "matplotlib",
            time_viewer = False,
            colorbar = True,
            size = (900, 400),
            show = False,
        )
        brain.suptitle(f"{sub_str} GA sLORETA @ {t_sec:.3f} s")
        brain.savefig(path_png, dpi = 150)
        plt.close(brain)
        src_logmsg(logf, "[FIG] %s", path_png)
    except Exception as e:
        src_logmsg(logf, "[WARN] Brain snapshot failed: %s", str(e))