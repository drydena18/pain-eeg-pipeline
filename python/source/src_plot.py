"""
src_plot.py - Plotting helpers for the source localization pipeline
V 2.5.0

Current stage: pre-stats. Two QC figures per subject (lh + rh heatmaps).

Per-ROI figures are deferred until after R GAMMs identify significant ROIs.

Figures produced
-----------------
    src_plot_ga_timecourse    ROI × time heatmap, one file per hemisphere:
                                  <prefix>_ga_timecourse_lh.png
                                  <prefix>_ga_timecourse_rh.png
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from src_io import src_logmsg


def src_plot_ga_timecourse(
        path_prefix: str,
        tc: np.ndarray,
        times: np.ndarray,
        roi_names: list[str],
        sub_str: str,
        n2_window: tuple[float, float],
        p2_window: tuple[float, float],
        logf,
):
    """
    Save two heatmaps of GA sLORETA amplitude — one per hemisphere.

    Each file is a clean ROI × time heatmap (z-scored per ROI) with
    stimulus onset, N2, and P2 windows marked.

    Args:
        path_prefix : Path without extension or hemisphere suffix, e.g.
                      '.../figures/sub-001_source_GA_timecourse'
        tc          : (n_epochs, n_rois, n_times)
        times       : (n_times,) in seconds
        roi_names   : Ordered ROI name strings
        sub_str     : e.g. "sub-001"
        n2_window   : (t_lo, t_hi) in seconds
        p2_window   : (t_lo, t_hi) in seconds
        logf        : Log file handle
    """
    ga   = np.mean(tc, axis=0)   # (n_rois, n_times)
    t_ms = times * 1000.0

    hemi_indices = {
        "lh": sorted([i for i, n in enumerate(roi_names) if n.endswith("-lh")]),
        "rh": sorted([i for i, n in enumerate(roi_names) if n.endswith("-rh")]),
    }

    legend_elements = [
        Line2D([0], [0], color="k",        lw=1.2, linestyle="--", label="Stimulus"),
        Line2D([0], [0], color="steelblue", lw=1.0, linestyle=":",  label="N2 window"),
        Line2D([0], [0], color="tomato",    lw=1.0, linestyle=":",  label="P2 window"),
    ]

    for hemi, idx in hemi_indices.items():
        if not idx:
            continue
        path_png = f"{path_prefix}_{hemi}.png"
        try:
            data  = ga[idx, :]
            names = [roi_names[i].replace(f"-{hemi}", "") for i in idx]
            n_roi = len(idx)

            # Z-score each ROI row independently
            data_z = np.zeros_like(data)
            for ri in range(n_roi):
                row = data[ri, :]
                std = np.std(row)
                data_z[ri, :] = (row - np.mean(row)) / std if std > 0 else row

            fig_h = max(5, n_roi * 0.22)
            fig, ax = plt.subplots(figsize=(11, fig_h))

            im = ax.imshow(
                data_z,
                aspect="auto",
                origin="upper",
                extent=[t_ms[0], t_ms[-1], n_roi - 0.5, -0.5],
                cmap="RdBu_r",
                vmin=-3, vmax=3,
            )

            ax.set_yticks(range(n_roi))
            ax.set_yticklabels(names, fontsize=7)
            ax.set_xlabel("Time (ms)", fontsize=9)
            ax.set_title(
                f"{sub_str}  GA sLORETA  {hemi}  (z-scored per ROI)",
                fontsize=10,
            )

            ax.axvline(0,                   color="k",         lw=1.2, linestyle="--")
            ax.axvline(n2_window[0] * 1000, color="steelblue", lw=1.0, linestyle=":")
            ax.axvline(n2_window[1] * 1000, color="steelblue", lw=1.0, linestyle=":")
            ax.axvline(p2_window[0] * 1000, color="tomato",    lw=1.0, linestyle=":")
            ax.axvline(p2_window[1] * 1000, color="tomato",    lw=1.0, linestyle=":")

            ax.legend(handles=legend_elements, fontsize=7,
                      loc="upper right", framealpha=0.7)

            plt.colorbar(im, ax=ax, label="Z-score", shrink=0.6, pad=0.01)

            fig.tight_layout()
            fig.savefig(path_png, dpi=150, bbox_inches="tight")
            plt.close(fig)
            src_logmsg(logf, "[FIG] %s", path_png)

        except Exception as e:
            src_logmsg(logf, "[WARN] Heatmap failed (%s): %s", hemi, str(e))