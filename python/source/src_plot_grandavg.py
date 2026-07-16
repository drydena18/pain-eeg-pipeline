"""
src_plot_grandavg.py  –  Cross-subject grand-average sLORETA heatmaps.
V 1.0.0

Loads per-subject GA timecourse .npy files written by source_core.py (step 6b)
and produces two heatmaps (lh, rh) of the grand-average sLORETA amplitude
across all available subjects and experiments.

Output
------
    <out_dir>/grandavg_source_GA_timecourse_lh.png
    <out_dir>/grandavg_source_GA_timecourse_rh.png

Usage
-----
    python plot_grandavg_source.py

Edit the CONFIG block below before running.
"""

from __future__ import annotations

import os
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# =============================================================================
# CONFIG — edit before running
# =============================================================================

DA_ROOT = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis"

EXPERIMENTS = [
    "26ByBiosemi",
    "142ByBiosemi",
    "29ByANT",
]

# Output directory for the grand-average figures
OUT_DIR = os.path.join(DA_ROOT, "grandavg_source_figures")

# Desikan-Killiany ROI names (aparc, 68 bilateral labels).
# These must match the roi_names order used when source_core.py was run.
# If you used use_custom_rois=false, this is the full parcellation order
# as loaded by MNE from fsaverage. The script will infer hemisphere split
# from the -lh / -rh suffix on each name.
# Leave as None to infer from the first subject's CSV roi column.
ROI_NAMES = None   # or provide an explicit list of strings

# Colormap limits (z-score)
VMIN, VMAX = -3, 3

# =============================================================================
# HELPERS
# =============================================================================

def _collect_npy_paths(da_root: str, experiments: list[str]) -> list[tuple[str, str]]:
    """
    Glob for all per-subject GA timecourse .npy files across experiments.

    Returns list of (sub_str, npy_path) tuples.
    """
    found = []
    for exp in experiments:
        pattern = os.path.join(
            da_root, exp, "source", "sub-*", "csv", "*_source_ga_timecourse.npy"
        )
        for path in sorted(glob.glob(pattern)):
            sub_str = os.path.basename(path).replace("_source_ga_timecourse.npy", "")
            found.append((sub_str, path))
    return found


def _load_times(npy_path: str) -> np.ndarray:
    """Load the times array saved alongside the GA timecourse."""
    times_path = npy_path.replace("_source_ga_timecourse.npy", "_source_times.npy")
    if not os.path.exists(times_path):
        raise FileNotFoundError(f"Times file not found: {times_path}")
    return np.load(times_path)


def _infer_roi_names(da_root: str, experiments: list[str]) -> list[str]:
    """
    Read roi names from the first available source_ga.csv.
    Falls back to raising an error if none found.
    """
    for exp in experiments:
        pattern = os.path.join(
            da_root, exp, "source", "sub-*", "csv", "*_source_ga.csv"
        )
        matches = sorted(glob.glob(pattern))
        if matches:
            import pandas as pd
            df = pd.read_csv(matches[0])
            if "roi" in df.columns:
                return df["roi"].tolist()
    raise RuntimeError(
        "Could not infer ROI names from any source_ga.csv. "
        "Set ROI_NAMES explicitly in the CONFIG block."
    )


def _zscore_rows(data: np.ndarray) -> np.ndarray:
    """Z-score each ROI row independently. Shape: (n_rois, n_times)."""
    out = np.zeros_like(data)
    for ri in range(data.shape[0]):
        row = data[ri, :]
        std = np.std(row)
        out[ri, :] = (row - np.mean(row)) / std if std > 0 else row
    return out


def _plot_heatmap(
        data_z: np.ndarray,
        roi_labels: list[str],
        t_ms: np.ndarray,
        hemi: str,
        n_subjects: int,
        out_path: str,
):
    """
    Save a single hemisphere heatmap.

    Args:
        data_z     : (n_rois, n_times) — already z-scored
        roi_labels : ROI name strings for y-axis (hemisphere suffix stripped)
        t_ms       : Time axis in milliseconds
        hemi       : 'lh' or 'rh'
        n_subjects : Number of subjects contributing (for title)
        out_path   : Full output .png path
    """
    n_roi = data_z.shape[0]
    fig_h = max(5, n_roi * 0.22)
    fig, ax = plt.subplots(figsize=(11, fig_h))

    im = ax.imshow(
        data_z,
        aspect="auto",
        origin="upper",
        extent=[t_ms[0], t_ms[-1], n_roi - 0.5, -0.5],
        cmap="RdBu_r",
        vmin=VMIN, vmax=VMAX,
    )

    ax.set_yticks(range(n_roi))
    ax.set_yticklabels(roi_labels, fontsize=7)
    ax.set_xlabel("Time (ms)", fontsize=9)
    ax.set_title(
        f"Grand-Average sLORETA  {hemi}  (z-scored per ROI)  N={n_subjects}",
        fontsize=10,
    )

    ax.axvline(0, color="k", lw=1.2, linestyle="--")

    legend_elements = [
        Line2D([0], [0], color="k", lw=1.2, linestyle="--", label="Stimulus onset"),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc="upper right", framealpha=0.7)

    plt.colorbar(im, ax=ax, label="Z-score", shrink=0.6, pad=0.01)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[FIG] {out_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Collect subject .npy paths ────────────────────────────────────────────
    subject_paths = _collect_npy_paths(DA_ROOT, EXPERIMENTS)
    if not subject_paths:
        raise RuntimeError(
            "No _source_ga_timecourse.npy files found. "
            "Make sure source_core.py step 6b has been run for all subjects."
        )
    print(f"[LOAD] Found {len(subject_paths)} subjects across {len(EXPERIMENTS)} experiments")

    # ── Infer ROI names if not provided ───────────────────────────────────────
    roi_names = ROI_NAMES or _infer_roi_names(DA_ROOT, EXPERIMENTS)
    print(f"[ROI] {len(roi_names)} ROIs")

    # ── Load and stack GA timecourses ─────────────────────────────────────────
    # Validate that all subjects share the same shape before stacking.
    arrays = []
    times  = None

    for sub_str, npy_path in subject_paths:
        try:
            ga_tc = np.load(npy_path)           # (n_rois, n_times)
            t     = _load_times(npy_path)

            if times is None:
                times = t
                expected_shape = ga_tc.shape
            else:
                if ga_tc.shape != expected_shape:
                    print(f"[WARN] {sub_str} shape {ga_tc.shape} != expected "
                          f"{expected_shape} — skipping")
                    continue

            arrays.append(ga_tc)
            print(f"[OK]  {sub_str}  shape={ga_tc.shape}")

        except Exception as e:
            print(f"[WARN] {sub_str} — could not load: {e}")

    if not arrays:
        raise RuntimeError("No valid subject arrays loaded.")

    n_subjects = len(arrays)
    stack = np.stack(arrays, axis=0)            # (n_subjects, n_rois, n_times)
    grand_avg = np.mean(stack, axis=0)          # (n_rois, n_times)
    t_ms = times * 1000.0

    print(f"[GA]  Grand-average shape: {grand_avg.shape}  N={n_subjects}")

    # ── Split by hemisphere and plot ──────────────────────────────────────────
    for hemi in ("lh", "rh"):
        idx = [i for i, n in enumerate(roi_names) if n.endswith(f"-{hemi}")]
        if not idx:
            print(f"[WARN] No ROIs found for {hemi} — skipping")
            continue

        data    = grand_avg[idx, :]
        labels  = [roi_names[i].replace(f"-{hemi}", "") for i in idx]
        data_z  = _zscore_rows(data)

        out_path = os.path.join(OUT_DIR, f"grandavg_source_GA_timecourse_{hemi}.png")
        _plot_heatmap(data_z, labels, t_ms, hemi, n_subjects, out_path)

    print("[DONE]")


if __name__ == "__main__":
    main()