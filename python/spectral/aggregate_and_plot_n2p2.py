"""
aggregate_and_plot_n2p2.py  –  Load per-subject checkpoints from
run_one_subject.py, compute the grand-average signed ERP, and plot.
V 1.0.0

Run this after driver_run_all_subjects.py has completed (or partially
completed — it will use whatever checkpoints exist and report N).

Output
------
    <out_dir>/grandavg_n2p2_erp_signed_postcentral.png
    <out_dir>/grandavg_n2p2_erp_signed_supramarginal.png

Usage
-----
    python aggregate_and_plot_n2p2.py
"""

from __future__ import annotations

import os
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# CONFIG
# =============================================================================

DA_ROOT = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis"
CHECKPOINT_DIR = os.path.join(DA_ROOT, "n2p2_signed_checkpoints")
OUT_DIR = os.path.join(DA_ROOT, "grandavg_source_figures")

ROIS_OF_INTEREST = ["postcentral", "supramarginal"]
TARGET_LABEL_NAMES = [f"{roi}-{hemi}" for roi in ROIS_OF_INTEREST for hemi in ("lh", "rh")]

N2_WINDOW = (0.10, 0.30)
P2_WINDOW = (0.25, 0.45)

PLOT_TMIN = -0.2
PLOT_TMAX = 0.8

COLOR_LH = "#1f77b4"
COLOR_RH = "#aec7e8"

# =============================================================================
# HELPERS
# =============================================================================

def _find_peak(waveform: np.ndarray, t: np.ndarray,
               t_lo: float, t_hi: float, polarity: str):
    mask = (t >= t_lo) & (t <= t_hi)
    if not np.any(mask):
        return float("nan"), float("nan")
    w = waveform[mask]
    tw = t[mask]
    idx = int(np.argmin(w)) if polarity == "neg" else int(np.argmax(w))
    return float(w[idx]), float(tw[idx])


def _plot_roi_figure(roi: str, wave_lh: np.ndarray, wave_rh: np.ndarray,
                      times: np.ndarray, n_subjects: int, out_path: str):
    t_ms = times * 1000.0
    fig, ax = plt.subplots(figsize=(9, 5.5))

    for wave, color, label in [(wave_lh, COLOR_LH, f"{roi}-lh"),
                                 (wave_rh, COLOR_RH, f"{roi}-rh")]:
        ax.plot(t_ms, wave, label=label, color=color, lw=1.8)

        n2_amp, n2_lat = _find_peak(wave, times, N2_WINDOW[0], N2_WINDOW[1], "neg")
        p2_amp, p2_lat = _find_peak(wave, times, P2_WINDOW[0], P2_WINDOW[1], "pos")

        if not np.isnan(n2_amp):
            ax.plot(n2_lat * 1000, n2_amp, marker="v", color=color,
                     markersize=12, markeredgecolor="k", markeredgewidth=0.8, zorder=5)
        if not np.isnan(p2_amp):
            ax.plot(p2_lat * 1000, p2_amp, marker="^", color=color,
                     markersize=12, markeredgecolor="k", markeredgewidth=0.8, zorder=5)

    ax.axvspan(N2_WINDOW[0] * 1000, N2_WINDOW[1] * 1000,
               color="steelblue", alpha=0.08, label="N2 window")
    ax.axvspan(P2_WINDOW[0] * 1000, P2_WINDOW[1] * 1000,
               color="tomato", alpha=0.08, label="P2 window")

    ax.axvline(0, color="k", lw=1.0, linestyle="--")
    ax.axhline(0, color="grey", lw=0.6)

    ax.set_xlim(PLOT_TMIN * 1000, PLOT_TMAX * 1000)
    ax.set_xlabel("Time (ms)", fontsize=13)
    ax.set_ylabel("Signed source amplitude (a.u.)", fontsize=13)
    ax.set_title(f"Grand-Average N2-P2 Waveform (signed) — {roi}  (N={n_subjects})", fontsize=14)
    ax.tick_params(axis="both", labelsize=11)
    ax.legend(fontsize=11, loc="upper right", framealpha=0.85, markerscale=1.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[FIG] {out_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    checkpoint_paths = sorted(glob.glob(os.path.join(CHECKPOINT_DIR, "*_n2p2_signed.npz")))
    if not checkpoint_paths:
        raise RuntimeError(
            f"No checkpoints found in {CHECKPOINT_DIR}. "
            "Run driver_run_all_subjects.py first."
        )
    print(f"[LOAD] Found {len(checkpoint_paths)} checkpoints")

    arrays = {name: [] for name in TARGET_LABEL_NAMES}
    times_ref = None
    n_done = 0

    for path in checkpoint_paths:
        try:
            data = np.load(path)
            ga_tc = data["ga_tc"]     # (n_rois, n_times) — order matches TARGET_LABEL_NAMES
            times = data["times"]

            if times_ref is None:
                times_ref = times
            elif ga_tc.shape[1] != len(times_ref):
                print(f"[WARN] {os.path.basename(path)} — time axis mismatch, skipping")
                continue

            for i, name in enumerate(TARGET_LABEL_NAMES):
                arrays[name].append(ga_tc[i, :])

            n_done += 1

        except Exception as e:
            print(f"[WARN] {os.path.basename(path)} — could not load: {e}")

    print(f"[GA] N={n_done} subjects contributing")

    grand_avg = {
        name: np.mean(np.stack(vals, axis=0), axis=0)
        for name, vals in arrays.items()
    }

    mask = (times_ref >= PLOT_TMIN) & (times_ref <= PLOT_TMAX)
    times_plot = times_ref[mask]
    for name in grand_avg:
        grand_avg[name] = grand_avg[name][mask]

    for roi in ROIS_OF_INTEREST:
        wave_lh = grand_avg[f"{roi}-lh"]
        wave_rh = grand_avg[f"{roi}-rh"]
        out_path = os.path.join(OUT_DIR, f"grandavg_n2p2_erp_signed_{roi}.png")
        _plot_roi_figure(roi, wave_lh, wave_rh, times_plot, n_done, out_path)

    print("[DONE]")


if __name__ == "__main__":
    main()