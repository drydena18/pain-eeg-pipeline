#!/usr/bin/env python3
"""
plot_grandavg_heatmaps.py
─────────────────────────────────────────────────────────────────────────────
Grand-average spectral heatmaps across ALL experiments.

Produces one standalone PNG per metric (18 total), each showing:

    Y-axis  →  subjects, ordered ascending by global_subjid,
               labelled with subjid_uid  (e.g. E01_S003)
    X-axis  →  trial number
    Colour  →  channel-grand-averaged metric value for (subject, trial)

Subjects present in the behavioural master but with no spectral CSV appear
as an all-NaN (light grey) row, preserving the global_subjid ordering gap.

Pipeline context
────────────────
Reads:  *_spectral_chan_by_trial.csv   (one per subject)
        behavioural_demo_master.csv    (produced by merge_participants_into_behavioural.R)

Writes: <out-dir>/<metric_key>.png   (16 files)

Usage
─────
python plot_grandavg_heatmaps.py \\
    --spectral-root  /cifs/seminowicz/eegPainDatasets/CNED/da-analysis \\
    --behav-master   /path/to/behavioural_demo_master.csv \\
    --out-dir        /path/to/figures/grandavg \\
    [--colormap      viridis]   \\
    [--clim-pct      2 98]      \\
    [--dpi           300]

Dependencies: numpy, pandas, matplotlib  (all present in ~/env)
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# METRIC DEFINITIONS
# Each tuple: (csv_column_name, display_title, colormap_hint)
# colormap_hint: 'seq' = sequential (viridis), 'div' = diverging (RdBu_r)
# ─────────────────────────────────────────────────────────────────────────────


METRICS: list[tuple[str, str, str]] = [
    # ── Pre-stimulus interaction metrics ─────────────────────────────────────
    ("bi_pre",         "BI_pre  (Balance Index)",            "div"),
    ("lr_pre",         "LR_pre  (Log-Ratio)",                "div"),
    ("cog_pre",        "CoG_pre  (Hz)",                      "seq"),
    ("psi_cog",        "ψ_CoG  (BI_pre × [CoG_pre − 10])",  "div"),
    # ── Pre-stimulus phase ────────────────────────────────────────────────────
    ("phase_slow_rad", "Slow α Phase at Onset  (rad)",       "div"),
    # ── Whole-epoch: absolute power ──────────────────────────────────────────
    ("pow_slow_alpha", "Slow α Power  (8–10 Hz)",            "seq"),
    ("pow_fast_alpha", "Fast α Power  (10–12 Hz)",           "seq"),
    ("pow_alpha_total","Total α Power  (8–12 Hz)",           "seq"),
    # ── Whole-epoch: relative / ratio metrics ────────────────────────────────
    ("rel_slow_alpha", "Relative Slow α",                    "seq"),
    ("rel_fast_alpha", "Relative Fast α",                    "seq"),
    ("sf_ratio",       "Slow/Fast Ratio",                    "div"),
    ("sf_logratio",    "Slow/Fast Log-Ratio  (whole epoch)", "div"),
    ("sf_balance",     "Slow/Fast Balance  (whole epoch)",   "div"),
    ("slow_alpha_frac","Slow α Fraction",                    "seq"),
    # ── Whole-epoch: PAF ─────────────────────────────────────────────────────
    ("paf_cog_hz",     "PAF CoG  (Hz)",                      "seq"),
    # ── Post-stimulus: ERD ───────────────────────────────────────────────────
    ("erd_slow",       "ERD_slow  (8–10 Hz)",                "div"),
    ("erd_fast",       "ERD_fast  (10–12 Hz)",               "div"),
    ("delta_erd",      "ΔERD  (ERD_slow − ERD_fast)",        "div"),
]

N_METRICS = len(METRICS)

BG_COLOR  = "white"
NAN_COLOR = "#D0D0D0"
FONT_FAMILY = ["Calibri", "DejaVu Sans", "sans-serif"]

# Diverging colormaps for metrics that have a meaningful zero
DIV_CMAP = "RdBu_r"


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def normalize_colnames(columns: list[str]) -> list[str]:
    out = []
    for c in columns:
        c = re.sub(r"[^A-Za-z0-9_]", "_", c)
        c = re.sub(r"_+", "_", c)
        c = c.strip("_").lower()
        out.append(c)
    return out


def resolve_experiment_name(csv_path: Path) -> str:
    parts = csv_path.parts
    for i, part in enumerate(parts):
        if part.lower() == "preproc" and i > 0:
            return parts[i - 1]
    for i, part in enumerate(parts):
        if re.match(r"^sub-\d+$", part) and i >= 2:
            return parts[i - 2]
    return ""


def load_behav_master(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = normalize_colnames(df.columns.tolist())

    required = {"subjid", "subjid_uid", "global_subjid"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"behavioural_demo_master.csv missing columns: {missing}")

    df["subjid"]        = df["subjid"].astype(int)
    df["global_subjid"] = df["global_subjid"].astype(int)

    return (
        df[["experiment_name", "subjid", "subjid_uid", "global_subjid"]]
        .drop_duplicates()
        .sort_values("global_subjid")
        .reset_index(drop=True)
    )


def read_chan_trial_csv(
    csv_path: Path,
    exp_name: str,
    subj_df: pd.DataFrame,
) -> pd.DataFrame | None:
    """
    Read one *_spectral_chan_by_trial.csv, average across channels per trial,
    return a tidy DataFrame with columns:
        [subjid_uid, global_subjid, trial, <metric_col>, ...]

    Only metric columns that actually exist in the CSV are averaged —
    missing columns (e.g. pre/post metrics not yet computed) are silently
    skipped; they will remain NaN in the heatmap.
    """
    try:
        df = pd.read_csv(csv_path)
        df.columns = normalize_colnames(df.columns.tolist())
    except Exception as exc:
        print(f"[WARN] Could not read {csv_path.name}: {exc}", file=sys.stderr)
        return None

    if not {"subjid", "trial"}.issubset(df.columns):
        print(f"[WARN] {csv_path.name} missing subjid/trial — skipping.", file=sys.stderr)
        return None

    unique_subjids = df["subjid"].unique()
    if len(unique_subjids) != 1:
        print(f"[WARN] {csv_path.name} has {len(unique_subjids)} subjid values — skipping.", file=sys.stderr)
        return None
    this_subjid = int(unique_subjids[0])

    # Resolve subjid_uid
    if exp_name and "experiment_name" in subj_df.columns:
        match = subj_df[
            (subj_df["experiment_name"] == exp_name) &
            (subj_df["subjid"] == this_subjid)
        ]
    else:
        match = subj_df[subj_df["subjid"] == this_subjid]

    if match.empty:
        print(
            f"[WARN] No master entry for subjid={this_subjid}, exp={exp_name!r} — skipping.",
            file=sys.stderr,
        )
        return None
    if len(match) > 1:
        match = match.iloc[[0]]

    subjid_uid    = match["subjid_uid"].iloc[0]
    global_subjid = int(match["global_subjid"].iloc[0])

    # Average across channels — only for columns that exist in this CSV
    metric_cols_present = [m[0] for m in METRICS if m[0] in df.columns]
    if not metric_cols_present:
        print(f"[WARN] {csv_path.name} has no recognised metric columns — skipping.", file=sys.stderr)
        return None

    for col in metric_cols_present:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    ga = (
        df.groupby("trial")[metric_cols_present]
        .mean(numeric_only=True)
        .reset_index()
    )

    # Add NaN columns for any metrics not present in this CSV
    for col, _, _ in METRICS:
        if col not in ga.columns:
            ga[col] = np.nan

    ga["subjid_uid"]    = subjid_uid
    ga["global_subjid"] = global_subjid

    return ga


def build_heatmap_matrices(
    all_data: pd.DataFrame,
    subj_order: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Returns:
        matrices  : float [n_subj × n_trials × N_METRICS], NaN for gaps
        trial_ids : int   [n_trials]
        y_labels  : list  [n_subj]  subjid_uid strings, row 0 = bottom
    """
    trial_ids = np.sort(all_data["trial"].unique()).astype(int)
    trial2col = {int(t): i for i, t in enumerate(trial_ids)}

    uid_list = subj_order["subjid_uid"].tolist()
    uid2row  = {uid: i for i, uid in enumerate(uid_list)}

    matrices = np.full((len(uid_list), len(trial_ids), N_METRICS), np.nan)

    for _, row in all_data.iterrows():
        uid = row["subjid_uid"]
        if uid not in uid2row:
            continue
        r = uid2row[uid]
        c = trial2col.get(int(row["trial"]))
        if c is None:
            continue
        for mi, (col, _, _) in enumerate(METRICS):
            v = row.get(col, np.nan)
            matrices[r, c, mi] = float(v) if pd.notna(v) else np.nan

    return matrices, trial_ids, uid_list


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE — one standalone PNG per metric
# ─────────────────────────────────────────────────────────────────────────────

def render_single_heatmap(
    M:           np.ndarray,
    trial_ids:   np.ndarray,
    y_labels:    list[str],
    metric_key:  str,
    display_label: str,
    cmap_hint:   str,
    out_path:    Path,
    default_cmap: str = "viridis",
    clim_pct:    tuple[float, float] | None = (2.0, 98.0),
    dpi:         int = 300,
) -> None:
    """
    Render and save one heatmap PNG for a single metric.

    Parameters
    ──────────
    M             : [n_subj × n_trials] float array
    trial_ids     : int array of actual trial numbers (x-axis)
    y_labels      : subjid_uid list (y-axis, row 0 = bottom)
    metric_key    : snake_case column name (used for filename)
    display_label : human-readable title string
    cmap_hint     : 'seq' or 'div' — controls colormap choice
    out_path      : full Path to write the PNG
    default_cmap  : sequential colormap name (CLI --colormap)
    clim_pct      : (lo, hi) percentile clip; None = full range
    dpi           : export resolution
    """
    matplotlib.rcParams.update({
        "font.family":      FONT_FAMILY,
        "figure.facecolor": BG_COLOR,
        "axes.facecolor":   BG_COLOR,
        "text.color":       "#1A1A1A",
        "axes.labelcolor":  "#1A1A1A",
        "xtick.color":      "#444444",
        "ytick.color":      "#444444",
        "axes.edgecolor":   "#CCCCCC",
    })

    n_subj, n_trials = M.shape

    # ── Colormap ─────────────────────────────────────────────────────────────
    cmap_name = DIV_CMAP if cmap_hint == "div" else default_cmap
    cmap_obj  = plt.get_cmap(cmap_name).copy()
    cmap_obj.set_bad(NAN_COLOR)

    # ── Colour limits ─────────────────────────────────────────────────────────
    finite_vals = M[np.isfinite(M)]
    if len(finite_vals) == 0:
        vmin, vmax = 0.0, 1.0
    elif clim_pct is not None:
        vmin = float(np.percentile(finite_vals, clim_pct[0]))
        vmax = float(np.percentile(finite_vals, clim_pct[1]))
        if vmin == vmax:
            vmax = vmin + 1e-9
    else:
        vmin, vmax = float(finite_vals.min()), float(finite_vals.max())

    # ── Figure sizing ─────────────────────────────────────────────────────────
    # Width: trial axis — aim for ~8 px per trial, min 8 inches
    fig_w = max(8.0, n_trials * 8 / 100)
    # Height: subject axis — aim for 18 px per subject, min 5 inches
    fig_h = max(5.0, n_subj * 18 / 100 + 1.5)   # +1.5 for title + xlabel

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor=BG_COLOR)

    M_masked = np.ma.masked_invalid(M)
    im = ax.imshow(
        M_masked,
        aspect="auto",
        origin="lower",       # row 0 = bottom = lowest global_subjid
        cmap=cmap_obj,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
    )

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor("#CCCCCC")

    ax.set_title(display_label, fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel("Trial", fontsize=10, labelpad=5)
    ax.set_ylabel("Subject  (global_uid)", fontsize=10, labelpad=5)

    # Y-axis ticks
    y_step = 1
    if n_subj > 30: y_step = 2
    if n_subj > 60: y_step = 4
    yt = np.arange(0, n_subj, y_step)
    ax.set_yticks(yt)
    ax.set_yticklabels([y_labels[i] for i in yt], fontsize=8)

    # X-axis ticks
    x_step = 1
    if n_trials > 40:  x_step = 5
    if n_trials > 100: x_step = 10
    xt = np.arange(0, n_trials, x_step)
    ax.set_xticks(xt)
    ax.set_xticklabels(
        [str(trial_ids[i]) for i in xt],
        fontsize=8,
        rotation=45 if n_trials > 20 else 0,
        ha="right" if n_trials > 20 else "center",
    )

    for spine in ax.spines.values():
        spine.set_edgecolor("#CCCCCC")
        spine.set_linewidth(0.5)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor=BG_COLOR)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Grand-average spectral heatmaps — one PNG per metric.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--spectral-root", required=True,
        help="Root directory to search recursively for *_spectral_chan_by_trial.csv files.")
    p.add_argument("--behav-master", required=True,
        help="Path to behavioural_demo_master.csv.")
    p.add_argument("--out-dir", default=None,
        help="Output directory. Defaults to <spectral-root>/../grandavg_heatmaps.")
    p.add_argument("--colormap", default="viridis",
        help="Sequential colormap name (diverging metrics always use RdBu_r).")
    p.add_argument("--clim-pct", nargs=2, type=float, default=[2.0, 98.0],
        metavar=("LO", "HI"),
        help="Percentile clip for colour limits. Pass '0 100' for full range.")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    spectral_root = Path(args.spectral_root)
    behav_path    = Path(args.behav_master)
    out_dir = (
        Path(args.out_dir) if args.out_dir
        else spectral_root.parent / "grandavg_heatmaps"
    )
    clim_pct = tuple(args.clim_pct) if args.clim_pct != [0.0, 100.0] else None

    # ── 1. Discover CSVs ─────────────────────────────────────────────────────
    csv_paths = sorted(spectral_root.rglob("*_spectral_chan_by_trial.csv"))
    if not csv_paths:
        print(f"[ERROR] No CSVs found under {spectral_root}", file=sys.stderr)
        sys.exit(1)
    print(f"[GRANDAVG] Found {len(csv_paths)} subject CSV(s).")

    # ── 2. Load master ───────────────────────────────────────────────────────
    print(f"[GRANDAVG] Loading: {behav_path}")
    subj_df = load_behav_master(behav_path)
    print(f"[GRANDAVG] {len(subj_df)} subjects in master.")

    # ── 3. Read + channel-average ────────────────────────────────────────────
    frames: list[pd.DataFrame] = []
    for csv_path in csv_paths:
        exp_name = resolve_experiment_name(csv_path)
        ga = read_chan_trial_csv(csv_path, exp_name, subj_df)
        if ga is not None:
            frames.append(ga)
            n_metrics_present = sum(
                1 for col, _, _ in METRICS if ga[col].notna().any()
            )
            print(
                f"[GRANDAVG]   {ga['subjid_uid'].iloc[0]:12s}  "
                f"exp={exp_name:20s}  "
                f"trials={len(ga)}  "
                f"metrics_present={n_metrics_present}/{N_METRICS}"
            )

    if not frames:
        print("[ERROR] No valid data loaded.", file=sys.stderr)
        sys.exit(1)

    all_data = pd.concat(frames, ignore_index=True)
    print(
        f"[GRANDAVG] Loaded {all_data['subjid_uid'].nunique()} subjects, "
        f"{all_data['trial'].nunique()} unique trials."
    )

    # ── 4. Build matrices ─────────────────────────────────────────────────────
    matrices, trial_ids, y_labels = build_heatmap_matrices(all_data, subj_df)
    n_subj, n_trials, _ = matrices.shape
    print(f"[GRANDAVG] Matrix: {n_subj} subjects × {n_trials} trials × {N_METRICS} metrics.")

    # ── 5. Render one PNG per metric ──────────────────────────────────────────
    print(f"[GRANDAVG] Writing PNGs to: {out_dir}")
    for mi, (col, label, hint) in enumerate(METRICS):
        M = matrices[:, :, mi]
        n_finite = int(np.isfinite(M).sum())

        out_path = out_dir / f"{col}.png"

        if n_finite == 0:
            print(f"[GRANDAVG]   SKIP  {col:30s} (all NaN — column not yet computed in pipeline)")
            continue

        render_single_heatmap(
            M             = M,
            trial_ids     = trial_ids,
            y_labels      = y_labels,
            metric_key    = col,
            display_label = label,
            cmap_hint     = hint,
            out_path      = out_path,
            default_cmap  = args.colormap,
            clim_pct      = clim_pct,
            dpi           = args.dpi,
        )
        print(f"[GRANDAVG]   OK    {col:30s}  ({n_finite} finite cells) → {out_path.name}")

    print("[GRANDAVG] Done.")


if __name__ == "__main__":
    main()