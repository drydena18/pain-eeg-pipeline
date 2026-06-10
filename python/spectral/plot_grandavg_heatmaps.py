#!/usr/bin/env python3
"""
plot_grandavg_heatmaps.py
─────────────────────────────────────────────────────────────────────────────
Grand-average spectral heatmaps, one PNG per metric per experiment.

For each experiment found under --spectral-root, and for each of the 18
spectral metrics, produces one heatmap where:

    Y-axis  →  subjects (one row each), ordered ascending by global_subjid,
               labelled with subjid_uid  (e.g. E01_S003)
    X-axis  →  trial number
    Colour  →  channel-grand-averaged metric value for that (subject, trial)

Output layout:
    <out-dir>/
        26ByBiosemi/
            bi_pre.png
            cog_pre.png
            ...
        142ByBiosemi/
            bi_pre.png
            ...

Subjects in the behavioural master with no CSV for that experiment appear
as an all-NaN (light grey) row, preserving the global_subjid ordering gap.

Usage
─────
python plot_grandavg_heatmaps.py \\
    --spectral-root  /cifs/seminowicz/eegPainDatasets/CNED/da-analysis \\
    --behav-master   /path/to/behavioural_demo_master.csv \\
    --out-dir        /path/to/figures/grandavg \\
    [--colormap      viridis] \\
    [--clim-pct      2 98]    \\
    [--dpi           300]

Dependencies: numpy, pandas, matplotlib
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# METRIC DEFINITIONS
# (csv_column_name, display_title, cmap_hint)
# cmap_hint: 'seq' = sequential (--colormap arg), 'div' = diverging (RdBu_r)
# ─────────────────────────────────────────────────────────────────────────────
METRICS: list[tuple[str, str, str]] = [
    # Pre-stimulus interaction metrics
    ("bi_pre",          "BI_pre  (Balance Index)",             "div"),
    ("lr_pre",          "LR_pre  (Log-Ratio)",                 "div"),
    ("cog_pre",         "CoG_pre  (Hz)",                       "seq"),
    ("psi_cog",         "ψ_CoG  (BI_pre × [CoG_pre − 10])",   "div"),
    # Pre-stimulus phase
    ("phase_slow_rad",  "Slow α Phase at Onset  (rad)",        "div"),
    # Whole-epoch power
    ("pow_slow_alpha",  "Slow α Power  (8–10 Hz)",             "seq"),
    ("pow_fast_alpha",  "Fast α Power  (10–12 Hz)",            "seq"),
    ("pow_alpha_total", "Total α Power  (8–12 Hz)",            "seq"),
    # Whole-epoch relative / ratio
    ("rel_slow_alpha",  "Relative Slow α",                     "seq"),
    ("rel_fast_alpha",  "Relative Fast α",                     "seq"),
    ("sf_ratio",        "Slow/Fast Ratio",                     "div"),
    ("sf_logratio",     "Slow/Fast Log-Ratio  (whole epoch)",  "div"),
    ("sf_balance",      "Slow/Fast Balance  (whole epoch)",    "div"),
    ("slow_alpha_frac", "Slow α Fraction",                     "seq"),
    # Whole-epoch PAF
    ("paf_cog_hz",      "PAF CoG  (Hz)",                       "seq"),
    # Post-stimulus ERD
    ("erd_slow",        "ERD_slow  (8–10 Hz)",                 "div"),
    ("erd_fast",        "ERD_fast  (10–12 Hz)",                "div"),
    ("delta_erd",       "ΔERD  (ERD_slow − ERD_fast)",         "div"),
]

N_METRICS   = len(METRICS)
DIV_CMAP    = "RdBu_r"
NAN_COLOR   = "#D0D0D0"
FONT_FAMILY = ["sans-serif"]


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
    """Extract experiment name from path: segment immediately above 'preproc/'."""
    parts = csv_path.parts
    for i, part in enumerate(parts):
        if part.lower() == "preproc" and i > 0:
            return parts[i - 1]
    for i, part in enumerate(parts):
        if re.match(r"^sub-\d+$", part) and i >= 2:
            return parts[i - 2]
    return "unknown"


def load_behav_master(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    df.columns = normalize_colnames(df.columns.tolist())
    missing = {"subjid", "subjid_uid", "global_subjid"} - set(df.columns)
    if missing:
        raise ValueError(f"behavioural_demo_master missing columns: {missing}")
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
    Read one CSV, average across channels per trial.
    Returns tidy DataFrame with columns [subjid_uid, global_subjid, trial, <metrics>]
    or None on failure.
    """
    try:
        df = pd.read_csv(csv_path)
        df.columns = normalize_colnames(df.columns.tolist())
    except Exception as exc:
        print(f"[WARN] Could not read {csv_path.name}: {exc}", file=sys.stderr)
        return None

    if not {"subjid", "trial"}.issubset(df.columns):
        print(f"[WARN] {csv_path.name}: missing subjid/trial — skipping.", file=sys.stderr)
        return None

    unique_sids = df["subjid"].unique()
    if len(unique_sids) != 1:
        print(f"[WARN] {csv_path.name}: {len(unique_sids)} subjid values — skipping.", file=sys.stderr)
        return None
    this_sid = int(unique_sids[0])

    # Resolve to subjid_uid via experiment_name + subjid
    if exp_name and "experiment_name" in subj_df.columns:
        match = subj_df[
            (subj_df["experiment_name"] == exp_name) &
            (subj_df["subjid"] == this_sid)
        ]
    else:
        match = subj_df[subj_df["subjid"] == this_sid]

    if match.empty:
        print(f"[WARN] No master entry for subjid={this_sid}, exp={exp_name!r} — skipping.", file=sys.stderr)
        return None
    if len(match) > 1:
        match = match.iloc[[0]]

    subjid_uid    = match["subjid_uid"].iloc[0]
    global_subjid = int(match["global_subjid"].iloc[0])

    # Average across channels per trial for present metric columns only
    metric_cols_present = [m[0] for m in METRICS if m[0] in df.columns]
    for col in metric_cols_present:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    ga = df.groupby("trial")[metric_cols_present].mean(numeric_only=True).reset_index()

    # Fill any missing metric columns with NaN
    for col, _, _ in METRICS:
        if col not in ga.columns:
            ga[col] = np.nan

    ga["subjid_uid"]    = subjid_uid
    ga["global_subjid"] = global_subjid
    return ga


def build_heatmap_matrices(
    exp_data: pd.DataFrame,
    exp_subj_df: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Build [n_subj × n_trials × N_METRICS] array for one experiment.

    Only subjects who have at least one row in exp_data are included,
    ordered by global_subjid ascending.
    """
    trial_ids = np.sort(exp_data["trial"].unique()).astype(int)
    trial2col = {int(t): i for i, t in enumerate(trial_ids)}

    # Keep only subjects that actually have data, preserving global_subjid order
    uids_with_data = set(exp_data["subjid_uid"].unique())
    exp_subj_df = exp_subj_df[
        exp_subj_df["subjid_uid"].isin(uids_with_data)
    ].reset_index(drop=True)

    uid_list = exp_subj_df["subjid_uid"].tolist()
    uid2row  = {uid: i for i, uid in enumerate(uid_list)}

    matrices = np.full((len(uid_list), len(trial_ids), N_METRICS), np.nan)

    for _, row in exp_data.iterrows():
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
# FIGURE
# ─────────────────────────────────────────────────────────────────────────────

def render_heatmap(
    M:             np.ndarray,
    trial_ids:     np.ndarray,
    y_labels:      list[str],
    display_label: str,
    exp_name:      str,
    out_path:      Path,
    seq_cmap:      str   = "viridis",
    cmap_hint:     str   = "seq",
    clim_pct:      tuple | None = (2.0, 98.0),
    dpi:           int   = 300,
) -> None:
    """Render and save one heatmap PNG."""
    matplotlib.rcParams.update({
        "font.family":       FONT_FAMILY,
        "figure.facecolor":  "white",
        "axes.facecolor":    "white",
        "savefig.facecolor": "white",
        "text.color":        "#1A1A1A",
        "axes.labelcolor":   "#1A1A1A",
        "xtick.color":       "#444444",
        "ytick.color":       "#444444",
        "axes.edgecolor":    "#CCCCCC",
    })

    n_subj, n_trials = M.shape

    # Colormap
    cmap_obj = plt.get_cmap(DIV_CMAP if cmap_hint == "div" else seq_cmap).copy()
    cmap_obj.set_bad(NAN_COLOR)

    # Colour limits
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

    # For diverging maps, centre on zero when the data spans zero
    if cmap_hint == "div" and vmin < 0 < vmax:
        bound = max(abs(vmin), abs(vmax))
        vmin, vmax = -bound, bound

    # Figure size: aim for square-ish cells
    # Width: at least 6 in, scale with trial count; height: scale with subject count
    cell_w   = max(0.18, 8.0 / max(n_trials, 1))   # inches per cell
    cell_h   = max(0.28, 6.0 / max(n_subj,   1))
    fig_w    = max(7.0,  cell_w * n_trials + 2.0)   # +2 for colorbar + margins
    fig_h    = max(4.0,  cell_h * n_subj   + 1.5)   # +1.5 for title + xlabel

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")

    im = ax.imshow(
        np.ma.masked_invalid(M),
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

    ax.set_title(f"{exp_name}  —  {display_label}", fontsize=11,
                 fontweight="bold", pad=8, color="#1A1A1A")
    ax.set_xlabel("Trial", fontsize=9, labelpad=4)
    ax.set_ylabel("Subject  (subjid_uid)", fontsize=9, labelpad=4)

    # Y-axis: one label per subject, thin if many
    y_step = 1
    if n_subj > 20: y_step = 2
    if n_subj > 40: y_step = 4
    yt = np.arange(0, n_subj, y_step)
    ax.set_yticks(yt)
    ax.set_yticklabels([y_labels[i] for i in yt], fontsize=7,
                       fontfamily="monospace")

    # X-axis: actual trial numbers, thin if many
    x_step = 1
    if n_trials > 30:  x_step = 5
    if n_trials > 80:  x_step = 10
    if n_trials > 150: x_step = 20
    xt = np.arange(0, n_trials, x_step)
    ax.set_xticks(xt)
    ax.set_xticklabels([str(trial_ids[i]) for i in xt], fontsize=7,
                       rotation=45 if n_trials > 20 else 0,
                       ha="right" if n_trials > 20 else "center")

    for spine in ax.spines.values():
        spine.set_edgecolor("#CCCCCC")
        spine.set_linewidth(0.5)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Grand-average spectral heatmaps — one PNG per metric per experiment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--spectral-root", required=True,
        help="Root dir to search recursively for *_spectral_chan_by_trial.csv files.")
    p.add_argument("--behav-master", required=True,
        help="Path to behavioural_demo_master.csv.")
    p.add_argument("--out-dir", default=None,
        help="Output root. Defaults to <spectral-root>/../grandavg_heatmaps.")
    p.add_argument("--colormap", default="viridis",
        help="Sequential colormap name. Diverging metrics always use RdBu_r.")
    p.add_argument("--clim-pct", nargs=2, type=float, default=[2.0, 98.0],
        metavar=("LO", "HI"),
        help="Percentile clip for colour limits. Pass '0 100' for full range.")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def main() -> None:
    args   = parse_args()
    root   = Path(args.spectral_root)
    out_dir = Path(args.out_dir) if args.out_dir else root.parent / "grandavg_heatmaps"
    clim_pct = tuple(args.clim_pct) if args.clim_pct != [0.0, 100.0] else None

    # 1. Discover all CSVs, group by experiment name
    csv_paths = sorted(root.rglob("*_spectral_chan_by_trial.csv"))
    if not csv_paths:
        print(f"[ERROR] No CSVs found under {root}", file=sys.stderr)
        sys.exit(1)

    by_exp: dict[str, list[Path]] = defaultdict(list)
    for p in csv_paths:
        by_exp[resolve_experiment_name(p)].append(p)

    print(f"[GRANDAVG] Found {len(csv_paths)} CSV(s) across "
          f"{len(by_exp)} experiment(s): {sorted(by_exp)}")

    # 2. Load behavioural master
    print(f"[GRANDAVG] Loading: {args.behav_master}")
    subj_df = load_behav_master(Path(args.behav_master))

    # 3. Process each experiment independently
    for exp_name, exp_csvs in sorted(by_exp.items()):
        print(f"\n[GRANDAVG] ── Experiment: {exp_name}  ({len(exp_csvs)} subject(s)) ──")

        # Subject rows for this experiment only, sorted by global_subjid
        if "experiment_name" in subj_df.columns:
            exp_subj_df = (
                subj_df[subj_df["experiment_name"] == exp_name]
                .sort_values("global_subjid")
                .reset_index(drop=True)
            )
        else:
            exp_subj_df = subj_df.copy()

        if exp_subj_df.empty:
            print(f"[WARN] No master rows for experiment {exp_name!r} — skipping.")
            continue

        # Read + channel-average each subject CSV
        frames: list[pd.DataFrame] = []
        for csv_path in exp_csvs:
            ga = read_chan_trial_csv(csv_path, exp_name, subj_df)
            if ga is not None:
                frames.append(ga)
                n_metrics_present = sum(
                    1 for col, _, _ in METRICS if ga[col].notna().any()
                )
                print(f"  {ga['subjid_uid'].iloc[0]:15s}  "
                      f"trials={len(ga):3d}  "
                      f"metrics={n_metrics_present}/{N_METRICS}")

        if not frames:
            print(f"[WARN] No valid data for {exp_name} — skipping.")
            continue

        exp_data = pd.concat(frames, ignore_index=True)
        n_trials_exp = exp_data["trial"].nunique()
        print(f"  → {exp_data['subjid_uid'].nunique()} subjects, "
              f"{n_trials_exp} unique trial numbers")

        # Build matrices for this experiment
        matrices, trial_ids, y_labels = build_heatmap_matrices(exp_data, exp_subj_df)
        n_subj, n_trials, _ = matrices.shape
        print(f"  → Matrix: {n_subj} subjects × {n_trials} trials × {N_METRICS} metrics")

        # Render one PNG per metric
        exp_out = out_dir / exp_name
        n_written = 0
        for mi, (col, label, hint) in enumerate(METRICS):
            M = matrices[:, :, mi]
            if not np.isfinite(M).any():
                print(f"  SKIP  {col} (all NaN)")
                continue
            out_path = exp_out / f"{col}.png"
            render_heatmap(
                M             = M,
                trial_ids     = trial_ids,
                y_labels      = y_labels,
                display_label = label,
                exp_name      = exp_name,
                out_path      = out_path,
                seq_cmap      = args.colormap,
                cmap_hint     = hint,
                clim_pct      = clim_pct,
                dpi           = args.dpi,
            )
            n_written += 1

        print(f"  → {n_written} PNGs written to {exp_out}/")

    print("\n[GRANDAVG] Done.")


if __name__ == "__main__":
    main()