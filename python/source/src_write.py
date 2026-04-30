"""
src_write.py - CSV writers for the source localization pipeline

All writers follow the same pattern:
    - Accept pre-assembled list-of-dicts (rows)
    - Resolve roi_idx -> roi_name using a shared roi_names list
    - Prepend subject ID
    - Write via pandas, log path via logf

Output files produced
----------------------
Per-subject, per-trial:
    sub-XXX_source_trial.csv        pre+post spectral + LEP + phase metrics

Per-subject, grand-average:
    sub-XXX_source_ga.csv           GA pre+post spectral + LEP + TVI_alpha _ITC
    sub-XXX_source_ga_fooof.csv     GA FOOOF metrics (optional)
"""

from __future__ import annotations

import pandas as pd

from src_io import src_logmsg

# ====================================================================
# ROI NAME RESOLVER
# ====================================================================
def _add_roi_name(rows: list[dict], roi_names: list[str]) -> list[dict]:
    """
    Replace roi_idx with the string roi name in a copy of each row.
    """
    out = []
    for row in rows:
        r = dict(row)
        ri = r.pop("roi_idx", None)
        if ri is not None and ri < len(roi_names):
            r["roi"] = roi_names[ri]
        out.append(r)
    return out

# ====================================================================
# TRIAL CSV (pre-stim + post-stim + LEP merged)
# ====================================================================
def src_write_trial_csv(
        path: str,
        sub: int,
        roi_names: list[str],
        prestim_rows: list[dict],
        poststim_rows: list[dict],
        lep_rows: list[dict],
        logf,
):
    """
    Merge pre-stim, post-stim, and LEP rows into a single trial-level CSV

    All three lists have the same (trial, roi_idx) key pairs. They are merged
    on those keys so each output row represents one (subject, trial, ROI)

    Args:
        path            : Output CSV path
        sub             : Subject ID (integer)
        roi_names       : List mapping roi_idx -> roi_name
        prestim_rows    : List of dicts with pre-stim metrics, each with keys
                          'trial', 'roi_idx', and metric fields
        poststim_rows   : List of dicts with post-stim metrics, each with keys
                          'trial', 'roi_idx', and metric fields
        lep_rows        : List of dicts with LEP metrics, each with keys
                          'trial', 'roi_idx', and metric fields
        logf            : Open file handle for logging
    """
    def _to_df(rows, roi_names):
        named = _add_roi_name(rows, roi_names)
        return pd.DataFrame(named)
    
    df_pre = _to_df(prestim_rows, roi_names)
    df_post = _to_df(poststim_rows, roi_names)
    df_lep = _to_df(lep_rows, roi_names)

    # Merge on (trial, roi)
    df = df_pre.merge(df_post, on = ["trial", "roi"], how = "outer", suffixes = ("", "_post_dup"))
    df = df.merge(df_lep, on = ["trial", "roi"], how = "outer", suffixes = ("", "_lep_dup"))

    # Drop any accidential duplicate columns from suffix collisions
    dup_cols = [c for c in df.columns if c.endswith("_post_dup") or c.endswith("_lep_dup")]
    df.drop(columns = dup_cols, inplace = True)

    df.insert(0, "subject", sub)
    df.to_csv(path, index = False)
    src_logmsg(logf, "[CSV] %s (%d rows x %d cols)", path, len(df), len(df.columns))

# ====================================================================
# GRAND-AVERAGE CSV (pre + post + LEP GA + TVI_alpha + ITC)
# ====================================================================
def src_write_ga_csv(
        path: str,
        sub: int,
        roi_names: list[str],
        ga_prestim_rows: list[dict],
        ga_poststim_rows: list[dict],
        ga_lep_rows: list[dict],
        tvi_by_roi: dict,
        itc_rows: list[dict],
        logf,
):
    """
    Write the grand-average CSV combining pre-stim, post-stim, LEP, TVI_alpha,
    and ITC metrics, one row per ROI

    Args:
        path              : Output CSV file path.
        sub               : Integer subject ID.
        roi_names         : Ordered list of ROI name strings.
        ga_prestim_rows   : From src_compute_ga_prestim_metrics().
        ga_poststim_rows  : From src_compute_ga_poststim_metrics().
        ga_lep_rows       : From src_compute_lep_ga().
        tvi_by_roi        : Dict mapping roi_idx -> TVI_alpha scalar.
        itc_rows          : From src_compute_itc().
        logf              : Log file handle.
    """
    def _named(rows):
        return pd.DataFrame(_add_roi_name(rows, roi_names))
    
    df_pre = _named(ga_prestim_rows)
    df_post = _named(ga_poststim_rows)
    df_lep = _named(ga_lep_rows)
    df_itc = _named(itc_rows)

    # Add TVI_alpha into the pre-stim frame
    df_pre["TVI_alpha"] = df_pre["roi"].map(
        {roi_names[ri]: v for ri, v in tvi_by_roi.items()}
    )

    df = df_pre \
        .merge(df_post, on = "roi", how = "outer", suffixes = ("", "_post_dup")) \
        .merge(df_lep, on = "roi", how = "outer", suffixes = ("", "_lep_dup")) \
        .merge(df_itc, on = "roi", how = "outer", suffixes = ("", "_itc_dup"))
    
    dup_cols = [c for c in df.columns if c.endswith(("_post_dup", "_lep_dup", "_itc_dup"))]
    df.drop(columns = dup_cols, inplace = True)

    df.insert(0, "subjects", sub)
    df.to_csv(path, index = False)
    src_logmsg(logf, "[CSV] %s (%d rows x %d cols)", path, len(df), len(df.columns))


# ====================================================================
# FOOOF GA CSV
# ====================================================================
def src_write_fooof_csv(
        path: str,
        roi_names: list[str],
        fooof_rows: list[dict],
        logf,
):
    """
    Write the grand-average FOOOF metrics CSV

    Args:
        path        : Output CSV file path
        roi_names   : Ordered list of ROI name strings
        fooof_rows  : From src_compute_fooof_ga()
        logf        : Log file handle
    """
    df = pd.DataFrame(_add_roi_name(fooof_rows, roi_names))
    df.to_csv(path, index = False)
    src_logmsg(logf, "[CSV] %s (%d rows x %d cols)", path, len(df), len(df.columns))