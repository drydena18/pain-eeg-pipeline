"""
src_write.py - CSV writers for the source localization pipeline.
V 2.0.0

V2.0.0 changes vs V1.x:
    - src_write_trial_csv and src_write_ga_csv now accept a SINGLE
      already-merged list of rows (built with src_merge_rows.py in
      source_core.py) instead of several separate row-lists merged here
      via pandas .merge(..., suffixes=(...)). The old pattern silently
      dropped a column if two row-lists happened to define the same
      non-key field name (the "_dup" columns were dropped without a
      warning) — exactly the kind of silent-failure mode this codebase
      tries to avoid elsewhere. With merging now done explicitly via
      src_merge_rows.py before these functions are called, a field
      collision is visible in the merge step itself rather than silently
      resolved by a suffix-matching rule here.
    - TVI_alpha is no longer special-cased into the writer; it's expected
      to already be a column on the GA rows passed in (added via
      src_merge_rows in source_core.py), same as every other metric.

All writers follow the same pattern:
    - Accept pre-assembled, pre-merged list-of-dicts (rows)
    - Resolve roi_idx -> roi_name using a shared roi_names list
    - Prepend subject ID
    - Write via pandas, log path via logf

Output files produced
----------------------
Per-subject, per-trial:
    sub-XXX_source_trial.csv        whole+pre+post+delta+ERD+LEP+phase metrics

Per-subject, grand-average:
    sub-XXX_source_ga.csv           GA whole+pre+post+delta+ERD+LEP+TVI_alpha+ITC

Per-subject, FOOOF (optional):
    sub-XXX_source_ga_fooof.csv     GA FOOOF metrics
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
# TRIAL CSV (already merged: whole + pre + post + delta + ERD + LEP + phase)
# ====================================================================
def src_write_trial_csv(
        path: str,
        sub: int,
        roi_names: list[str],
        trial_rows: list[dict],
        logf,
):
    """
    Write the per-subject trial-level CSV from an already-merged row list.

    Args:
        path        : Output CSV path
        sub         : Subject ID (integer)
        roi_names   : List mapping roi_idx -> roi_name
        trial_rows  : List of dicts, one per (trial, roi), already merged
                      via src_merge_rows.py (whole_/pre_/post_/delta_/erd_
                      metrics, LEP, phase, p5_flag, etc.)
        logf        : Open file handle for logging
    """
    df = pd.DataFrame(_add_roi_name(trial_rows, roi_names))
    df.insert(0, "subject", sub)
    df.to_csv(path, index=False)
    src_logmsg(logf, "[CSV] %s (%d rows x %d cols)", path, len(df), len(df.columns))

# ====================================================================
# GRAND-AVERAGE CSV (already merged: whole + pre + post + delta + ERD +
# LEP + TVI_alpha + ITC)
# ====================================================================
def src_write_ga_csv(
        path: str,
        sub: int,
        roi_names: list[str],
        ga_rows: list[dict],
        logf,
):
    """
    Write the grand-average CSV from an already-merged row list, one row
    per ROI.

    Args:
        path      : Output CSV file path.
        sub       : Integer subject ID.
        roi_names : Ordered list of ROI name strings.
        ga_rows   : List of dicts, one per ROI, already merged via
                    src_merge_rows.py (whole_/pre_/post_/delta_/erd_
                    metrics, LEP, TVI_alpha, ITC, etc.)
        logf      : Log file handle.
    """
    df = pd.DataFrame(_add_roi_name(ga_rows, roi_names))
    df.insert(0, "subject", sub)
    df.to_csv(path, index=False)
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
    Write the grand-average FOOOF metrics CSV. Unchanged from V1.x.

    Args:
        path        : Output CSV file path
        roi_names   : Ordered list of ROI name strings
        fooof_rows  : From src_compute_fooof_ga()
        logf        : Log file handle
    """
    df = pd.DataFrame(_add_roi_name(fooof_rows, roi_names))
    df.to_csv(path, index=False)
    src_logmsg(logf, "[CSV] %s (%d rows x %d cols)", path, len(df), len(df.columns))