"""
src_write.py - CSV writers for the source localization pipeline.
V 2.1.0

All writers:
    - accept pre-merged list-of-dicts rows (merged in source_core.py with
      src_merge_rows.py, so a field collision is visible there and never
      silently resolved here)
    - resolve roi_idx -> roi name
    - put 'subject' first

Output files (per subject)
--------------------------
    sub-XXX_source_trial.csv     trial x ROI: whole/pre/post/delta, ERD, p5_flag,
                                 phase, LEP                 (read by R merge)
    sub-XXX_source_ga.csv        subject x ROI: same families from trial-mean
                                 powers/PSD/waveform + TVI_alpha + ITC + FOOOF
    sub-XXX_source_ga_fooof.csv  subject x ROI FOOOF only (convenience copy)

V2.1.0: the FOOOF CSV now carries 'subject' like the other two files
(src_fooof V3.0.0 no longer embeds it in each row).
"""

from __future__ import annotations

import pandas as pd

from src_io import src_logmsg


def _add_roi_name(rows: list[dict], roi_names: list[str]) -> list[dict]:
    """Replace roi_idx with the ROI name in a copy of each row."""
    out = []
    for row in rows:
        r = dict(row)
        ri = r.pop("roi_idx", None)
        if ri is not None and ri < len(roi_names):
            r["roi"] = roi_names[ri]
        out.append(r)
    return out


def _write(path: str, sub: int, roi_names: list[str], rows: list[dict], logf):
    df = pd.DataFrame(_add_roi_name(rows, roi_names))
    if "subject" in df.columns:
        del df["subject"]
    df.insert(0, "subject", sub)
    if "roi" in df.columns:
        df.insert(1, "roi", df.pop("roi"))
    df.to_csv(path, index = False)
    src_logmsg(logf, "[CSV] %s (%d rows x %d cols)", path, len(df), len(df.columns))


def src_write_trial_csv(path: str, sub: int, roi_names: list[str], trial_rows: list[dict], logf):
    """Per-subject trial x ROI CSV."""
    _write(path, sub, roi_names, trial_rows, logf)


def src_write_ga_csv(path: str, sub: int, roi_names: list[str], ga_rows: list[dict], logf):
    """Per-subject subject x ROI (grand-average) CSV."""
    _write(path, sub, roi_names, ga_rows, logf)


def src_write_fooof_csv(path: str, sub: int, roi_names: list[str], fooof_rows: list[dict], logf):
    """Per-subject subject x ROI FOOOF CSV."""
    _write(path, sub, roi_names, fooof_rows, logf)
