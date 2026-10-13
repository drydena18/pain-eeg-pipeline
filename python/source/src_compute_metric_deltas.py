"""
src_compute_metric_deltas.py - delta_<name> = post_<name> - pre_<name>

Python counterpart to spec_compute_metric_deltas.m. Generic raw pre->post
change, computed for every field common to matched (by key_cols) pre_rows
and post_rows entries, for BOTH unprefixed inputs (caller passes the plain
src_compute_alpha_features output, before src_add_prefix is applied).

Works for both trial-level (key_cols=("trial", "roi_idx)) and GA-level
(key_cols=("roi_idx",)) row lists.
"""

from __future__ import annotations

import math

def src_compute_metric_deltas(
        pre_rows: list[dict],
        post_rows: list[dict],
        key_cols: tuple[str, ...] = ("trial", "roi_idx"),
) -> list[dict]:
    """
    Compute delta_<name> = post_<name> - pre_<name> for every field common
    to matched pre_rows/post_rows entries (matched by key_cols).

    Args:
        pre_rows, post_rows : list[dict], unprefixed
        key_cols : fields to match rows on

    Returns:
        List of dicts, one per matched (pre, post) pair: key_cols plus
        delta_<name> for every field common to both inpute (excluding
        key_cols). Rows present in only one input are skipped silently.
    """
    pre_by_key = {tuple(r[k] for k in key_cols): r for r in pre_rows}
    post_by_key = {tuple(r[k] for k in key_cols): r for r in post_rows}

    common_keys = [k for k in pre_by_key if k in post_by_key]

    delta_rows: list[dict] = []
    for key in common_keys:
        pre_r = pre_by_key[key]
        post_r = post_by_key[key]

        common_fields = (set(pre_r.keys()) & set(post_r.keys())) - set(key_cols)

        row = {k: pre_r[k] for k in key_cols}
        for fn in common_fields:
            pv, qv = pre_r[fn], post_r[fn]
            try:
                if isinstance(pv, (int, float)) and isinstance(qv, (int, float)):
                    row[f"delta_{fn}"] = float(qv) - float(pv)
                    if math.isnan(row[f"delta_{fn}"]):
                        row[f"delta_{fn}"] = float("nan")
            except (TypeError, ValueError):
                continue

        delta_rows.append(row)

    return delta_rows