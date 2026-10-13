"""
src_add_prefix.py - Prefix dict keys, skipping key/join columns.

Python counterpart to spec_add_prefix.m. Used to tag whole_/pre_/post_
feature dicts before merging into a single per-(trial, roi) row.
"""

from __future__ import annotations

def src_add_prefix(d: dict, prefix: str, skip_keys: tuple[str, ...] = ("trial", "roi_idx", "roi", "subject")) -> dict:
    """
    Return a new dict with every key except skip_keys renamed to
    f"{prefix}{key}". Key/join columns (trial, roi_idx, roi, subject) pass
    through unprefixed so rows can still be merged/joined on them afterwards.

    Args:
        d : source dict
        prefix : e.g. "whole_", "pre_", "post_"
        skip_keys : keys to leave unprefixed

    Returns:
        New dict, same values, renamed keys.
    """
    out = {}
    for k, v in d.items():
        out[k if k in skip_keys else f"{prefix}{k}"] = v
    return out


def src_add_prefix_rows(
        rows: list[dict],
        prefix: str,
        skip_keys: tuple[str, ...] = ("trial", "roi_idx", "roi", "subject"), 
) -> list[dict]:
    """
    Apply src_add_prefix to every dict in a list of rows.
    """
    return [src_add_prefix(r, prefix, skip_keys) for r in rows]