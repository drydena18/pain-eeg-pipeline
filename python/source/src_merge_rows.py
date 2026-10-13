"""
src_merge_rows.py - Merge multiple row-lists (list[dict]) sharing key
columns into one row per key.
"""

from __future__ import annotations

def src_merge_rows(
        *row_lists: list[dict],
        key_cols: tuple[str, ...] = ("trial", "roi_idx"),
) -> list[dict]:
    """
    Merge any number of row-lists sharing key_cols into one list of dicts,
    one per unique key combination.

    Fields from later row_lists overwrite same-named fields from earlier
    ones (including key_cols themselves). A row missing from a later list
    simply keeps whatever fields it already has from earlier lists.

    Args:
        *row_lists : any number of list[dict], each dict containing at
                     least the key_cols fields
        key_cols : tuple of dict keys to merge on

    Returns:
        List of merged dicts, one per unique key_cols combination, in
        first-seen order.
    """
    merged: dict[tuple, dict] = {}
    order: list[tuple] = []

    for rows in row_lists:
        for r in rows:
            key = tuple(r[k] for k in key_cols)
            if key not in merged:
                merged[key] = {}
                order.append(key)
            merged[key].update(r)

    return [merged[k] for k in order]