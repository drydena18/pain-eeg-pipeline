"""
src_erd.py - ERD family (fractional pre -> post change) + p5_flag QC gate.
V 3.1.0

MATLAB twin: spec_noise_stats.m + spec_erd_metrics.m. This module owns ONLY
the fractional (ERD-style) pre -> post change, which is not produced
generically by src_compute_metric_deltas.py because it is not meaningful for
bounded/ratio metrics (sf_balance, sf_logratio, psi_cog, ...).

ERD is always numeric: the denominator carries a floor eps0, and low-power
trials are flagged via p5_flag rather than dropped.

Per-unit statistics (identical in both pipelines; unit = ROI here, channel
in MATLAB), pooled across trials within the unit:
    eps0     = eps0_frac * median(pre-stim slow and fast power), >= 1e-12
    thr_slow = p-th percentile of pre-stim slow power  (numpy 'linear')
    thr_fast = p-th percentile of pre-stim fast power
    p5_flag  = pre slow < thr_slow  OR  pre fast < thr_fast

V3.1.0: delta_erd renamed erd_asym (= erd_slow - erd_fast). The delta_
prefix is reserved for post - pre change scores.

V3.0.0 changes vs V2.0.0:
    - REMOVED the 45-55 Hz quiet-band noise floor. The preprocessing
      low-pass is 40 Hz, so that band only measured filter stopband residue.
      eps0 is now a data-relative guard (eps0_frac x median pre-stim power).
    - Percentile is configurable (p5_percentile, default 5).
"""

from __future__ import annotations

import numpy as np

_EPS0 = 1e-12


# ==================================================================
# PER-UNIT NOISE FLOOR / PERCENTILE THRESHOLDS
# ==================================================================
def src_compute_noise_stats(
        pre_rows: list[dict],
        eps0_frac: float = 1e-3,
        p_pct: float = 5.0,
) -> dict:
    """
    Per-ROI eps0 and percentile thresholds from the pre-window trial rows.

    Args:
        pre_rows  : unprefixed pre-window rows (roi_idx, pow_slow_alpha,
                    pow_fast_alpha)
        eps0_frac : eps0 as a fraction of the median pre-stim power
        p_pct     : percentile for the p5_flag thresholds

    Returns:
        Dict roi_idx -> {"eps0", "thr_slow", "thr_fast"}
    """
    roi_ids = sorted({r["roi_idx"] for r in pre_rows})
    stats: dict = {}

    for ri in roi_ids:
        slow_arr = np.asarray([r["pow_slow_alpha"] for r in pre_rows if r["roi_idx"] == ri], dtype = float)
        fast_arr = np.asarray([r["pow_fast_alpha"] for r in pre_rows if r["roi_idx"] == ri], dtype = float)

        thr_slow = float(np.nanpercentile(slow_arr, p_pct)) if np.any(~np.isnan(slow_arr)) else float("nan")
        thr_fast = float(np.nanpercentile(fast_arr, p_pct)) if np.any(~np.isnan(fast_arr)) else float("nan")

        pooled = np.concatenate([slow_arr, fast_arr])
        pooled = pooled[~np.isnan(pooled)]
        eps0 = float(np.median(pooled)) * eps0_frac if pooled.size > 0 else _EPS0
        eps0 = max(eps0, _EPS0)

        stats[ri] = {"eps0": eps0, "thr_slow": thr_slow, "thr_fast": thr_fast}

    return stats


# ==================================================================
# ERD FAMILY
# ==================================================================
def src_compute_erd_metrics(
        pre_rows: list[dict],
        post_rows: list[dict],
        noise_stats: dict,
        key_cols: tuple[str, ...] = ("trial", "roi_idx"),
        use_p5_flag: bool = True,
) -> list[dict]:
    """
    erd_slow, erd_fast, erd_pow_alpha_total, erd_asym and (optionally)
    p5_flag from matched unprefixed pre/post rows.

    Args:
        pre_rows, post_rows : rows with pow_slow_alpha, pow_fast_alpha,
                              pow_alpha_total and key_cols
        noise_stats         : from src_compute_noise_stats()
        key_cols            : fields to match rows on
        use_p5_flag         : False for the subject-level (GA) rows, where a
                              percentile of one value is degenerate

    Returns:
        List of dicts: key_cols + erd_slow, erd_fast, erd_pow_alpha_total,
        erd_asym [, p5_flag]
    """
    pre_by_key = {tuple(r[k] for k in key_cols): r for r in pre_rows}

    rows: list[dict] = []
    for post_row in post_rows:
        key = tuple(post_row[k] for k in key_cols)
        pre_row = pre_by_key.get(key)
        if pre_row is None:
            continue

        st = noise_stats.get(post_row["roi_idx"],
                             {"eps0": _EPS0, "thr_slow": float("nan"), "thr_fast": float("nan")})
        eps0 = st["eps0"]

        ps0, pf0, pa0 = pre_row["pow_slow_alpha"], pre_row["pow_fast_alpha"], pre_row["pow_alpha_total"]
        ps1, pf1, pa1 = post_row["pow_slow_alpha"], post_row["pow_fast_alpha"], post_row["pow_alpha_total"]

        erd_slow = (ps1 - ps0) / (ps0 + eps0)
        erd_fast = (pf1 - pf0) / (pf0 + eps0)
        erd_alpha = (pa1 - pa0) / (pa0 + eps0)

        row = {k: post_row[k] for k in key_cols}
        row.update({
            "erd_slow": float(erd_slow),
            "erd_fast": float(erd_fast),
            "erd_pow_alpha_total": float(erd_alpha),
            "erd_asym": float(erd_slow - erd_fast),
        })

        if use_p5_flag:
            thr_slow, thr_fast = st["thr_slow"], st["thr_fast"]
            flagged = (
                (not np.isnan(ps0) and not np.isnan(thr_slow) and ps0 < thr_slow)
                or (not np.isnan(pf0) and not np.isnan(thr_fast) and pf0 < thr_fast)
            )
            row["p5_flag"] = int(flagged)

        rows.append(row)

    return rows
