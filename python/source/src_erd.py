"""
src_erd.py - ERD family (fractional pre -> post change) + p5_flag QC gate.

Python counterpart to spec_compute_interaction_metrics.m V2.1.0. This
module now owns ONLY the fractional (ERD-style) pre -> post change, which is
NOT produced generically by src_compute_metric_deltas.py because it is not
meaningful for bounded/ratio metrics (sf_balance, sf_logratio, psi_cog...).

Behaviour matches MATLAB exactly; always computes a numeric ERD value using
an epsilon floor on the denominator (never NaN due to low pre-stim power),
and flags low-power trials separately via p5_flag rather than dropping the
value.
"""

from __future__ import annotations

import numpy as np

from src_spectral import src_psd_welch, src_bandpower

_EPS0 = 1e-12

# ==================================================================
# NOISE-FLOOR / 5TH-PERCENTILE THRESHOLDS (per ROI, pooled across trials)
# ==================================================================
def src_compute_noise_stats(
        tc_pre: np.ndarray,
        sfreq: float,
        slow: tuple[float, float],
        fast: tuple[float, float],
        fmin: float,
        fmax: float,
        psd_window_sec: float = 0.5,
) -> dict:
    """
    Per-ROI noise floor (eps0) and 5th-percentile pre-stim power thresholds
    for the ERD denominator guard / p5_flag, pooled across trials within
    each ROI (see module dosctring).

    To match MATLAB's literal pooling-across-channels convention instead,
    concatenate slow_vals/fast_vals/quiet_vals across ALL ROIs before
    taking percentiles/median, and use the same eps0/thresholds for every
    ROI.

    Args:
        tc_pre : (n_epochs, n_rois, n_times_pre)
        sfreq : Sampling frequency in Hz
        slow, fast : (lo, hi) band bounds in Hz
        fmin, fmax : PSD frequency range for Welch estimation
        psd_window_sec : Welch segment length (s)

    Returns:
        Dict roi_idx -> {"eps0": float, "thr_slow": float, "thr_fast": float}
    """
    n_epochs, n_rois, _ = tc_pre.shape
    stats: dict = {}

    for ri in range(n_rois):
        slow_vals: list[float] = []
        fast_vals: list[float] = []
        quiet_vals: list[float] = []

        for ei in range(n_epochs):
            freqs, psd = src_psd_welch(tc_pre[ei, ri, :], sfreq, fmin, fmax, psd_window_sec)
            slow_vals.append(src_bandpower(freqs, psd, slow[0], slow[1]))
            fast_vals.append(src_bandpower(freqs, psd, fast[0], fast[1]))
            idxQ = (freqs >= 45.0) & (freqs <= 55.0)
            if np.any(idxQ):
                quiet_vals.append(float(np.nanmedian(psd[idxQ])))

        slow_arr = np.asaray(slow_vals, dtype = float)
        fast_arr = np.asarray(fast_vals, dtype = float)

        thr_slow = float(np.nanpercentile(slow_arr, 5)) if np.any(~np.isnan(slow_arr)) else float("nan")
        thr_fast = float(np.nanpercentile(fast_arr, 5)) if np.any(~np.isnan(fast_arr)) else float("nan")

        if len(quiet_vals) > 0 and np.nanmedian(queit_vals) > 0:
            eps0 = float(np.nanmedian(quiet_vals))
        else:
            pooled = np.concatenate([slow_arr, fast_arr])
            pooled = pooled[~np.isnan(pooled)]
            eps0 = float(np.nanmedian(pooled)) * 1e-3 if pooled.size > 0 else _EPS0
        eps0 = max(eps0, _EPS0)

        stats[ri] = {"eps0": eps0, "thr_slow": thr_slow, "thr_fast": thr_fast}

    return stats


# ==================================================================
# PER-TRIAL ERD FAMILY
# ==================================================================
def src_compute_erd_metrics(
        pre_rows: list[dict],
        post_rows: list[dict],
        noise_stats: dict,
        key_cols: tuple[str, ...] = ("trial", "roi_idx"),
        use_p5_flag: bool = True,
) -> list[dict]:
    """
    Compute erd_slow, erd_fast, erd_pow_alpha_total, delta_erd, and
    (optionally) p5_flag from matched pre_/post_ rows.

    Reuses pow_slow_alpha/pow_fast_alpha/pow_alpha_total already computed
    by src_compute_window_alpha_features rather than recomputing PSDs a
    third time.

    Args:
        pre_rows, post_rows : unprefixed rows from src_alpha_features.py
                              (must contain pow_slow_alpha, pow_fast_alpha,
                              pow_alpha_total, and key_cols)
        noise_stats : from src_compute_noise_stats(), roi_idx -> stats
        key_cols : fields to match rows on
        use_p5_flag : if False, skips the percentile-threshold
                      flag entirely (used for the GA path, where
                      the percentile of a single value is
                      degenerate - matches
                      spec_compute_interaction_metrics.m's GA
                      behaviour of using only an epsilon guard)
    
    Returns:
        List of dicts: key_cols + erd_slow, erd_fast, erd_pow_alpha_total,
        delta_erd, p5_flag (p5_flag omitted if use_p5_flag = False)
    """
    pre_by_key = {tuple(r[k] for k in key_cols): r for r in pre_rows}

    rows: list[dict] = []
    for post_row in post_rows:
        key = tuple(post_row[k] for k in key_cols)
        pre_row = pre_by_key.get(key)
        if pre_row is None:
            continue

        ri = post_row["roi_idx"]
        st = noise_stats.get(ri, {"eps0": _EPS0, "thr_slow": float("nan"), "thr_fast": float("nan")})
        eps0 = st["eps0"]

        pow_pre_slow = pre_row["pow_slow_alpha"]
        pow_pre_fast = pre_row["pow_fast_alpha"]
        pow_pre_alpha = pre_row["pow_alpha_total"]
        pow_post_slow = post_row["pow_slow_alpha"]
        pow_post_fast = post_row["pow_fast_alpha"]
        pow_post_alpha = post_row["pow_alpha_total"]

        erd_slow = (pow_post_slow - pow_pre_slow) / (pow_pre_slow + eps0)
        erd_fast = (pow_post_fast - pow_pre_fast) / (pow_pre_fast + eps0)
        erd_pow_alpha_total = (pow_post_alpha - pow_pre_alpha) / (pow_pre_alpha + eps0)
        delta_erd = erd_slow - erd_fast

        row = {k: post_row[k] for k in key_cols}
        row.update({
            "erd_slow": float(erd_slow),
            "erd_fast": float(erd_fast),
            "erd_pow_alpha_total": float(erd_pow_alpha_total),
            "delta_erd": float(delta_erd),
        })

        if use_p5_flag:
            thr_slow, thr_fast = st["thr_slow"], st["thr_fast"]
            flagged = (
                (not np.isnan(pow_pre_slow) and not np.isnan(thr_slow) and pow_pre_slow < thr_slow)
                or (not np.isnan(pow_pre_fast) and not np.isnan(thr_fast) and pow_pre_fast < thr_fast)
                    )
            
            row["p5_flag"] = int(flagged)

        rows.append(row)

    return rows